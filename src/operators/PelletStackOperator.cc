
#include "operators/PelletStackOperator.h"
#include "operators/map/NodeToNodeMap.h"
#include "vectors/VectorSelector.h"

namespace AMP {
  namespace Operator {

    PelletStackOperator :: PelletStackOperator(const boost::shared_ptr<PelletStackOperatorParameters> & params)
      : Operator(params) {
        d_totalNumberOfPellets = (params->d_db)->getInteger("TOTAL_NUMBER_OF_PELLETS");
        d_useSerial = (params->d_db)->getBoolWithDefault("USE_SERIAL", false);
        d_onlyZcorrection = (params->d_db)->getBoolWithDefault("ONLY_Z_CORRECTION", false);
        AMP_ASSERT(!(d_useSerial && d_onlyZcorrection));
        d_masterId = (params->d_db)->getInteger("MASTER");
        d_slaveId = (params->d_db)->getInteger("SLAVE");
        if((params->d_db)->keyExists("SCALING_FACTOR")) {
          d_useScaling = true;
          d_scalingFactor = (params->d_db)->getDouble("SCALING_FACTOR");
        } else {
          d_useScaling = false;
        }
        d_frozenVectorSet = false;
        std::string varName = (params->d_db)->getString("Variable");
        d_var.reset(new AMP::LinearAlgebra::Variable(varName));       
        std::string meshNamePrefix = (params->d_db)->getString("MeshNamePrefix");
        d_currentPellet = params->d_currentPellet;
        d_pelletStackComm = params->d_pelletStackComm;
        d_n2nMaps = params->d_n2nMaps;
        for(unsigned int pellId = 0; pellId < d_totalNumberOfPellets; pellId++) {
          char pellId2Str[256];
          sprintf(pellId2Str, "%u", (pellId + 1));
          std::string meshName = meshNamePrefix + "_" + pellId2Str;
          AMP::Mesh::Mesh::shared_ptr currMesh = d_Mesh->Subset(meshName);
          if(currMesh == NULL) {
            continue;
          }
          currMesh->setName(meshName);      // This is needed since subset may change the name and we rely on it later
          d_pelletIds.push_back(pellId);
          d_meshes.push_back(currMesh);
        }//end for pellId
      }

    void PelletStackOperator :: reset(const boost::shared_ptr<OperatorParameters>& params) {
      boost::shared_ptr<PelletStackOperatorParameters> myParams = boost::dynamic_pointer_cast<PelletStackOperatorParameters>(params);
      d_currentPellet = myParams->d_currentPellet;
    }

    std::vector<AMP::Mesh::Mesh::shared_ptr> PelletStackOperator :: getLocalMeshes() {
      return d_meshes;
    }

    std::vector<unsigned int> PelletStackOperator :: getLocalPelletIds() {
      return d_pelletIds;
    }

    bool PelletStackOperator :: useSerial() {
      return d_useSerial;
    }

    bool PelletStackOperator :: onlyZcorrection() {
      return d_onlyZcorrection;
    }

    bool PelletStackOperator :: useScaling() {
      return d_useScaling;
    }

    unsigned int PelletStackOperator :: getTotalNumberOfPellets() {
      return d_totalNumberOfPellets;
    }

    int PelletStackOperator :: getLocalIndexForPellet(unsigned int pellId) {
      for(size_t i = 0; i < d_pelletIds.size(); i++) {
        if(d_pelletIds[i] == pellId) {
          return i;
        }
      }//end for i
      return -1;
    }

    void PelletStackOperator :: applyUnscaling(AMP::LinearAlgebra::Vector::shared_ptr f) {
      AMP::LinearAlgebra::Vector::shared_ptr subF = f->subsetVectorForVariable(d_var);
      AMP::Discretization::DOFManager::shared_ptr dof_map = subF->getDOFManager();
      AMP::Mesh::MeshIterator bnd = d_Mesh->getBoundaryIDIterator(AMP::Mesh::Vertex, d_slaveId, 0);
      AMP::Mesh::MeshIterator end_bnd = bnd.end();
      for( ; bnd != end_bnd; ++bnd) {
        std::vector<size_t> bndGlobalIds;
        dof_map->getDOFs(bnd->globalID(), bndGlobalIds);
        for(size_t j = 0; j < bndGlobalIds.size(); ++j) {
          double val = subF->getLocalValueByGlobalID( bndGlobalIds[j] );
          subF->setLocalValueByGlobalID(bndGlobalIds[j], val/d_scalingFactor);
        }//end for j
      }//end for bnd
    }

    void PelletStackOperator :: apply(const AMP::LinearAlgebra::Vector::shared_ptr &f,
        const AMP::LinearAlgebra::Vector::shared_ptr &u, AMP::LinearAlgebra::Vector::shared_ptr &r,
        const double a, const double b) {
      if(d_onlyZcorrection) {
        applyOnlyZcorrection(r);
      } else {
        if(!d_frozenVectorSet) {
          d_frozenVectorForMaps = d_n2nMaps->getFrozenVector();
          d_frozenVectorSet = true;
        }
        if(d_useSerial) {
          applySerial(f, u, r);
        } else {
          applyXYZcorrection(f, u, r);
        }
      }
    }

    void PelletStackOperator :: applyOnlyZcorrection(AMP::LinearAlgebra::Vector::shared_ptr &u) {
      std::vector<double> finalMaxZdispsList;
      computeZscan(u, finalMaxZdispsList); 
      AMP::LinearAlgebra::Vector::shared_ptr subU = u->subsetVectorForVariable(d_var);
      AMP::Discretization::DOFManager::shared_ptr dof_map = subU->getDOFManager();
      for(size_t i = 0; i < d_pelletIds.size(); ++i) {
        if(d_pelletIds[i] > 0) {
          AMP::Mesh::MeshIterator nd = d_meshes[i]->getIterator(AMP::Mesh::Vertex, 0);
          AMP::Mesh::MeshIterator end_nd = nd.end();
          for(; nd != end_nd; ++nd) {
            std::vector<size_t> dofIds;
            dof_map->getDOFs(nd->globalID(), dofIds);       
            subU->addLocalValueByGlobalID(dofIds[2], finalMaxZdispsList[d_pelletIds[i] - 1]);
          }//end for nd
        }
      }//end for i
    }

    void PelletStackOperator :: applyXYZcorrection(const AMP::LinearAlgebra::Vector::shared_ptr &f,
        const AMP::LinearAlgebra::Vector::shared_ptr &u, AMP::LinearAlgebra::Vector::shared_ptr  &r) {
      AMP_ASSERT(d_frozenVectorSet);
      AMP::LinearAlgebra::Vector::shared_ptr nullVec;
      r->copyVector(f);
      d_frozenVectorForMaps->zero();
      d_n2nMaps->apply(nullVec, u, nullVec, 1.0, 0.0);
      AMP::LinearAlgebra::Vector::shared_ptr subU = d_frozenVectorForMaps->subsetVectorForVariable(d_var);
      AMP::LinearAlgebra::Vector::shared_ptr subR = r->subsetVectorForVariable(d_var);
      AMP::Discretization::DOFManager::shared_ptr dof_map = subR->getDOFManager();
      AMP::Mesh::MeshIterator bnd = d_Mesh->getBoundaryIDIterator(AMP::Mesh::Vertex, d_slaveId, 0);
      AMP::Mesh::MeshIterator end_bnd = bnd.end();
      for( ; bnd != end_bnd; ++bnd) {
        std::vector<size_t> bndGlobalIds;
        dof_map->getDOFs(bnd->globalID(), bndGlobalIds);
        for(unsigned int j = 0; j < bndGlobalIds.size(); j++) {
          double val = subU->getLocalValueByGlobalID( bndGlobalIds[j] );
          subR->addLocalValueByGlobalID(bndGlobalIds[j], val);
        }//end for j
      }//end for bnd
      std::vector<double> finalMaxZdispsList;
      computeZscan(u, finalMaxZdispsList); 
      for(size_t i = 0; i < d_pelletIds.size(); ++i) {
        if(d_pelletIds[i] > 1) {
          AMP::Mesh::MeshIterator bnd = d_meshes[i]->getBoundaryIDIterator(AMP::Mesh::Vertex, d_slaveId, 0);
          AMP::Mesh::MeshIterator end_bnd = bnd.end();
          for( ; bnd != end_bnd; ++bnd) {
            std::vector<size_t> bndGlobalIds;
            dof_map->getDOFs(bnd->globalID(), bndGlobalIds);
            subR->addLocalValueByGlobalID(bndGlobalIds[2], finalMaxZdispsList[d_pelletIds[i] - 2]);
          }//end for bnd
        }
      }//end for i
    }

    void PelletStackOperator :: applySerial(const AMP::LinearAlgebra::Vector::shared_ptr &f,
        const AMP::LinearAlgebra::Vector::shared_ptr &u, AMP::LinearAlgebra::Vector::shared_ptr &r) {
      AMP_ASSERT(d_frozenVectorSet);
      AMP_ASSERT(d_currentPellet > 0);
      AMP::LinearAlgebra::Vector::shared_ptr nullVec;
      int currPellIdx = getLocalIndexForPellet(d_currentPellet);
      int prevPellIdx = getLocalIndexForPellet(d_currentPellet - 1);
      size_t numMaps = d_n2nMaps->getNumberOfOperators();
      if(currPellIdx != -1) {
        for(size_t m = 0; m < numMaps; m++) {
          boost::shared_ptr<AMP::Operator::NodeToNodeMap> currMap = boost::dynamic_pointer_cast<
            AMP::Operator::NodeToNodeMap>(d_n2nMaps->getOperator(m));
          if(currMap->getMesh(2) == d_meshes[currPellIdx]) {
            currMap->applyStart(nullVec, u, nullVec, 1.0, 0.0);
            break;
          }
        }//end for m
      }
      if(prevPellIdx != -1) {
        for(size_t m = 0; m < numMaps; m++) {
          boost::shared_ptr<AMP::Operator::NodeToNodeMap> currMap = boost::dynamic_pointer_cast<
            AMP::Operator::NodeToNodeMap>(d_n2nMaps->getOperator(m));
          if(currMap->getMesh(1) == d_meshes[prevPellIdx]) {
            currMap->applyStart(nullVec, u, nullVec, 1.0, 0.0);
            break;
          }
        }//end for m
      }
      if(currPellIdx != -1) {
        for(size_t m = 0; m < numMaps; m++) {
          boost::shared_ptr<AMP::Operator::NodeToNodeMap> currMap = boost::dynamic_pointer_cast<
            AMP::Operator::NodeToNodeMap>(d_n2nMaps->getOperator(m));
          if((currMap->getMesh(2)) == d_meshes[currPellIdx]) {
            currMap->applyFinish(nullVec, u, nullVec, 1.0, 0.0);
            break;
          }
        }//end for m
      }
      if(prevPellIdx != -1) {
        for(size_t m = 0; m < numMaps; m++) {
          boost::shared_ptr<AMP::Operator::NodeToNodeMap> currMap = boost::dynamic_pointer_cast<
            AMP::Operator::NodeToNodeMap>(d_n2nMaps->getOperator(m));
          if((currMap->getMesh(1)) == d_meshes[prevPellIdx]) {
            currMap->applyFinish(nullVec, u, nullVec, 1.0, 0.0);
            break;
          }
        }//end for m
      }
      AMP::LinearAlgebra::Vector::shared_ptr subF = f->subsetVectorForVariable(d_var);
      AMP::LinearAlgebra::Vector::shared_ptr subR = r->subsetVectorForVariable(d_var);
      AMP::LinearAlgebra::Vector::shared_ptr subU = d_frozenVectorForMaps->subsetVectorForVariable(d_var);
      AMP::Discretization::DOFManager::shared_ptr dof_map = subR->getDOFManager();
      if(currPellIdx != -1) {
        subR->copyVector(subF);
        AMP::Mesh::MeshIterator bnd = d_meshes[currPellIdx]->getBoundaryIDIterator(AMP::Mesh::Vertex, d_slaveId, 0);
        AMP::Mesh::MeshIterator end_bnd = bnd.end();
        for( ; bnd != end_bnd; ++bnd) {
          std::vector<size_t> bndGlobalIds;
          dof_map->getDOFs(bnd->globalID(), bndGlobalIds);
          for(size_t j = 0; j < bndGlobalIds.size(); ++j) {
            double val = subU->getLocalValueByGlobalID( bndGlobalIds[j] );
            subR->addLocalValueByGlobalID(bndGlobalIds[j], val);
          }//end for j
        }//end for bnd
      }
    }

    void PelletStackOperator :: computeZscan(const AMP::LinearAlgebra::Vector::shared_ptr &u, 
        std::vector<double> &finalMaxZdispsList) {
      AMP::LinearAlgebra::Vector::shared_ptr subU = u->subsetVectorForVariable(d_var);
      AMP::Discretization::DOFManager::shared_ptr dof_map = subU->getDOFManager();
      std::vector<double> myMaxZdisps(d_pelletIds.size(), 0.0);
      for(size_t i = 0; i < d_pelletIds.size(); i++) {
        AMP::Mesh::MeshIterator bnd = d_meshes[i]->getBoundaryIDIterator(AMP::Mesh::Vertex, d_masterId, 0);
        AMP::Mesh::MeshIterator end_bnd = bnd.end();
        for( ; bnd != end_bnd; ++bnd) {
          std::vector<size_t> bndGlobalIds;
          dof_map->getDOFs(bnd->globalID(), bndGlobalIds);
          double val = subU->getLocalValueByGlobalID( bndGlobalIds[2] );
          if(fabs(myMaxZdisps[i]) < fabs(val)) {
            myMaxZdisps[i] = val;
          }
        }//end for bnd
      }//end for i

      std::vector<int> recvCnts(d_pelletStackComm.getSize()); 
      d_pelletStackComm.allGather<int>(d_pelletIds.size(), &(recvCnts[0]));

      std::vector<int> recvDisps(recvCnts.size());
      recvDisps[0] = 0;
      for(size_t i = 1; i < recvDisps.size(); i++) {
        recvDisps[i] = recvDisps[i - 1] + recvCnts[i - 1];
      }//end for i

      std::vector<unsigned int> allPelletIds((*(recvDisps.end() - 1)) + (*(recvCnts.end() - 1)));
      d_pelletStackComm.allGather<unsigned int>(&(d_pelletIds[0]), d_pelletIds.size(), &(allPelletIds[0]), 
          &(recvCnts[0]), &(recvDisps[0]), true);

      std::vector<double> allPelletMaxZdisps(allPelletIds.size());
      d_pelletStackComm.allGather<double>(&(myMaxZdisps[0]), d_pelletIds.size(), &(allPelletMaxZdisps[0]), 
          &(recvCnts[0]), &(recvDisps[0]), true);

      finalMaxZdispsList.resize(d_totalNumberOfPellets, 0.0);
      for(size_t i = 0; i < allPelletIds.size(); i++) {
        if(fabs(allPelletMaxZdisps[i]) > fabs(finalMaxZdispsList[allPelletIds[i]])) {
          finalMaxZdispsList[allPelletIds[i]] = allPelletMaxZdisps[i];
        }
      }//end for i

      for(size_t i = 1; i < d_totalNumberOfPellets; i++) {
        finalMaxZdispsList[i] += finalMaxZdispsList[i - 1];
      }//end for i
    }

  }
}


