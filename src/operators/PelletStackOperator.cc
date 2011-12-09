
#include "operators/PelletStackOperator.h"
#include "operators/map/NodeToNodeMap.h"

namespace AMP {
  namespace Operator {

    PelletStackOperator :: PelletStackOperator(const boost::shared_ptr<OperatorParameters> & params)
      : Operator(params) {
        d_totalNumberOfPellets = (params->d_db)->getInteger("TOTAL_NUMBER_OF_PELLETS");
        d_useSerial = (params->d_db)->getBool("USE_SERIAL");
        d_onlyZcorrection = (params->d_db)->getBool("ONLY_Z_CORRECTION");
        d_masterId = (params->d_db)->getInteger("MASTER");
        d_slaveId = (params->d_db)->getInteger("SLAVE");
        if((params->d_db)->keyExists("SCALING_FACTOR")) {
          d_useScaling = true;
          d_scalingFactor = (params->d_db)->getDouble("SCALING_FACTOR");
        } else {
          d_useScaling = false;
        }
        d_currentPellet = static_cast<unsigned int>(-1);
      }

    void PelletStackOperator :: setFrozenVectorForMaps(AMP::LinearAlgebra::Vector::shared_ptr vec) {
      d_frozenVectorForMaps = vec;
    }

    void PelletStackOperator :: setMaps(boost::shared_ptr<AMP::Operator::AsyncMapColumnOperator> maps) {
      d_n2nMaps = maps;
    }

    void PelletStackOperator :: setPelletStackComm(AMP_MPI comm) {
      d_pelletStackComm = comm;
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

    void PelletStackOperator :: setCurrentPellet(unsigned int pellId) {
      d_currentPellet = pellId;
    }

    unsigned int PelletStackOperator :: getTotalNumberOfPellets() {
      return d_totalNumberOfPellets;
    }

    void PelletStackOperator :: setVariables(AMP::LinearAlgebra::Variable::shared_ptr rhs, 
        AMP::LinearAlgebra::Variable::shared_ptr sol) {
      d_rhsVar = boost::dynamic_pointer_cast<AMP::LinearAlgebra::MultiVariable>(rhs);
      d_solVar = boost::dynamic_pointer_cast<AMP::LinearAlgebra::MultiVariable>(sol);
    }

    void PelletStackOperator :: setLocalMeshes(std::vector<AMP::Mesh::MeshManager::Adapter::shared_ptr> inp) {
      d_meshes = inp;
    }

    void PelletStackOperator :: setLocalPelletIds(std::vector<unsigned int> inp) {
      d_pelletIds = inp;
    }

    int PelletStackOperator :: getLocalIndexForPellet(unsigned int pellId) {
      for(size_t i = 0; i < d_pelletIds.size(); i++) {
        if(d_pelletIds[i] == pellId) {
          return i;
        }
      }//end for i
      return -1;
    }

    void PelletStackOperator :: apply(const AMP::LinearAlgebra::Vector::shared_ptr &f,
        const AMP::LinearAlgebra::Vector::shared_ptr &u, AMP::LinearAlgebra::Vector::shared_ptr &r,
        const double a, const double b) {
      if(d_useSerial) {
        applySerial(f, u, r);
      } else if(d_onlyZcorrection) {
        applyOnlyZcorrection(r);
      } else {
        applyXYZcorrection(f, u, r);
      }
    }

    void PelletStackOperator :: applyUnscaling(AMP::LinearAlgebra::Vector::shared_ptr f) {
      for(size_t i = 0; i < d_pelletIds.size(); i++) {
        if(d_pelletIds[i] > 0) {
          AMP::LinearAlgebra::Variable::shared_ptr currVar = d_rhsVar->getVariable(i);
          AMP::LinearAlgebra::Vector::shared_ptr subF = f->subsetVectorForVariable(currVar);
          AMP::Mesh::DOFMap::shared_ptr dof_map = d_meshes[i]->getDOFMap(currVar);
          AMP::Mesh::MeshManager::Adapter::OwnedBoundaryNodeIterator bnd = d_meshes[i]->beginOwnedBoundary( d_slaveId );
          AMP::Mesh::MeshManager::Adapter::OwnedBoundaryNodeIterator end_bnd = d_meshes[i]->endOwnedBoundary( d_slaveId );
          std::vector<unsigned int> dofIds(3);
          dofIds[0] = 0; dofIds[1] = 1; dofIds[2] = 2;
          for( ; bnd != end_bnd; ++bnd) {
            std::vector<unsigned int> bndGlobalIds;
            dof_map->getDOFs(*bnd, bndGlobalIds, dofIds);
            for(unsigned int j = 0; j < bndGlobalIds.size(); j++) {
              double val = subF->getLocalValueByGlobalID( bndGlobalIds[j] );
              subF->setLocalValueByGlobalID(bndGlobalIds[j], val/d_scalingFactor);
            }//end for j
          }//end for bnd
        }
      }//end for i
    }

    void PelletStackOperator :: applyOnlyZcorrection(AMP::LinearAlgebra::Vector::shared_ptr &u) {
      std::vector<double> finalMaxZdispsList;
      computeZscan(u, finalMaxZdispsList); 
      for(size_t i = 0; i < d_pelletIds.size(); i++) {
        if(d_pelletIds[i] > 0) {
          AMP::LinearAlgebra::Variable::shared_ptr currVar = d_solVar->getVariable(i);
          AMP::LinearAlgebra::Vector::shared_ptr subU = u->subsetVectorForVariable(currVar);
          AMP::LinearAlgebra::Vector::shared_ptr zVec = subU->select( AMP::LinearAlgebra::VS_Stride("Z", 2, 3) , "Z" );
          zVec->addScalar(zVec, finalMaxZdispsList[d_pelletIds[i] - 1]);
        }
      }//end for i
    }

    void PelletStackOperator :: applyXYZcorrection(const AMP::LinearAlgebra::Vector::shared_ptr &f,
        const AMP::LinearAlgebra::Vector::shared_ptr &u, AMP::LinearAlgebra::Vector::shared_ptr  &r) {
      AMP::LinearAlgebra::Vector::shared_ptr nullVec;
      r->copyVector(f);
      d_n2nMaps->apply(nullVec, u, nullVec, 1.0, 0.0);
      for(size_t i = 0; i < d_pelletIds.size(); i++) {
        if(d_pelletIds[i] > 0) {
          AMP::LinearAlgebra::Variable::shared_ptr currVar = d_rhsVar->getVariable(i);
          AMP::LinearAlgebra::Vector::shared_ptr subU = d_frozenVectorForMaps->subsetVectorForVariable(currVar);
          AMP::LinearAlgebra::Vector::shared_ptr subR = r->subsetVectorForVariable(currVar);
          AMP::Mesh::DOFMap::shared_ptr dof_map = d_meshes[i]->getDOFMap(currVar);
          AMP::Mesh::MeshManager::Adapter::OwnedBoundaryNodeIterator bnd = d_meshes[i]->beginOwnedBoundary( d_slaveId );
          AMP::Mesh::MeshManager::Adapter::OwnedBoundaryNodeIterator end_bnd = d_meshes[i]->endOwnedBoundary( d_slaveId );
          std::vector<unsigned int> dofIds(3);
          dofIds[0] = 0; dofIds[1] = 1; dofIds[2] = 2;
          for( ; bnd != end_bnd; ++bnd) {
            std::vector<unsigned int> bndGlobalIds;
            dof_map->getDOFs(*bnd, bndGlobalIds, dofIds);
            for(unsigned int j = 0; j < bndGlobalIds.size(); j++) {
              double val = subU->getLocalValueByGlobalID( bndGlobalIds[j] );
              subR->addLocalValueByGlobalID(bndGlobalIds[j], val);
            }//end for j
          }//end for bnd
        }
      }//end for i
      std::vector<double> finalMaxZdispsList;
      computeZscan(u, finalMaxZdispsList); 
      for(size_t i = 0; i < d_pelletIds.size(); i++) {
        if(d_pelletIds[i] > 1) {
          AMP::LinearAlgebra::Variable::shared_ptr currVar = d_rhsVar->getVariable(i);
          AMP::LinearAlgebra::Vector::shared_ptr subR = r->subsetVectorForVariable(currVar);
          AMP::Mesh::DOFMap::shared_ptr dof_map = d_meshes[i]->getDOFMap(currVar);
          AMP::Mesh::MeshManager::Adapter::OwnedBoundaryNodeIterator bnd = d_meshes[i]->beginOwnedBoundary( d_slaveId );
          AMP::Mesh::MeshManager::Adapter::OwnedBoundaryNodeIterator end_bnd = d_meshes[i]->endOwnedBoundary( d_slaveId );
          std::vector<unsigned int> dofIds(1);
          dofIds[0] = 2;
          for( ; bnd != end_bnd; ++bnd) {
            std::vector<unsigned int> bndGlobalIds;
            dof_map->getDOFs(*bnd, bndGlobalIds, dofIds);
            for(unsigned int j = 0; j < bndGlobalIds.size(); j++) {
              subR->addLocalValueByGlobalID(bndGlobalIds[j], finalMaxZdispsList[d_pelletIds[i] - 2]);
            }//end for j
          }//end for bnd
        }
      }//end for i
    }

    void PelletStackOperator :: computeZscan(const AMP::LinearAlgebra::Vector::shared_ptr &u, 
        std::vector<double> &finalMaxZdispsList) {
      std::vector<double> myMaxZdisps(d_pelletIds.size(), 0.0);
      for(size_t i = 0; i < d_pelletIds.size(); i++) {
        AMP::LinearAlgebra::Variable::shared_ptr currVar = d_solVar->getVariable(i);
        AMP::LinearAlgebra::Vector::shared_ptr subU = u->subsetVectorForVariable(currVar);
        AMP::Mesh::DOFMap::shared_ptr dof_map = d_meshes[i]->getDOFMap(currVar);
        AMP::Mesh::MeshManager::Adapter::OwnedBoundaryNodeIterator bnd = d_meshes[i]->beginOwnedBoundary( d_masterId );
        AMP::Mesh::MeshManager::Adapter::OwnedBoundaryNodeIterator end_bnd = d_meshes[i]->endOwnedBoundary( d_masterId );
        std::vector<unsigned int> dofIds(1);
        dofIds[0] = 2; 
        for( ; bnd != end_bnd; ++bnd) {
          std::vector<unsigned int> bndGlobalIds;
          dof_map->getDOFs(*bnd, bndGlobalIds, dofIds);
          for(unsigned int j = 0; j < bndGlobalIds.size(); j++) {
            double val = subU->getLocalValueByGlobalID( bndGlobalIds[j] );
            if(fabs(myMaxZdisps[i]) < fabs(val)) {
              myMaxZdisps[i] = val;
            }
          }//end for j
        }//end for bnd
      }//end for i

      std::vector<int> recvCnts(d_pelletStackComm.getSize()); 
      d_pelletStackComm.allGather<int>(d_pelletIds.size(), &(recvCnts[0]));

      std::vector<int> recvDisps(recvCnts.size());
      recvDisps[0] = 0;
      for(size_t i = 1; i < recvDisps.size(); i++) {
        recvDisps[i] = recvDisps[i - 1] + recvCnts[i - 1];
      }//end for i

      std::vector<unsigned int> allPelletIds((*(recvDisps.end())) + (*(recvCnts.end())));
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

    void PelletStackOperator :: applySerial(const AMP::LinearAlgebra::Vector::shared_ptr &f,
        const AMP::LinearAlgebra::Vector::shared_ptr &u, AMP::LinearAlgebra::Vector::shared_ptr &r) {
      AMP_ASSERT(d_currentPellet > 0);
      AMP::LinearAlgebra::Vector::shared_ptr nullVec;
      int currPellIdx = getLocalIndexForPellet(d_currentPellet);
      int prevPellIdx = getLocalIndexForPellet(d_currentPellet - 1);
      int numMaps = d_n2nMaps->getNumberOfOperators();
      if(currPellIdx != -1) {
        AMP::LinearAlgebra::Variable::shared_ptr currVar = d_solVar->getVariable(currPellIdx);
        for(int m = 0; m < numMaps; m++) {
          boost::shared_ptr<AMP::Operator::NodeToNodeMap> currMap = boost::dynamic_pointer_cast<
            AMP::Operator::NodeToNodeMap>(d_n2nMaps->getOperator(m));
          AMP::LinearAlgebra::Variable::shared_ptr currMapVar = currMap->getInputVariable();
          if((*currMapVar) == (*currVar)) {
            if(currMap->isMaster() == false) {
              currMap->applyStart(nullVec, u, nullVec, 1.0, 0.0);
              break;
            }
          }
        }//end for m
      }
      if(prevPellIdx != -1) {
        AMP::LinearAlgebra::Variable::shared_ptr currVar = d_solVar->getVariable(prevPellIdx);
        for(int m = 0; m < numMaps; m++) {
          boost::shared_ptr<AMP::Operator::NodeToNodeMap> currMap = boost::dynamic_pointer_cast<
            AMP::Operator::NodeToNodeMap>(d_n2nMaps->getOperator(m));
          AMP::LinearAlgebra::Variable::shared_ptr currMapVar = currMap->getInputVariable();
          if((*currMapVar) == (*currVar)) {
            if(currMap->isMaster() == true) {
              currMap->applyStart(nullVec, u, nullVec, 1.0, 0.0);
              break;
            }
          }
        }//end for m
      }
      if(currPellIdx != -1) {
        AMP::LinearAlgebra::Variable::shared_ptr currVar = d_solVar->getVariable(currPellIdx);
        for(int m = 0; m < numMaps; m++) {
          boost::shared_ptr<AMP::Operator::NodeToNodeMap> currMap = boost::dynamic_pointer_cast<
            AMP::Operator::NodeToNodeMap>(d_n2nMaps->getOperator(m));
          AMP::LinearAlgebra::Variable::shared_ptr currMapVar = currMap->getInputVariable();
          if((*currMapVar) == (*currVar)) {
            if(currMap->isMaster() == false) {
              currMap->applyFinish(nullVec, u, nullVec, 1.0, 0.0);
              break;
            }
          }
        }//end for m
      }
      if(prevPellIdx != -1) {
        AMP::LinearAlgebra::Variable::shared_ptr currVar = d_solVar->getVariable(prevPellIdx);
        for(int m = 0; m < numMaps; m++) {
          boost::shared_ptr<AMP::Operator::NodeToNodeMap> currMap = boost::dynamic_pointer_cast<
            AMP::Operator::NodeToNodeMap>(d_n2nMaps->getOperator(m));
          AMP::LinearAlgebra::Variable::shared_ptr currMapVar = currMap->getInputVariable();
          if((*currMapVar) == (*currVar)) {
            if(currMap->isMaster() == true) {
              currMap->applyFinish(nullVec, u, nullVec, 1.0, 0.0);
              break;
            }
          }
        }//end for m
      }
      if(currPellIdx != -1) {
        AMP::LinearAlgebra::Variable::shared_ptr currVar = d_rhsVar->getVariable(currPellIdx);
        AMP::LinearAlgebra::Vector::shared_ptr subF = f->subsetVectorForVariable(currVar);
        AMP::LinearAlgebra::Vector::shared_ptr subR = r->subsetVectorForVariable(currVar);
        AMP::LinearAlgebra::Vector::shared_ptr subU = d_frozenVectorForMaps->subsetVectorForVariable(currVar);
        subR->copyVector(subF);
        AMP::Mesh::DOFMap::shared_ptr dof_map = d_meshes[currPellIdx]->getDOFMap(currVar);
        AMP::Mesh::MeshManager::Adapter::OwnedBoundaryNodeIterator bnd = d_meshes[currPellIdx]->beginOwnedBoundary( d_slaveId );
        AMP::Mesh::MeshManager::Adapter::OwnedBoundaryNodeIterator end_bnd = d_meshes[currPellIdx]->endOwnedBoundary( d_slaveId );
        std::vector<unsigned int> dofIds(3);
        dofIds[0] = 0; dofIds[1] = 1; dofIds[2] = 2;
        for( ; bnd != end_bnd; ++bnd) {
          std::vector<unsigned int> bndGlobalIds;
          dof_map->getDOFs(*bnd, bndGlobalIds, dofIds);
          for(unsigned int j = 0; j < bndGlobalIds.size(); j++) {
            double val = subU->getLocalValueByGlobalID( bndGlobalIds[j] );
            subR->addLocalValueByGlobalID(bndGlobalIds[j], val);
          }//end for j
        }//end for bnd
      }
    }

  }
}


