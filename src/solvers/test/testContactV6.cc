
#include "utils/InputManager.h"
#include "utils/AMPManager.h"
#include "utils/UnitTest.h"
#include "utils/Utilities.h"

#include "boost/shared_ptr.hpp"

#include <iostream>
#include <cstdio>
#include <cstring>
#include <vector>
#include <algorithm>
#include <cmath>

#include "ampmesh/MeshManager.h"
#include "ampmesh/MeshAdapter.h"
#include "ampmesh/SiloIO.h"

#include "operators/ColumnOperator.h"
#include "operators/OperatorBuilder.h"
#include "operators/mechanics/MechanicsLinearFEOperator.h"
#include "operators/boundary/DirichletVectorCorrection.h"
#include "operators/mechanics/ConstructLinearMechanicsRHSVector.h"

#include "fe_interface.h"
#include "cell_hex8.h"

#include "ml_include.h"
#include "solvers/MLoptions.h"
#include "indexHolder.h"

#include "ContactSearchUtils.h"

void deformMesh(AMP::Mesh::MeshManager::Adapter::shared_ptr meshAdapter,
    AMP::LinearAlgebra::Vector::shared_ptr mechSolVec) {
  AMP::Mesh::DOFMap::shared_ptr dof_map = meshAdapter->getDOFMap( mechSolVec->getVariable() );

  AMP::Mesh::MeshManager::Adapter::OwnedNodeIterator nd = meshAdapter->beginOwnedNode();
  AMP::Mesh::MeshManager::Adapter::OwnedNodeIterator end_nd = meshAdapter->endOwnedNode();

  std::vector <unsigned int> dofIds(3);
  dofIds[0] = 0; dofIds[1] = 1; dofIds[2] = 2;

  for( ; nd != end_nd; ++nd) {
    std::vector<unsigned int> ndGlobalIds;
    dof_map->getDOFs(*nd, ndGlobalIds, dofIds);

    double xDisp = mechSolVec->getValueByGlobalID(ndGlobalIds[0]);
    double yDisp = mechSolVec->getValueByGlobalID(ndGlobalIds[1]);
    double zDisp = mechSolVec->getValueByGlobalID(ndGlobalIds[2]);

    nd->translate(xDisp, yDisp, zDisp);
  }//end for nd
}

ML_Aggregate *agg_object;
ML* ml_object;

struct MyData {
  AMP::Mesh::MeshManager::Adapter::shared_ptr slaveMeshAdapter;
  AMP::Mesh::MeshManager::Adapter::shared_ptr masterMeshAdapter;
  boost::shared_ptr<AMP::Operator::MechanicsLinearFEOperator> masterOperator;
  boost::shared_ptr<AMP::Operator::MechanicsLinearFEOperator> slaveOperator;
  boost::shared_ptr<AMP::Operator::ColumnOperator> columnOperator;
  AMP::LinearAlgebra::Variable::shared_ptr slaveVar;
  AMP::LinearAlgebra::Variable::shared_ptr masterVar;
  AMP::LinearAlgebra::Vector::shared_ptr matVecInVec;
  AMP::LinearAlgebra::Vector::shared_ptr matVecOutVec;
  std::vector<double> pcSolVec;
  std::vector<double> pcRhsVec;
  std::vector<unsigned int> ml2MasterOrSlave;
  std::vector<unsigned int> master2Ml;
  std::vector<unsigned int> slave2Ml;
  std::vector<std::vector<unsigned int> > master2SlaveIds;
  std::vector<std::vector<double> > master2SlaveFactors;
  std::vector<unsigned int> slaveId2SlaveListIndex;
  std::vector<unsigned int> slaveId2Dof;
  std::vector<double> slaveOffsets;
  std::vector<size_t> slaveNodes;
  std::vector<std::vector<size_t> > slave2MasterNodes;
  std::vector<std::vector<double> > slave2MasterFactors;
  std::vector<unsigned int> masterBndIds;
  std::vector<std::vector<unsigned int> > masterBndDofIds;
  std::vector<std::vector<double> > masterBndDofVals;
  std::vector<unsigned int> slaveBndIds;
  std::vector<std::vector<unsigned int> > slaveBndDofIds;
  std::vector<std::vector<double> > slaveBndDofVals;
  unsigned int mlSize;
  size_t numMasterDirichletDofs;
  size_t numSlaveDirichletDofs;
} obj;

void freeMyData() {
  obj.slaveMeshAdapter.reset();
  obj.masterMeshAdapter.reset();
  obj.masterOperator.reset();
  obj.slaveOperator.reset();
  obj.columnOperator.reset();
  obj.slaveVar.reset();
  obj.masterVar.reset();
  obj.matVecInVec.reset();
  obj.matVecOutVec.reset();
}

void computeSlaveId2SlaveListIndex() {
  size_t totalSlaveDofs = (3*(obj.slaveMeshAdapter->numLocalNodes()));

  obj.slaveId2SlaveListIndex.resize(totalSlaveDofs);
  obj.slaveId2Dof.resize(totalSlaveDofs);

  for(size_t i = 0; i < totalSlaveDofs; i++) {
    obj.slaveId2SlaveListIndex[i] = static_cast<unsigned int>(-1);
    obj.slaveId2Dof[i] = static_cast<unsigned int>(-1);
  }//end for i

  AMP::Mesh::DOFMap::shared_ptr slave_dof_map = obj.slaveMeshAdapter->getDOFMap(obj.slaveVar);

  std::vector<unsigned int> dofs(3);
  dofs[0] = 0; dofs[1] = 1; dofs[2] = 2;

  for(size_t i = 0; i < obj.slaveNodes.size(); i++) {
    AMP::Mesh::LibMeshNode slaveNd =  obj.slaveMeshAdapter->getNode( obj.slaveNodes[i] );
    std::vector<unsigned int> slaveGlobalIds;
    slave_dof_map->getDOFs(slaveNd, slaveGlobalIds, dofs);
    for(size_t k = 0; k < dofs.size(); k++) {
      obj.slaveId2SlaveListIndex[slaveGlobalIds[k]] = i;
      obj.slaveId2Dof[slaveGlobalIds[k]] = dofs[k];
    }//end for k
  }//end for i
}

void computeMaster2Slaves() {
  size_t totalMasterDofs = (3*(obj.masterMeshAdapter->numLocalNodes()));

  obj.master2SlaveIds.resize(totalMasterDofs);
  obj.master2SlaveFactors.resize(totalMasterDofs);

  AMP::Mesh::DOFMap::shared_ptr slave_dof_map = obj.slaveMeshAdapter->getDOFMap(obj.slaveVar);
  AMP::Mesh::DOFMap::shared_ptr master_dof_map = obj.masterMeshAdapter->getDOFMap(obj.masterVar);

  std::vector<unsigned int> dofs(3);
  dofs[0] = 0; dofs[1] = 1; dofs[2] = 2;

  for(size_t i = 0; i < obj.slaveNodes.size(); i++) {
    AMP::Mesh::LibMeshNode slaveNd =  obj.slaveMeshAdapter->getNode( obj.slaveNodes[i] );
    std::vector<unsigned int> slaveGlobalIds;
    slave_dof_map->getDOFs(slaveNd, slaveGlobalIds, dofs);
    for(size_t j = 0; j < obj.slave2MasterNodes[i].size(); j++) {
      AMP::Mesh::LibMeshNode masterNd =  obj.masterMeshAdapter->getNode( obj.slave2MasterNodes[i][j] );
      std::vector<unsigned int> masterGlobalIds;
      master_dof_map->getDOFs(masterNd, masterGlobalIds, dofs);
      for(size_t k = 0; k < dofs.size(); k++) {
        obj.master2SlaveIds[masterGlobalIds[k]].push_back(slaveGlobalIds[k]);
        obj.master2SlaveFactors[masterGlobalIds[k]].push_back(obj.slave2MasterFactors[i][j]);
      }//end for k
    }//end for j
  }//end for i
}

void masterToSlaveCorrection(AMP::LinearAlgebra::Vector::shared_ptr vec) {
  AMP::Mesh::DOFMap::shared_ptr slave_dof_map = obj.slaveMeshAdapter->getDOFMap(obj.slaveVar);

  AMP::Mesh::DOFMap::shared_ptr master_dof_map = obj.masterMeshAdapter->getDOFMap(obj.masterVar);

  AMP::LinearAlgebra::Vector::shared_ptr slaveVec = vec->subsetVectorForVariable(obj.slaveVar);

  AMP::LinearAlgebra::Vector::shared_ptr masterVec = vec->subsetVectorForVariable(obj.masterVar);

  std::vector<unsigned int> dofs(3);
  dofs[0] = 0; dofs[1] = 1; dofs[2] = 2;

  for(size_t i = 0; i < obj.slaveNodes.size(); i++) {
    AMP::Mesh::LibMeshNode slaveNd =  obj.slaveMeshAdapter->getNode( obj.slaveNodes[i] );
    std::vector<unsigned int> slaveGlobalIds;
    slave_dof_map->getDOFs(slaveNd, slaveGlobalIds, dofs);
    double slaveVal[] = {0.0, 0.0, 0.0};
    for(size_t j = 0; j < obj.slave2MasterNodes[i].size(); j++) {
      AMP::Mesh::LibMeshNode masterNd =  obj.masterMeshAdapter->getNode( obj.slave2MasterNodes[i][j] );
      std::vector<unsigned int> masterGlobalIds;
      master_dof_map->getDOFs(masterNd, masterGlobalIds, dofs);
      for(size_t k = 0; k < dofs.size(); k++) {
        double masterVal = masterVec->getValueByGlobalID(masterGlobalIds[k]);
        slaveVal[k] += (obj.slave2MasterFactors[i][j]*masterVal);
      }//end for k
    }//end for j
    for(size_t k = 0; k < dofs.size(); k++) {
      slaveVec->setLocalValueByGlobalID(slaveGlobalIds[k], slaveVal[k]);
    }//end for k
  }//end for i

  vec->makeConsistent( AMP::LinearAlgebra::Vector::CONSISTENT_SET );
}

void slaveToMasterCorrection(AMP::LinearAlgebra::Vector::shared_ptr vec) {
  AMP::Mesh::DOFMap::shared_ptr slave_dof_map = obj.slaveMeshAdapter->getDOFMap(obj.slaveVar);

  AMP::Mesh::DOFMap::shared_ptr master_dof_map = obj.masterMeshAdapter->getDOFMap(obj.masterVar);

  AMP::LinearAlgebra::Vector::shared_ptr slaveVec = vec->subsetVectorForVariable(obj.slaveVar);

  AMP::LinearAlgebra::Vector::shared_ptr masterVec = vec->subsetVectorForVariable(obj.masterVar);

  std::vector<unsigned int> dofs(3);
  dofs[0] = 0; dofs[1] = 1; dofs[2] = 2;

  for(size_t i = 0; i < obj.slaveNodes.size(); i++) {
    AMP::Mesh::LibMeshNode slaveNd =  obj.slaveMeshAdapter->getNode( obj.slaveNodes[i] );
    std::vector<unsigned int> slaveGlobalIds;
    slave_dof_map->getDOFs(slaveNd, slaveGlobalIds, dofs);
    for(size_t j = 0; j < obj.slave2MasterNodes[i].size(); j++) {
      AMP::Mesh::LibMeshNode masterNd =  obj.masterMeshAdapter->getNode( obj.slave2MasterNodes[i][j] );
      std::vector<unsigned int> masterGlobalIds;
      master_dof_map->getDOFs(masterNd, masterGlobalIds, dofs);
      for(size_t k = 0; k < dofs.size(); k++) {
        double slaveVal = slaveVec->getValueByGlobalID(slaveGlobalIds[k]);
        double masterVal = (obj.slave2MasterFactors[i][j]*slaveVal);
        masterVec->addValueByGlobalID(masterGlobalIds[k], masterVal);
      }//end for k
    }//end for j
  }//end for i

  vec->makeConsistent( AMP::LinearAlgebra::Vector::CONSISTENT_ADD );
}

void addSlaveOffsets(AMP::LinearAlgebra::Vector::shared_ptr vec) {
  AMP::Mesh::DOFMap::shared_ptr slave_dof_map = obj.slaveMeshAdapter->getDOFMap(obj.slaveVar);

  AMP::LinearAlgebra::Vector::shared_ptr slaveVec = vec->subsetVectorForVariable(obj.slaveVar);

  std::vector<unsigned int> dofs(3);
  dofs[0] = 0; dofs[1] = 1; dofs[2] = 2;

  for(size_t i = 0; i < obj.slaveNodes.size(); i++) {
    AMP::Mesh::LibMeshNode slaveNd =  obj.slaveMeshAdapter->getNode( obj.slaveNodes[i] );
    std::vector<unsigned int> slaveGlobalIds;
    slave_dof_map->getDOFs(slaveNd, slaveGlobalIds, dofs);
    for(size_t j = 0; j < dofs.size(); j++) {
      slaveVec->addValueByGlobalID(slaveGlobalIds[j], obj.slaveOffsets[(3*i) + j]);
    }//end for j
  }//end for i

  vec->makeConsistent( AMP::LinearAlgebra::Vector::CONSISTENT_ADD );
}

void setSlaveToZero(AMP::LinearAlgebra::Vector::shared_ptr vec) {
  AMP::Mesh::DOFMap::shared_ptr slave_dof_map = obj.slaveMeshAdapter->getDOFMap(obj.slaveVar);

  AMP::LinearAlgebra::Vector::shared_ptr slaveVec = vec->subsetVectorForVariable(obj.slaveVar);

  std::vector<unsigned int> dofs(3);
  dofs[0] = 0; dofs[1] = 1; dofs[2] = 2;

  for(size_t i = 0; i < obj.slaveNodes.size(); i++) {
    AMP::Mesh::LibMeshNode slaveNd =  obj.slaveMeshAdapter->getNode( obj.slaveNodes[i] );
    std::vector<unsigned int> slaveGlobalIds;
    slave_dof_map->getDOFs(slaveNd, slaveGlobalIds, dofs);
    for(size_t j = 0; j < dofs.size(); j++) {
      slaveVec->setLocalValueByGlobalID(slaveGlobalIds[j], 0.0);
    }//end for j
  }//end for i

  vec->makeConsistent( AMP::LinearAlgebra::Vector::CONSISTENT_SET );
}

void computeSlaveOffsets() {
  obj.slaveOffsets.resize(3*(obj.slaveNodes.size()));
  for(size_t i = 0; i < obj.slaveNodes.size(); i++) {
    AMP::Mesh::LibMeshNode slaveNd =  obj.slaveMeshAdapter->getNode( obj.slaveNodes[i] );
    double xTrue = slaveNd.x();
    double yTrue = slaveNd.y();
    double zTrue = slaveNd.z();
    double xProj = 0.0;
    double yProj = 0.0;
    double zProj = 0.0;
    for(size_t j = 0; j < obj.slave2MasterNodes[i].size(); j++) {
      AMP::Mesh::LibMeshNode masterNd =  obj.masterMeshAdapter->getNode( obj.slave2MasterNodes[i][j] );
      xProj += (obj.slave2MasterFactors[i][j]*(masterNd.x()));
      yProj += (obj.slave2MasterFactors[i][j]*(masterNd.y()));
      zProj += (obj.slave2MasterFactors[i][j]*(masterNd.z()));
    }//end for j
    obj.slaveOffsets[(3*i)] = (xProj - xTrue);
    obj.slaveOffsets[(3*i) + 1] = (yProj - yTrue);
    obj.slaveOffsets[(3*i) + 2] = (zProj - zTrue);
  }//end for i
}

void computeMl2MasterOrSlave() {
  std::vector<unsigned int> tmpVec;

  size_t totalMasterDofs = (3*(obj.masterMeshAdapter->numLocalNodes()));
  size_t totalSlaveDofs = (3*(obj.slaveMeshAdapter->numLocalNodes()));

  tmpVec.resize(totalMasterDofs + totalSlaveDofs);

  for(size_t i = 0; i < totalMasterDofs; i++) {
    tmpVec[i] = i;
  }//end for i
  for(size_t i = 0; i < totalSlaveDofs; i++) {
    tmpVec[totalMasterDofs + i] = i;
  }//end for i

  AMP::Mesh::DOFMap::shared_ptr master_dof_map = obj.masterMeshAdapter->getDOFMap(obj.masterVar);

  for(size_t i = 0; i < obj.masterBndIds.size(); i++) {
    AMP::Mesh::MeshManager::Adapter::OwnedBoundaryNodeIterator bnd = 
      obj.masterMeshAdapter->beginOwnedBoundary( obj.masterBndIds[i] );
    AMP::Mesh::MeshManager::Adapter::OwnedBoundaryNodeIterator end_bnd = 
      obj.masterMeshAdapter->endOwnedBoundary( obj.masterBndIds[i] );

    for( ; bnd != end_bnd; ++bnd) {
      std::vector<unsigned int> bndGlobalIds;
      master_dof_map->getDOFs(*bnd, bndGlobalIds, obj.masterBndDofIds[i]);
      for(size_t j = 0; j < bndGlobalIds.size(); j++) {
        tmpVec[bndGlobalIds[j]] = static_cast<unsigned int>(-1);
      }//end for j
    }//end for bnd
  }//end for i

  AMP::Mesh::DOFMap::shared_ptr slave_dof_map = obj.slaveMeshAdapter->getDOFMap(obj.slaveVar);

  for(size_t i = 0; i < obj.slaveBndIds.size(); i++) {
    AMP::Mesh::MeshManager::Adapter::OwnedBoundaryNodeIterator bnd = 
      obj.slaveMeshAdapter->beginOwnedBoundary( obj.slaveBndIds[i] );
    AMP::Mesh::MeshManager::Adapter::OwnedBoundaryNodeIterator end_bnd = 
      obj.slaveMeshAdapter->endOwnedBoundary( obj.slaveBndIds[i] );

    for( ; bnd != end_bnd; ++bnd) {
      std::vector<unsigned int> bndGlobalIds;
      slave_dof_map->getDOFs(*bnd, bndGlobalIds, obj.slaveBndDofIds[i]);
      for(size_t j = 0; j < bndGlobalIds.size(); j++) {
        tmpVec[totalMasterDofs + bndGlobalIds[j]] = static_cast<unsigned int>(-1);
      }//end for j
    }//end for bnd
  }//end for i

  std::vector<unsigned int> dofs(3);
  dofs[0] = 0; dofs[1] = 1; dofs[2] = 2;

  for(size_t i = 0; i < obj.slaveNodes.size(); i++) {
    AMP::Mesh::LibMeshNode slaveNd =  obj.slaveMeshAdapter->getNode( obj.slaveNodes[i] );
    std::vector<unsigned int> slaveGlobalIds;
    slave_dof_map->getDOFs(slaveNd, slaveGlobalIds, dofs);
    for(size_t j = 0; j < dofs.size(); j++) {
      tmpVec[totalMasterDofs + slaveGlobalIds[j]] = static_cast<unsigned int>(-1);
    }//end for j
  }//end for i

  obj.ml2MasterOrSlave.clear();
  for(size_t i = 0; i < tmpVec.size(); i++) {
    if(tmpVec[i] != static_cast<unsigned int>(-1)) {
      obj.ml2MasterOrSlave.push_back(tmpVec[i]);
    }
  }//end for i
}

void computeSlave2Ml() {
  size_t totalMasterDofs = (3*(obj.masterMeshAdapter->numLocalNodes()));
  size_t totalSlaveDofs = (3*(obj.slaveMeshAdapter->numLocalNodes()));

  obj.slave2Ml.resize(totalSlaveDofs);

  for(size_t i = 0; i < totalSlaveDofs; i++) {
    obj.slave2Ml[i] = static_cast<unsigned int>(-1);
  }//end for i

  size_t numMlMasterDofs = totalMasterDofs - obj.numMasterDirichletDofs;

  for(size_t i = numMlMasterDofs; i < obj.mlSize; i++) {
    obj.slave2Ml[obj.ml2MasterOrSlave[i]] = i;
  }//end for i
}

void computeMaster2Ml() {
  size_t totalMasterDofs = (3*(obj.masterMeshAdapter->numLocalNodes()));

  obj.master2Ml.resize(totalMasterDofs);

  for(size_t i = 0; i < totalMasterDofs; i++) {
    obj.master2Ml[i] = static_cast<unsigned int>(-1);
  }//end for i

  size_t numMlMasterDofs = totalMasterDofs - obj.numMasterDirichletDofs;

  for(size_t i = 0; i < numMlMasterDofs; i++) {
    obj.master2Ml[obj.ml2MasterOrSlave[i]] = i;
  }//end for i
}

void setSlaveDirichletDofsToZero(AMP::LinearAlgebra::Vector::shared_ptr vec) {
  AMP::Mesh::DOFMap::shared_ptr slave_dof_map = obj.slaveMeshAdapter->getDOFMap(obj.slaveVar);

  AMP::LinearAlgebra::Vector::shared_ptr slaveVec =  vec->subsetVectorForVariable(obj.slaveVar);

  for(size_t i = 0; i < obj.slaveBndIds.size(); i++) {
    AMP::Mesh::MeshManager::Adapter::OwnedBoundaryNodeIterator bnd = 
      obj.slaveMeshAdapter->beginOwnedBoundary( obj.slaveBndIds[i] );
    AMP::Mesh::MeshManager::Adapter::OwnedBoundaryNodeIterator end_bnd = 
      obj.slaveMeshAdapter->endOwnedBoundary( obj.slaveBndIds[i] );

    for( ; bnd != end_bnd; ++bnd) {
      std::vector<unsigned int> bndGlobalIds;
      slave_dof_map->getDOFs(*bnd, bndGlobalIds, obj.slaveBndDofIds[i]);

      for(size_t j = 0; j < bndGlobalIds.size(); j++) {
        slaveVec->setLocalValueByGlobalID(bndGlobalIds[j], 0.0);
      }//end for j
    }//end for bnd
  }//end for i

  vec->makeConsistent( AMP::LinearAlgebra::Vector::CONSISTENT_SET );
}

/*
   void setSlaveDirichletVals(AMP::LinearAlgebra::Vector::shared_ptr vec) {
   AMP::Mesh::DOFMap::shared_ptr slave_dof_map = obj.slaveMeshAdapter->getDOFMap(obj.slaveVar);

   AMP::LinearAlgebra::Vector::shared_ptr slaveVec =  vec->subsetVectorForVariable(obj.slaveVar);

   for(size_t i = 0; i < obj.slaveBndIds.size(); i++) {
   AMP::Mesh::MeshManager::Adapter::OwnedBoundaryNodeIterator bnd = 
   obj.slaveMeshAdapter->beginOwnedBoundary( obj.slaveBndIds[i] );
   AMP::Mesh::MeshManager::Adapter::OwnedBoundaryNodeIterator end_bnd = 
   obj.slaveMeshAdapter->endOwnedBoundary( obj.slaveBndIds[i] );

   for( ; bnd != end_bnd; ++bnd) {
   std::vector<unsigned int> bndGlobalIds;
   slave_dof_map->getDOFs(*bnd, bndGlobalIds, obj.slaveBndDofIds[i]);

   for(size_t j = 0; j < bndGlobalIds.size(); j++) {
   slaveVec->setLocalValueByGlobalID(bndGlobalIds[j], obj.slaveBndDofVals[i][j]);
   }//end for j
   }//end for bnd
   }//end for i

   vec->makeConsistent( AMP::LinearAlgebra::Vector::CONSISTENT_SET );
   }
   */

void setMasterDirichletDofsToZero(AMP::LinearAlgebra::Vector::shared_ptr vec) {
  AMP::Mesh::DOFMap::shared_ptr master_dof_map = obj.masterMeshAdapter->getDOFMap(obj.masterVar);

  AMP::LinearAlgebra::Vector::shared_ptr masterVec =  vec->subsetVectorForVariable(obj.masterVar);

  for(size_t i = 0; i < obj.masterBndIds.size(); i++) {
    AMP::Mesh::MeshManager::Adapter::OwnedBoundaryNodeIterator bnd = 
      obj.masterMeshAdapter->beginOwnedBoundary( obj.masterBndIds[i] );
    AMP::Mesh::MeshManager::Adapter::OwnedBoundaryNodeIterator end_bnd = 
      obj.masterMeshAdapter->endOwnedBoundary( obj.masterBndIds[i] );

    for( ; bnd != end_bnd; ++bnd) {
      std::vector<unsigned int> bndGlobalIds;
      master_dof_map->getDOFs(*bnd, bndGlobalIds, obj.masterBndDofIds[i]);

      for(size_t j = 0; j < bndGlobalIds.size(); j++) {
        masterVec->setLocalValueByGlobalID(bndGlobalIds[j], 0.0);
      }//end for j
    }//end for bnd
  }//end for i

  vec->makeConsistent( AMP::LinearAlgebra::Vector::CONSISTENT_SET );
}

/*
   void setMasterDirichletVals(AMP::LinearAlgebra::Vector::shared_ptr vec) {
   AMP::Mesh::DOFMap::shared_ptr master_dof_map = obj.masterMeshAdapter->getDOFMap(obj.masterVar);

   AMP::LinearAlgebra::Vector::shared_ptr masterVec =  vec->subsetVectorForVariable(obj.masterVar);

   for(size_t i = 0; i < obj.masterBndIds.size(); i++) {
   AMP::Mesh::MeshManager::Adapter::OwnedBoundaryNodeIterator bnd = 
   obj.masterMeshAdapter->beginOwnedBoundary( obj.masterBndIds[i] );
   AMP::Mesh::MeshManager::Adapter::OwnedBoundaryNodeIterator end_bnd = 
   obj.masterMeshAdapter->endOwnedBoundary( obj.masterBndIds[i] );

   for( ; bnd != end_bnd; ++bnd) {
   std::vector<unsigned int> bndGlobalIds;
   master_dof_map->getDOFs(*bnd, bndGlobalIds, obj.masterBndDofIds[i]);

   for(size_t j = 0; j < bndGlobalIds.size(); j++) {
   masterVec->setLocalValueByGlobalID(bndGlobalIds[j], obj.masterBndDofVals[i][j]);
   }//end for j
   }//end for bnd
   }//end for i

   vec->makeConsistent( AMP::LinearAlgebra::Vector::CONSISTENT_SET );
   }
   */

void copyML2Column(double* mlVec, AMP::LinearAlgebra::Vector::shared_ptr vec) {
  AMP::LinearAlgebra::Vector::shared_ptr slaveVec = vec->subsetVectorForVariable(obj.slaveVar);
  AMP::LinearAlgebra::Vector::shared_ptr masterVec = vec->subsetVectorForVariable(obj.masterVar);

  size_t totalMasterDofs = (3*(obj.masterMeshAdapter->numLocalNodes()));
  size_t numMlMasterDofs = totalMasterDofs - obj.numMasterDirichletDofs;

  for(size_t i = 0; i < numMlMasterDofs; i++) {
    masterVec->setValueByGlobalID(obj.ml2MasterOrSlave[i], mlVec[i]);
  }//end for i

  for(size_t i = numMlMasterDofs; i < obj.mlSize; i++) {
    slaveVec->setValueByGlobalID(obj.ml2MasterOrSlave[i], mlVec[i]);
  }//end for i

  vec->makeConsistent( AMP::LinearAlgebra::Vector::CONSISTENT_SET );
}

void copyColumn2ML(double* mlVec, AMP::LinearAlgebra::Vector::shared_ptr vec) {
  AMP::LinearAlgebra::Vector::shared_ptr slaveVec = vec->subsetVectorForVariable(obj.slaveVar);
  AMP::LinearAlgebra::Vector::shared_ptr masterVec = vec->subsetVectorForVariable(obj.masterVar);

  size_t totalMasterDofs = (3*(obj.masterMeshAdapter->numLocalNodes()));
  size_t numMlMasterDofs = totalMasterDofs - obj.numMasterDirichletDofs;

  for(size_t i = 0; i < numMlMasterDofs; i++) {
    mlVec[i] = masterVec->getValueByGlobalID(obj.ml2MasterOrSlave[i]);
  }//end for i

  for(size_t i = numMlMasterDofs; i < obj.mlSize; i++) {
    mlVec[i] = slaveVec->getValueByGlobalID(obj.ml2MasterOrSlave[i]);
  }//end for i
}

int myMatVec(ML_Operator *data, int in_length, double in[], int out_length, double out[]) {
  AMP::LinearAlgebra::Vector::shared_ptr nullVec;

  copyML2Column(in, obj.matVecInVec);

  setMasterDirichletDofsToZero(obj.matVecInVec);
  setSlaveDirichletDofsToZero(obj.matVecInVec);

  masterToSlaveCorrection(obj.matVecInVec);

  (obj.columnOperator)->apply(nullVec, obj.matVecInVec, obj.matVecOutVec, 1.0, 0.0);

  slaveToMasterCorrection(obj.matVecOutVec);

  copyColumn2ML(out, obj.matVecOutVec);

  return 0;
}

void myGetSingleRow(int row, std::vector<unsigned int> &cols, std::vector<double> &values) {
  cols.clear();
  values.clear();

}

int myGetMultipleRows(ML_Operator *data, int N_requested_rows, int requested_rows[],
    int allocated_space, int columns[], double values[], int row_lengths[] ) {
  int spaceRequired = 0;
  int cnt = 0;
  for(int i = 0; i < N_requested_rows; i++) {
    int row = requested_rows[i];
    std::vector<unsigned int> cols;
    std::vector<double> vals;

    myGetSingleRow(row, cols, vals);
    spaceRequired += cols.size();

    if(allocated_space >= spaceRequired) {
      for(size_t j = 0; j < cols.size(); j++) {
        columns[cnt] = cols[j];
        values[cnt] = vals[j];
        cnt++;
      }
      row_lengths[i] = cols.size();
    } else {
      return 0;
    }
  }

  return 1;
}

void pcMLapply(AMP::LinearAlgebra::Vector::shared_ptr fVec,
    AMP::LinearAlgebra::Vector::shared_ptr uVec) {
  double * solArr = &((obj.pcSolVec)[0]);
  double * rhsArr = &((obj.pcRhsVec)[0]);

  for(size_t i = 0; i < obj.mlSize; i++) {
    solArr[i] = 0;
  }

  copyColumn2ML(rhsArr, fVec);

  ML_Iterate(ml_object, solArr, rhsArr);

  copyML2Column(solArr, uVec);
}

void pcNoneApply(AMP::LinearAlgebra::Vector::shared_ptr fVec,
    AMP::LinearAlgebra::Vector::shared_ptr uVec) {
  //Identity PC
  uVec->copyVector(fVec);
}

/*
   void checkRadius(unsigned int id, AMP::Mesh::MeshManager::Adapter::shared_ptr meshAdapter) {
   AMP::Mesh::MeshManager::Adapter::OwnedBoundaryNodeIterator bnd =
   meshAdapter->beginOwnedBoundary( id );
   AMP::Mesh::MeshManager::Adapter::OwnedBoundaryNodeIterator end_bnd =
   meshAdapter->endOwnedBoundary( id );

   double minR;
   double maxR;
   double zAtMinR;
   double zAtMaxR;
   for(int i = 0; bnd != end_bnd; ++bnd, i++) {
   double x = bnd->x();
   double y = bnd->y();
   double z = bnd->z();
   double r = std::sqrt((x*x) + (y*y));
   if(i == 0) {
   minR = r;
   maxR = r;
   zAtMinR = z;
   zAtMaxR = z;
   } else {
   if(minR > r) {
   minR = r;
   zAtMinR = z;
   }
   if(maxR < r) {
   maxR = r; 
   zAtMaxR = z;
   }
   }
   }//end for bnd

   printf("Radius Min = %.15e Max = %.15e \n", minR, maxR);
   printf("Z at min R = %.15e \n",zAtMinR);
   printf("Z at max R = %.15e \n",zAtMaxR);

   }
   */

void myTest(AMP::UnitTest *ut, std::string exeName) {
  std::string input_file = "input_" + exeName;
  std::string log_file = "output_" + exeName;


  AMP::PIO::logOnlyNodeZero(log_file);

  boost::shared_ptr<AMP::InputDatabase> input_db(new AMP::InputDatabase("input_db"));
  AMP::InputManager::getManager()->parseInputFile(input_file, input_db);
  input_db->printClassData(AMP::plog);

  const double contactPrecision = input_db->getDouble("ContactPrecision");
  const std::string masterMeshName = input_db->getString("MasterMesh");
  const std::string slaveMeshName = input_db->getString("SlaveMesh");
  const unsigned int slaveId = input_db->getInteger("SlaveId");
  const unsigned int masterId = input_db->getInteger("MasterId");
  const unsigned int rgDim = input_db->getInteger("Num1Dcells");
  const bool useMLpc = input_db->getBool("USE_ML_PC");
  const unsigned int maxPicardIterations = input_db->getInteger("MaxPicardIterations");
  const unsigned int maxCGiterations = input_db->getInteger("MaxCGiterations");
  const double resTol = input_db->getDouble("tolerance");
  const double initialTemperature = input_db->getDouble("InitialTemperature");
  const double maxFinalTemperature = input_db->getDouble("MaxFinalTemperature");
  const double parabolaConstant = input_db->getDouble("ParabolaConstant"); 

  boost::shared_ptr<AMP::Database> slaveMeshDatabase = input_db->getDatabase("Mesh_2");
  const double slaveXoffset = slaveMeshDatabase->getDouble("x_offset");
  const double slaveYoffset = slaveMeshDatabase->getDouble("y_offset");

  AMP::Mesh::MeshManagerParameters::shared_ptr meshmgrParams ( new AMP::Mesh::MeshManagerParameters ( input_db ) );
  AMP::Mesh::MeshManager::shared_ptr manager ( new AMP::Mesh::MeshManager ( meshmgrParams ) );
  obj.masterMeshAdapter = manager->getMesh ( masterMeshName );
  obj.slaveMeshAdapter = manager->getMesh ( slaveMeshName );

  std::cout<<"Master volume has "<<(obj.masterMeshAdapter->numLocalNodes())<<" nodes."<<std::endl;
  std::cout<<"Slave volume has "<<(obj.slaveMeshAdapter->numLocalNodes())<<" nodes."<<std::endl;

  //std::cout<<"Master:"<<std::endl;
  //checkRadius(masterId, obj.masterMeshAdapter);

  //std::cout<<"Slave:"<<std::endl;
  //checkRadius(slaveId, obj.slaveMeshAdapter);

  boost::shared_ptr<AMP::Operator::ElementPhysicsModel> masterElementPhysicsModel;
  obj.masterOperator = boost::dynamic_pointer_cast<
    AMP::Operator::MechanicsLinearFEOperator>(AMP::Operator::OperatorBuilder::createOperator(
          obj.masterMeshAdapter, "MasterVolumeOperator", input_db, masterElementPhysicsModel));

  boost::shared_ptr<AMP::Operator::ElementPhysicsModel> slaveElementPhysicsModel;
  obj.slaveOperator = boost::dynamic_pointer_cast<
    AMP::Operator::MechanicsLinearFEOperator>(AMP::Operator::OperatorBuilder::createOperator(
          obj.slaveMeshAdapter, "SlaveVolumeOperator", input_db, slaveElementPhysicsModel));

  boost::shared_ptr<AMP::Operator::OperatorParameters> dummyParams;
  obj.columnOperator.reset(new AMP::Operator::ColumnOperator(dummyParams));
  obj.columnOperator->append(obj.masterOperator);
  obj.columnOperator->append(obj.slaveOperator);

  obj.masterVar = obj.masterOperator->getOutputVariable();
  obj.slaveVar = obj.slaveOperator->getOutputVariable();
  AMP::LinearAlgebra::Variable::shared_ptr columnVar = obj.columnOperator->getOutputVariable();

  AMP::LinearAlgebra::Vector::shared_ptr nullVec;
  AMP::LinearAlgebra::Vector::shared_ptr columnSolVec = manager->createVector(columnVar);
  AMP::LinearAlgebra::Vector::shared_ptr columnRhsVec = columnSolVec->cloneVector();
  AMP::LinearAlgebra::Vector::shared_ptr columnResVec = columnSolVec->cloneVector();

  AMP::LinearAlgebra::Vector::shared_ptr slaveSolVec = columnSolVec->subsetVectorForVariable(obj.slaveVar);
  AMP::LinearAlgebra::Vector::shared_ptr slaveRhsVec = columnRhsVec->subsetVectorForVariable(obj.slaveVar);
  AMP::LinearAlgebra::Vector::shared_ptr slaveOffsetCorrectionVec = slaveRhsVec->cloneVector();
  AMP::LinearAlgebra::Vector::shared_ptr masterSolVec = columnSolVec->subsetVectorForVariable(obj.masterVar);

  boost::shared_ptr<AMP::Database> masterBndDatabase = input_db->getDatabase("MasterDirichletBoundary");
  obj.numMasterDirichletDofs = 0;
  obj.masterBndIds.resize(masterBndDatabase->getInteger("number_of_ids"));
  obj.masterBndDofIds.resize(obj.masterBndIds.size());
  obj.masterBndDofVals.resize(obj.masterBndIds.size());
  for(unsigned int i = 0; i < obj.masterBndIds.size(); i++) {
    char tmp[200];
    sprintf(tmp, "id_%u", i);
    obj.masterBndIds[i] = masterBndDatabase->getInteger(tmp);
    sprintf(tmp, "number_of_dofs_%u", i); 
    obj.masterBndDofIds[i].resize(masterBndDatabase->getInteger(tmp));
    obj.masterBndDofVals[i].resize(obj.masterBndDofIds[i].size());
    AMP::Mesh::MeshManager::Adapter::OwnedBoundaryNodeIterator bnd =
      obj.masterMeshAdapter->beginOwnedBoundary( obj.masterBndIds[i] );
    AMP::Mesh::MeshManager::Adapter::OwnedBoundaryNodeIterator end_bnd = 
      obj.masterMeshAdapter->endOwnedBoundary( obj.masterBndIds[i] );
    for(; bnd != end_bnd; bnd++) {
      obj.numMasterDirichletDofs += obj.masterBndDofIds[i].size();
    }//end for bnd
    for(unsigned int j = 0; j < obj.masterBndDofIds[i].size(); j++) {
      char key[200];
      sprintf(key, "dof_%u_%u", i, j);
      obj.masterBndDofIds[i][j] = masterBndDatabase->getInteger(key);
      sprintf(key, "value_%u_%u", i, j);
      obj.masterBndDofVals[i][j] = masterBndDatabase->getDouble(key);
    }//end for j
  }//end for i

  boost::shared_ptr<AMP::Database> temperatureRhsDatabase = input_db->getDatabase("TemperatureRHS");

  AMP::LinearAlgebra::Variable::shared_ptr tempVar(new AMP::LinearAlgebra::VectorVariable<
      AMP::Mesh::NodalVariable, 1>("temp", obj.slaveMeshAdapter) );

  AMP::LinearAlgebra::Vector::shared_ptr currTempVec = (obj.slaveMeshAdapter)->createVector( tempVar );
  AMP::LinearAlgebra::Vector::shared_ptr prevTempVec = (obj.slaveMeshAdapter)->createVector( tempVar );

  prevTempVec->setToScalar(initialTemperature);

  AMP::Mesh::MeshManager::Adapter::OwnedNodeIterator nd = obj.slaveMeshAdapter->beginOwnedNode();
  AMP::Mesh::MeshManager::Adapter::OwnedNodeIterator end_nd = obj.slaveMeshAdapter->endOwnedNode();

  for( ; nd != end_nd; ++nd) {
    double xPos = nd->x();
    double yPos = nd->y(); 
    double xDelta = xPos - slaveXoffset;
    double yDelta = yPos - slaveYoffset;
    double val = maxFinalTemperature - (parabolaConstant*((xDelta*xDelta) + (yDelta*yDelta)));
    currTempVec->setValueByGlobalID(nd->globalID(), val);
  }//end for nd
  currTempVec->makeConsistent(AMP::LinearAlgebra::Vector::CONSISTENT_SET);

  if(currTempVec->min() < initialTemperature) {
    assert(false);
  }

  columnSolVec->zero();

  slaveOffsetCorrectionVec->zero();

  for(unsigned int picardIter = 0; picardIter < maxPicardIterations; picardIter++) {
    std::cout<<"Picard-Iter = "<<picardIter<<std::endl;

    std::cout<<"# slave contact nodes = "<<(obj.slaveNodes.size())<<std::endl;

    obj.numSlaveDirichletDofs = 0;
    if(obj.slaveNodes.size() > 2) {
      obj.slaveBndIds.clear();
      obj.slaveBndDofIds.clear();
      obj.slaveBndDofVals.clear();
    } else {
      boost::shared_ptr<AMP::Database> slaveBndDatabase = input_db->getDatabase("SlaveDirichletBoundary");
      obj.slaveBndIds.resize(slaveBndDatabase->getInteger("number_of_ids"));
      obj.slaveBndDofIds.resize(obj.slaveBndIds.size());
      obj.slaveBndDofVals.resize(obj.slaveBndIds.size());
      for(unsigned int i = 0; i < obj.slaveBndIds.size(); i++) {
        char tmp[200];
        sprintf(tmp, "id_%u", i);
        obj.slaveBndIds[i] = slaveBndDatabase->getInteger(tmp);
        sprintf(tmp, "number_of_dofs_%u", i); 
        obj.slaveBndDofIds[i].resize(slaveBndDatabase->getInteger(tmp));
        obj.slaveBndDofVals[i].resize(obj.slaveBndDofIds[i].size());
        AMP::Mesh::MeshManager::Adapter::OwnedBoundaryNodeIterator bnd =
          obj.slaveMeshAdapter->beginOwnedBoundary( obj.slaveBndIds[i] );
        AMP::Mesh::MeshManager::Adapter::OwnedBoundaryNodeIterator end_bnd = 
          obj.slaveMeshAdapter->endOwnedBoundary( obj.slaveBndIds[i] );
        for(; bnd != end_bnd; bnd++) {
          obj.numSlaveDirichletDofs += obj.slaveBndDofIds[i].size();
        }//end for bnd
        for(unsigned int j = 0; j < obj.slaveBndDofIds[i].size(); j++) {
          char key[200];
          sprintf(key, "dof_%u_%u", i, j);
          obj.slaveBndDofIds[i][j] = slaveBndDatabase->getInteger(key);
          sprintf(key, "value_%u_%u", i, j);
          obj.slaveBndDofVals[i][j] = slaveBndDatabase->getDouble(key);
        }//end for j
      }//end for i
    }

    if(useMLpc) {
      int maxMLiters = 1;
      int printInfoLevel = 10;
      int numGrids = 10;
      int numPDEs = 3;
      int maxCoarseSize = 128;

      computeMl2MasterOrSlave();
      obj.mlSize = obj.ml2MasterOrSlave.size();

      computeMaster2Ml();
      computeSlave2Ml();
      computeSlaveId2SlaveListIndex();
      computeMaster2Slaves();

      obj.matVecInVec = columnSolVec->cloneVector();
      obj.matVecOutVec = columnSolVec->cloneVector();
      obj.pcSolVec.resize(obj.mlSize);
      obj.pcRhsVec.resize(obj.mlSize);

      ML_Create (&ml_object, numGrids);
      ML_Init_Amatrix(ml_object, 0, obj.mlSize, obj.mlSize, NULL);
      ML_Set_Amatrix_Getrow(ml_object, 0, &myGetMultipleRows, NULL, obj.mlSize);
      ML_Set_Amatrix_Matvec(ml_object, 0, &myMatVec);
      ML_Set_MaxIterations(ml_object, maxMLiters);
      ML_Set_ResidualOutputFrequency(ml_object, 0);
      ML_Set_PrintLevel( printInfoLevel );
      ML_Set_OutputLevel(ml_object, printInfoLevel );

      ML_Aggregate_Create(&agg_object);
      agg_object->num_PDE_eqns = numPDEs;
      agg_object->nullspace_dim = numPDEs;
      ML_Aggregate_Set_MaxCoarseSize(agg_object, maxCoarseSize);
      ML_Aggregate_Set_CoarsenScheme_UncoupledMIS(agg_object);

      int nlevels = ML_Gen_MGHierarchy_UsingAggregation(ml_object, 0, ML_INCREASING, agg_object);
      std::cout<<"Number of actual levels : "<< nlevels <<std::endl;

      for(int lev = 0; lev < (nlevels - 1); lev++) {
        ML_Gen_Smoother_SymGaussSeidel(ml_object, lev, ML_BOTH, 2, 1.0);
      }
      ML_Gen_Smoother_Amesos(ml_object, (nlevels - 1), ML_AMESOS_KLU, -1, 0.0);

      ML_Gen_Solver(ml_object, ML_MGV, 0, (nlevels-1));
    }

    columnRhsVec->zero();
    computeTemperatureRhsVector(obj.slaveMeshAdapter, temperatureRhsDatabase, tempVar,
        obj.slaveVar, currTempVec, prevTempVec, columnRhsVec);

    /*
       columnSolVec->zero();
       setMasterDirichletVals(columnSolVec);
       setSlaveDirichletVals(columnSolVec);
       masterToSlaveCorrection(columnSolVec);
       obj.columnOperator->apply(nullVec, columnSolVec, columnResVec, 1.0, 0.0);
       columnRhsVec->subtract(columnRhsVec, columnResVec);
       */

    slaveRhsVec->subtract(slaveRhsVec, slaveOffsetCorrectionVec);

    {
      std::cout<<"MPC-MLCG: "<<std::endl;
      AMP::LinearAlgebra::Vector::shared_ptr MatOutVec = columnSolVec->cloneVector();
      AMP::LinearAlgebra::Vector::shared_ptr pVec = columnSolVec->cloneVector();
      AMP::LinearAlgebra::Vector::shared_ptr zVec = columnSolVec->cloneVector();

      columnSolVec->zero();

      obj.columnOperator->apply(nullVec, columnSolVec, MatOutVec, 1.0, 0.0);
      columnResVec->subtract(columnRhsVec, MatOutVec);
      slaveToMasterCorrection(columnResVec);
      setSlaveToZero(columnResVec);
      setMasterDirichletDofsToZero(columnResVec);
      setSlaveDirichletDofsToZero(columnResVec);

      if(useMLpc) {
        pcMLapply(columnResVec, zVec);
      } else {
        pcNoneApply(columnResVec, zVec);
      }
      setMasterDirichletDofsToZero(zVec);
      setSlaveDirichletDofsToZero(zVec);
      masterToSlaveCorrection(zVec);

      pVec->copyVector(zVec);

      bool cgConverged = false;
      for(unsigned int cgIter = 0; cgIter < maxCGiterations; cgIter++) {
        double resNorm = columnResVec->L2Norm();

        // std::cout<<"CG-Iter = "<<cgIter<<" ResNorm2 = "<<std::setprecision(15)<<resNorm<<std::endl;

        if(resNorm <= resTol) {
          std::cout<<"CG converged to tolerance = "<<std::setprecision(15)
            <<resTol<<" in "<<cgIter<<" iterations."<<std::endl;
          cgConverged = true;
          break;
        }

        obj.columnOperator->apply(nullVec, pVec, MatOutVec, 1.0, 0.0);
        slaveToMasterCorrection(MatOutVec);
        setSlaveToZero(MatOutVec);
        setMasterDirichletDofsToZero(MatOutVec);
        setSlaveDirichletDofsToZero(MatOutVec);

        double resOldDotZ = columnResVec->dot(zVec);

        double alphaDenom = MatOutVec->dot(pVec);

        double alpha = resOldDotZ/alphaDenom;

        columnSolVec->axpy(alpha, pVec, columnSolVec);

        columnResVec->axpy(-alpha, MatOutVec, columnResVec);

        if(useMLpc) {
          pcMLapply(columnResVec, zVec);
        } else {
          pcNoneApply(columnResVec, zVec);
        }
        setMasterDirichletDofsToZero(zVec);
        setSlaveDirichletDofsToZero(zVec);
        masterToSlaveCorrection(zVec);

        double resNewDotZ = columnResVec->dot(zVec);

        double beta = resNewDotZ/resOldDotZ;

        pVec->axpy(beta, pVec, zVec);
      }//end for cgIter
      assert(cgConverged);
    }

    //setMasterDirichletVals(columnSolVec);
    //setSlaveDirichletVals(columnSolVec);
    masterToSlaveCorrection(columnSolVec);
    addSlaveOffsets(columnSolVec);

    if(useMLpc) {
      ML_Aggregate_Destroy(&agg_object);
      ML_Destroy(&ml_object);
    }

    obj.slaveMeshAdapter->displaceMesh(slaveSolVec);
    obj.masterMeshAdapter->displaceMesh(masterSolVec);

    //std::cout<<"Master:"<<std::endl;
    //checkRadius(masterId, obj.masterMeshAdapter);

    //std::cout<<"Slave:"<<std::endl;
    //checkRadius(slaveId, obj.slaveMeshAdapter);

    double minXYZ[3];
    double maxXYZ[3];

    computeRGboundingBox(contactPrecision, obj.masterMeshAdapter, minXYZ, maxXYZ);

    double rgH[3];
    std::vector<std::vector<size_t> > rg2ElemMap;

    computeRG2ElemMap(contactPrecision, rgDim, obj.masterMeshAdapter, minXYZ, maxXYZ, rg2ElemMap, rgH);

    std::vector<size_t> slave2MasterElem;

    std::vector<size_t> tmpSlaveNodes;
    std::vector<std::vector<size_t> > tmpSlave2MasterNodes;
    std::vector<std::vector<double> > tmpSlave2MasterFactors;

    computeSlave2MasterElem(slaveId, masterId, contactPrecision, rgH, rgDim, 
        obj.slaveMeshAdapter, obj.masterMeshAdapter, minXYZ,
        maxXYZ, rg2ElemMap, tmpSlaveNodes, slave2MasterElem);

    computeSlave2MasterNodes(contactPrecision, slaveId, masterId, 
        obj.slaveMeshAdapter, obj.masterMeshAdapter,
        tmpSlaveNodes, slave2MasterElem, tmpSlave2MasterNodes, tmpSlave2MasterFactors);

    columnSolVec->scale(-1.0);
    obj.slaveMeshAdapter->displaceMesh(slaveSolVec);
    obj.masterMeshAdapter->displaceMesh(masterSolVec);
    columnSolVec->scale(-1.0);

    bool activeSetChanged = true;
    if(tmpSlaveNodes == obj.slaveNodes) {
      std::cout<<"Slave nodes has not changed!"<<std::endl;
      if(tmpSlave2MasterNodes == obj.slave2MasterNodes) {
        std::cout<<"Slave2MasterNodes has not changed!"<<std::endl;
        bool factorsEquals = true;
        for(size_t i = 0; i < tmpSlave2MasterFactors.size(); i++) {
          for(size_t j = 0; j < tmpSlave2MasterFactors[i].size(); j++) {
            if(fabs(tmpSlave2MasterFactors[i][j] - obj.slave2MasterFactors[i][j]) > contactPrecision) {
              factorsEquals = false;
              break;
            }
          }//end for j
          if(!factorsEquals) {
            break;
          }
        }//end for i
        if(factorsEquals) {
          activeSetChanged = false;
        }
      }
    }

    if(activeSetChanged) {
      std::cout<<"Active set changed!"<<std::endl;
      obj.slaveNodes = tmpSlaveNodes;
      obj.slave2MasterNodes = tmpSlave2MasterNodes;
      obj.slave2MasterFactors = tmpSlave2MasterFactors;
    } else {
      std::cout<<"Picard converged in "<<picardIter<<" iterations."<<std::endl;
      break;
    }

    computeSlaveOffsets();

    slaveSolVec->zero();
    addSlaveOffsets(slaveSolVec);
    obj.slaveOperator->apply(nullVec, slaveSolVec, slaveOffsetCorrectionVec, 1.0, 0.0);

  }//end for picardIter

#ifdef USE_SILO
  manager->registerVectorAsData ( columnSolVec, "Displacement" );
  manager->registerVectorAsData ( currTempVec, "FinalTemperature" );
  manager->writeFile<AMP::Mesh::SiloIO> ( exeName , 1 );

  deformMesh(obj.masterMeshAdapter, masterSolVec);
  deformMesh(obj.slaveMeshAdapter, slaveSolVec);

  manager->writeFile<AMP::Mesh::SiloIO> ( exeName + "_Deformed" , 1 );
#endif

  freeMyData();

  ut->passes(exeName);
}

int main(int argc, char *argv[])
{
  AMP::AMPManager::startup(argc, argv);
  AMP::UnitTest ut;

  std::string exeName = "testContactV6";

  try {
    myTest(&ut, exeName);
  } catch (std::exception &err) {
    std::cout << "ERROR: While testing "<<argv[0] << err.what() << std::endl;
    ut.failure("ERROR: While testing");
  } catch( ... ) {
    std::cout << "ERROR: While testing "<<argv[0] << "An unknown exception was thrown." << std::endl;
    ut.failure("ERROR: While testing");
  }

  ut.report();

  int num_failed = ut.NumFailGlobal();
  AMP::AMPManager::shutdown();
  return num_failed;
}  


