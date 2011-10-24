
#include "utils/InputManager.h"
#include "utils/AMPManager.h"
#include "utils/UnitTest.h"
#include "utils/Utilities.h"

#include "boost/shared_ptr.hpp"

#include <iostream>
#include <string>
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

#include "fe_interface.h"
#include "cell_hex8.h"

#include "ml_include.h"
#include "solvers/MLoptions.h"
#include "indexHolder.h"

#include "ContactSearchUtils.h"

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
  std::vector<size_t> slaveNodes;
  std::vector<std::vector<size_t> > slave2MasterNodes;
  std::vector<std::vector<double> > slave2MasterFactors;
  std::vector<unsigned int> bndIds;
  std::vector<std::vector<unsigned int> > bndDofIds;
  std::vector<std::vector<double> > bndDofVals;
  size_t mlSize;
  size_t totalDirichletDofs;
} obj;

void freeMyData() {
  obj.columnOperator.reset();
  obj.slaveMeshAdapter.reset();
  obj.masterMeshAdapter.reset();
  obj.masterOperator.reset();
  obj.slaveOperator.reset();
  obj.slaveVar.reset();
  obj.masterVar.reset();
  obj.matVecInVec.reset();
  obj.matVecOutVec.reset();
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

void setDirichletVals(AMP::LinearAlgebra::Vector::shared_ptr vec) {
  AMP::Mesh::DOFMap::shared_ptr master_dof_map = obj.masterMeshAdapter->getDOFMap(obj.masterVar);

  AMP::LinearAlgebra::Vector::shared_ptr masterVec =  vec->subsetVectorForVariable(obj.masterVar);

  for(size_t i = 0; i < obj.bndIds.size(); i++) {
    AMP::Mesh::MeshManager::Adapter::OwnedBoundaryNodeIterator bnd = 
      obj.masterMeshAdapter->beginOwnedBoundary( obj.bndIds[i] );
    AMP::Mesh::MeshManager::Adapter::OwnedBoundaryNodeIterator end_bnd = 
      obj.masterMeshAdapter->endOwnedBoundary( obj.bndIds[i] );

    for( ; bnd != end_bnd; ++bnd) {
      std::vector<unsigned int> bndGlobalIds;
      master_dof_map->getDOFs(*bnd, bndGlobalIds, obj.bndDofIds[i]);

      for(size_t j = 0; j < bndGlobalIds.size(); j++) {
        masterVec->setLocalValueByGlobalID(bndGlobalIds[j], obj.bndDofVals[i][j]);
      }//end for j
    }//end for bnd
  }//end for i

  vec->makeConsistent( AMP::LinearAlgebra::Vector::CONSISTENT_SET );
}

void setDirichletDofsToZero(AMP::LinearAlgebra::Vector::shared_ptr vec) {
  AMP::Mesh::DOFMap::shared_ptr master_dof_map = obj.masterMeshAdapter->getDOFMap(obj.masterVar);

  AMP::LinearAlgebra::Vector::shared_ptr masterVec =  vec->subsetVectorForVariable(obj.masterVar);

  for(size_t i = 0; i < obj.bndIds.size(); i++) {
    AMP::Mesh::MeshManager::Adapter::OwnedBoundaryNodeIterator bnd = 
      obj.masterMeshAdapter->beginOwnedBoundary( obj.bndIds[i] );
    AMP::Mesh::MeshManager::Adapter::OwnedBoundaryNodeIterator end_bnd = 
      obj.masterMeshAdapter->endOwnedBoundary( obj.bndIds[i] );

    for( ; bnd != end_bnd; ++bnd) {
      std::vector<unsigned int> bndGlobalIds;
      master_dof_map->getDOFs(*bnd, bndGlobalIds, obj.bndDofIds[i]);

      for(size_t j = 0; j < bndGlobalIds.size(); j++) {
        masterVec->setLocalValueByGlobalID(bndGlobalIds[j], 0.0);
      }//end for j
    }//end for bnd
  }//end for i

  vec->makeConsistent( AMP::LinearAlgebra::Vector::CONSISTENT_SET );
}

void copyML2Column(double* mlVec, AMP::LinearAlgebra::Vector::shared_ptr vec) {
  AMP::LinearAlgebra::Vector::shared_ptr slaveVec = vec->subsetVectorForVariable(obj.slaveVar);
  AMP::LinearAlgebra::Vector::shared_ptr masterVec = vec->subsetVectorForVariable(obj.masterVar);

  size_t totalMasterDofs = (3*(obj.masterMeshAdapter->numLocalNodes()));
  size_t numMlMasterDofs = totalMasterDofs - obj.totalDirichletDofs;

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
  size_t numMlMasterDofs = totalMasterDofs - obj.totalDirichletDofs;

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

  setDirichletDofsToZero(obj.matVecInVec);

  masterToSlaveCorrection(obj.matVecInVec);

  (obj.columnOperator)->apply(nullVec, obj.matVecInVec, obj.matVecOutVec, 1.0, 0.0);

  slaveToMasterCorrection(obj.matVecOutVec);

  copyColumn2ML(out, obj.matVecOutVec);

  return 0;
}

void myGetSingleRow(unsigned int row, std::vector<unsigned int> &cols, std::vector<double> &values) {
  size_t totalMasterDofs = (3*(obj.masterMeshAdapter->numLocalNodes()));
  size_t numMlMasterDofs = totalMasterDofs - obj.totalDirichletDofs;

  cols.clear();
  values.clear();

  if(row < numMlMasterDofs) {
    std::vector<unsigned int> masterCols;
    std::vector<double> masterVals;
    unsigned int masterRow = obj.ml2MasterOrSlave[row];
    AMP::LinearAlgebra::Matrix::shared_ptr masterMat = (obj.masterOperator)->getMatrix();
    masterMat->getRowByGlobalID(masterRow, masterCols, masterVals);

    for(std::vector<unsigned int>::iterator colIter = masterCols.begin(); colIter != masterCols.end();) {
      if(obj.master2Ml[*colIter] == static_cast<unsigned int>(-1)) {
        std::vector<double>::iterator valIter = masterVals.begin() + (colIter - masterCols.begin());
        colIter = masterCols.erase(colIter);
        valIter = masterVals.erase(valIter);
      } else {
        colIter++;
      }//end if dirichlet
    }//end for

    if(obj.master2SlaveIds[masterRow].empty()) {
      for(size_t i = 0; i < masterCols.size(); i++) {
        cols.push_back(obj.master2Ml[masterCols[i]]);
        values.push_back(masterVals[i]);
      }//end for i
    } else {
      std::vector<IndexHolder<unsigned int> > ihVec;
      std::vector<double> tmpVals;
      for(size_t i = 0; i < masterCols.size(); i++) {
        if(obj.master2SlaveIds[masterCols[i]].empty()) {
          cols.push_back(obj.master2Ml[masterCols[i]]);
          values.push_back(masterVals[i]);
        } else {
          IndexHolder<unsigned int> tmpElem;
          tmpElem.value = &(obj.master2Ml[masterCols[i]]);
          tmpElem.index = tmpVals.size();
          tmpVals.push_back(masterVals[i]);
          ihVec.push_back(tmpElem);
        }//end if contact col
      }//end for i

      AMP::Mesh::DOFMap::shared_ptr master_dof_map = obj.masterMeshAdapter->getDOFMap(obj.masterVar);

      AMP::LinearAlgebra::Matrix::shared_ptr slaveMat = (obj.slaveOperator)->getMatrix();

      for(size_t i = 0; i < obj.master2SlaveIds[masterRow].size(); i++) {
        std::vector<unsigned int> slaveCols;
        std::vector<double> slaveVals;
        unsigned int slaveRow = obj.master2SlaveIds[masterRow][i];
        slaveMat->getRowByGlobalID(slaveRow, slaveCols, slaveVals);
        for(size_t j = 0; j < slaveCols.size(); j++) {
          if(obj.slave2Ml[slaveCols[j]] != static_cast<unsigned int>(-1)) {
            IndexHolder<unsigned int> tmpElem;
            tmpElem.value = &(obj.slave2Ml[slaveCols[j]]);
            tmpElem.index = tmpVals.size();
            tmpVals.push_back(obj.master2SlaveFactors[masterRow][i]*slaveVals[j]);
            ihVec.push_back(tmpElem);
          } else {
            unsigned int slaveIndex = obj.slaveId2SlaveListIndex[slaveCols[j]];
            unsigned int dofId = obj.slaveId2Dof[slaveCols[j]];
            for(size_t k = 0; k < obj.slave2MasterNodes[slaveIndex].size(); k++) {
              AMP::Mesh::LibMeshNode masterNd =  obj.masterMeshAdapter->getNode(obj.slave2MasterNodes[slaveIndex][k]);
              std::vector<unsigned int> masterGlobalIds;
              std::vector<unsigned int> singleton(1);
              singleton[0] = dofId;
              master_dof_map->getDOFs(masterNd, masterGlobalIds, singleton);
              IndexHolder<unsigned int> tmpElem;
              tmpElem.value = &(obj.master2Ml[masterGlobalIds[0]]);
              if((*(tmpElem.value)) != static_cast<unsigned int>(-1)) {
                tmpElem.index = tmpVals.size();
                tmpVals.push_back(obj.master2SlaveFactors[masterRow][i]*
                    slaveVals[j]*obj.slave2MasterFactors[slaveIndex][k]);
                ihVec.push_back(tmpElem);
              }//end if dirichlet
            }//end for k
          }//end if contact col
        }//end for j
      }//end for i

      assert(!(ihVec.empty()));
      std::sort(ihVec.begin(), ihVec.end());
      cols.push_back(*(ihVec[0].value));
      values.push_back(tmpVals[ihVec[0].index]);
      for(size_t i = 1; i < ihVec.size(); i++) {
        if(ihVec[i] == ihVec[i - 1]) {
          values[values.size() - 1] += tmpVals[ihVec[i].index];
        } else {
          cols.push_back(*(ihVec[i].value));
          values.push_back(tmpVals[ihVec[i].index]);
        }
      }//end for i
    }//end if contact row
  } else {
    std::vector<unsigned int> slaveCols;
    std::vector<double> slaveVals;
    unsigned int slaveRow = obj.ml2MasterOrSlave[row];
    AMP::LinearAlgebra::Matrix::shared_ptr slaveMat = (obj.slaveOperator)->getMatrix();
    slaveMat->getRowByGlobalID(slaveRow, slaveCols, slaveVals);
    AMP::Mesh::DOFMap::shared_ptr master_dof_map = obj.masterMeshAdapter->getDOFMap(obj.masterVar);
    std::vector<IndexHolder<unsigned int> > ihVec;
    std::vector<double> tmpVals;
    for(size_t i = 0; i < slaveCols.size(); i++) {
      if(obj.slave2Ml[slaveCols[i]] != static_cast<unsigned int>(-1)) {
        cols.push_back(obj.slave2Ml[slaveCols[i]]);
        values.push_back(slaveVals[i]);
      } else {
        unsigned int slaveIndex = obj.slaveId2SlaveListIndex[slaveCols[i]];
        unsigned int dofId = obj.slaveId2Dof[slaveCols[i]];
        for(size_t j = 0; j < obj.slave2MasterNodes[slaveIndex].size(); j++) {
          AMP::Mesh::LibMeshNode masterNd =  obj.masterMeshAdapter->getNode(obj.slave2MasterNodes[slaveIndex][j]);
          std::vector<unsigned int> masterGlobalIds;
          std::vector<unsigned int> singleton(1);
          singleton[0] = dofId;
          master_dof_map->getDOFs(masterNd, masterGlobalIds, singleton);
          if(obj.master2Ml[masterGlobalIds[0]] != static_cast<unsigned int>(-1)) {
            IndexHolder<unsigned int> tmpElem;
            tmpElem.value = &(obj.master2Ml[masterGlobalIds[0]]);
            tmpElem.index = tmpVals.size();
            tmpVals.push_back(obj.slave2MasterFactors[slaveIndex][j]*slaveVals[i]);
            ihVec.push_back(tmpElem);
          }//end if dirichlet
        }//end for j
      }//end if contact col
    }//end for i

    if(!(ihVec.empty())) {
      std::sort(ihVec.begin(), ihVec.end());
      cols.push_back(*(ihVec[0].value));
      values.push_back(tmpVals[ihVec[0].index]);
      for(size_t i = 1; i < ihVec.size(); i++) {
        if(ihVec[i] == ihVec[i - 1]) {
          values[values.size() - 1] += tmpVals[ihVec[i].index];
        } else {
          cols.push_back(*(ihVec[i].value));
          values.push_back(tmpVals[ihVec[i].index]);
        }
      }//end for i
    } 
  }//end if master or slave
}

int myGetMultipleRows(ML_Operator *data, int N_requested_rows, int requested_rows[],
    int allocated_space, int columns[], double values[], int row_lengths[] ) {
  int spaceRequired = 0;
  int cnt = 0;
  for(int i = 0; i < N_requested_rows; i++) {
    unsigned int row = requested_rows[i];
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

void computeMaster2Ml() {
  size_t totalMasterDofs = (3*(obj.masterMeshAdapter->numLocalNodes()));

  obj.master2Ml.resize(totalMasterDofs);

  for(size_t i = 0; i < totalMasterDofs; i++) {
    obj.master2Ml[i] = static_cast<unsigned int>(-1);
  }//end for i

  size_t numMlMasterDofs = totalMasterDofs - obj.totalDirichletDofs;

  for(size_t i = 0; i < numMlMasterDofs; i++) {
    obj.master2Ml[obj.ml2MasterOrSlave[i]] = i;
  }//end for i
}

void computeSlave2Ml() {
  size_t totalMasterDofs = (3*(obj.masterMeshAdapter->numLocalNodes()));
  size_t totalSlaveDofs = (3*(obj.slaveMeshAdapter->numLocalNodes()));

  obj.slave2Ml.resize(totalSlaveDofs);

  for(size_t i = 0; i < totalSlaveDofs; i++) {
    obj.slave2Ml[i] = static_cast<unsigned int>(-1);
  }//end for i

  size_t numMlMasterDofs = totalMasterDofs - obj.totalDirichletDofs;

  for(size_t i = numMlMasterDofs; i < obj.mlSize; i++) {
    obj.slave2Ml[obj.ml2MasterOrSlave[i]] = i;
  }//end for i
}

void computeMl2MasterOrSlave() {
  std::vector<unsigned int> tmpVec;

  size_t totalMasterDofs = (3*(obj.masterMeshAdapter->numLocalNodes()));
  size_t totalSlaveDofs = (3*(obj.slaveMeshAdapter->numLocalNodes()));

  tmpVec.resize(totalMasterDofs + totalSlaveDofs);

  size_t cnt = 0;
  for(size_t i = 0; i < totalMasterDofs; i++) {
    tmpVec[cnt] = i;
    cnt++;
  }//end for i
  for(size_t i = 0; i < totalSlaveDofs; i++) {
    tmpVec[cnt] = i;
    cnt++;
  }//end for i

  AMP::Mesh::DOFMap::shared_ptr master_dof_map = obj.masterMeshAdapter->getDOFMap(obj.masterVar);

  for(size_t i = 0; i < obj.bndIds.size(); i++) {
    AMP::Mesh::MeshManager::Adapter::OwnedBoundaryNodeIterator bnd = 
      obj.masterMeshAdapter->beginOwnedBoundary( obj.bndIds[i] );
    AMP::Mesh::MeshManager::Adapter::OwnedBoundaryNodeIterator end_bnd = 
      obj.masterMeshAdapter->endOwnedBoundary( obj.bndIds[i] );

    for( ; bnd != end_bnd; ++bnd) {
      std::vector<unsigned int> bndGlobalIds;
      master_dof_map->getDOFs(*bnd, bndGlobalIds, obj.bndDofIds[i]);
      for(size_t j = 0; j < bndGlobalIds.size(); j++) {
        tmpVec[bndGlobalIds[j]] = static_cast<unsigned int>(-1);
      }//end for j
    }//end for bnd
  }//end for i

  std::vector<unsigned int> dofs(3);
  dofs[0] = 0; dofs[1] = 1; dofs[2] = 2;

  AMP::Mesh::DOFMap::shared_ptr slave_dof_map = obj.slaveMeshAdapter->getDOFMap(obj.slaveVar);

  for(size_t i = 0; i < obj.slaveNodes.size(); i++) {
    AMP::Mesh::LibMeshNode slaveNd =  obj.slaveMeshAdapter->getNode( obj.slaveNodes[i] );
    std::vector<unsigned int> slaveGlobalIds;
    slave_dof_map->getDOFs(slaveNd, slaveGlobalIds, dofs);
    for(size_t j = 0; j < dofs.size(); j++) {
      tmpVec[totalMasterDofs + slaveGlobalIds[j]] = static_cast<unsigned int>(-1);
    }//end for j
  }//end for i

  obj.ml2MasterOrSlave.resize(obj.mlSize);
  cnt = 0;
  for(size_t i = 0; i < tmpVec.size(); i++) {
    if(tmpVec[i] != static_cast<unsigned int>(-1)) {
      obj.ml2MasterOrSlave[cnt] = tmpVec[i];
      cnt++;
    }
  }//end for i
}

void myTest(AMP::UnitTest *ut, std::string exeName) {
  std::string input_file = "input_" + exeName;
  std::string log_file = "output_" + exeName;

  const double precision = 1.0e-12;

  AMP::PIO::logOnlyNodeZero(log_file);

  boost::shared_ptr<AMP::InputDatabase> input_db(new AMP::InputDatabase("input_db"));
  AMP::InputManager::getManager()->parseInputFile(input_file, input_db);
  input_db->printClassData(AMP::plog);

  const std::string masterMeshName = input_db->getString("MasterMesh");
  const std::string slaveMeshName = input_db->getString("SlaveMesh");
  const unsigned int slaveId = input_db->getInteger("SlaveId");
  const unsigned int masterId = input_db->getInteger("MasterId");
  const unsigned int rgDim = input_db->getInteger("Num1Dcells");
  const bool useMLpc = input_db->getBool("USE_ML_PC");

  AMP::Mesh::MeshManagerParameters::shared_ptr meshmgrParams ( new AMP::Mesh::MeshManagerParameters ( input_db ) );
  AMP::Mesh::MeshManager::shared_ptr manager ( new AMP::Mesh::MeshManager ( meshmgrParams ) );
  obj.masterMeshAdapter = manager->getMesh ( masterMeshName );
  obj.slaveMeshAdapter = manager->getMesh ( slaveMeshName );

  std::cout<<"Master volume has "<<(obj.masterMeshAdapter->numLocalNodes())<<" nodes."<<std::endl;
  std::cout<<"Slave volume has "<<(obj.slaveMeshAdapter->numLocalNodes())<<" nodes."<<std::endl;

  double minXYZ[3];
  double maxXYZ[3];

  computeRGboundingBox(precision, obj.masterMeshAdapter, minXYZ, maxXYZ);

  double rgH[3];
  std::vector<std::vector<size_t> > rg2ElemMap;

  computeRG2ElemMap(precision, rgDim, obj.masterMeshAdapter, minXYZ, maxXYZ, rg2ElemMap, rgH);

  std::vector<size_t> slave2MasterElem;

  computeSlave2MasterElem(slaveId, masterId, precision, rgH, rgDim, 
      obj.slaveMeshAdapter, obj.masterMeshAdapter, minXYZ, maxXYZ, rg2ElemMap, obj.slaveNodes, slave2MasterElem);

  std::cout<<"# slave contact nodes = "<<(obj.slaveNodes.size())<<std::endl;

  computeSlave2MasterNodes(precision, slaveId, masterId, obj.slaveMeshAdapter, obj.masterMeshAdapter,
      obj.slaveNodes, slave2MasterElem, obj.slave2MasterNodes, obj.slave2MasterFactors);

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

  boost::shared_ptr<AMP::Operator::ElementPhysicsModel> dummyModel;
  boost::shared_ptr<AMP::Operator::DirichletVectorCorrection> loadOperator = boost::dynamic_pointer_cast<
    AMP::Operator::DirichletVectorCorrection>(AMP::Operator::OperatorBuilder::createOperator(
          obj.slaveMeshAdapter, "SlaveLoadOperator", input_db, dummyModel));
  loadOperator->setVariable(obj.slaveVar);

  AMP::LinearAlgebra::Vector::shared_ptr nullVec;
  AMP::LinearAlgebra::Vector::shared_ptr columnSolVec = manager->createVector(columnVar);
  AMP::LinearAlgebra::Vector::shared_ptr columnRhsVec = columnSolVec->cloneVector();
  AMP::LinearAlgebra::Vector::shared_ptr columnResVec = columnSolVec->cloneVector();

  boost::shared_ptr<AMP::Database> bndDatabase = input_db->getDatabase("MasterDirichletBoundary");
  obj.totalDirichletDofs = 0;
  obj.bndIds.resize(bndDatabase->getInteger("number_of_ids"));
  obj.bndDofIds.resize(obj.bndIds.size());
  obj.bndDofVals.resize(obj.bndIds.size());
  for(unsigned int i = 0; i < obj.bndIds.size(); i++) {
    char tmp[200];
    sprintf(tmp, "id_%u", i);
    obj.bndIds[i] = bndDatabase->getInteger(tmp);
    sprintf(tmp, "number_of_dofs_%u", i); 
    obj.bndDofIds[i].resize(bndDatabase->getInteger(tmp));
    obj.bndDofVals[i].resize(obj.bndDofIds[i].size());
    AMP::Mesh::MeshManager::Adapter::OwnedBoundaryNodeIterator bnd =
      obj.masterMeshAdapter->beginOwnedBoundary( obj.bndIds[i] );
    AMP::Mesh::MeshManager::Adapter::OwnedBoundaryNodeIterator end_bnd = 
      obj.masterMeshAdapter->endOwnedBoundary( obj.bndIds[i] );
    for(; bnd != end_bnd; bnd++) {
      obj.totalDirichletDofs += obj.bndDofIds[i].size();
    }//end for bnd
    for(unsigned int j = 0; j < obj.bndDofIds[i].size(); j++) {
      char key[200];
      sprintf(key, "dof_%u_%u", i, j);
      obj.bndDofIds[i][j] = bndDatabase->getInteger(key);
      sprintf(key, "value_%u_%u", i, j);
      obj.bndDofVals[i][j] = bndDatabase->getDouble(key);
    }//end for j
  }//end for i

  if(useMLpc) {
    int maxMLiters = 1;
    int printInfoLevel = 10;
    int numGrids = 10;
    int numPDEs = 3;
    int maxCoarseSize = 128;
    obj.mlSize = (3*(obj.masterMeshAdapter->numLocalNodes() +
          obj.slaveMeshAdapter->numLocalNodes() - 
          obj.slaveNodes.size())) - obj.totalDirichletDofs;

    computeMl2MasterOrSlave();
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

  loadOperator->apply(nullVec, nullVec, columnRhsVec, 1.0, 0.0);

  columnSolVec->zero();

  setDirichletVals(columnSolVec);

  masterToSlaveCorrection(columnSolVec);

  obj.columnOperator->apply(nullVec, columnSolVec, columnResVec, 1.0, 0.0);

  columnRhsVec->subtract(columnRhsVec, columnResVec);

  int maxCGiters = input_db->getInteger("maxCGiters");

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
    setDirichletDofsToZero(columnResVec);

    if(useMLpc) {
      pcMLapply(columnResVec, zVec);
    } else {
      pcNoneApply(columnResVec, zVec);
    }
    setDirichletDofsToZero(zVec);
    masterToSlaveCorrection(zVec);

    pVec->copyVector(zVec);

    for(int iter = 0; iter <= maxCGiters; iter++) {
      double resNorm = columnResVec->L2Norm();

      std::cout<<"CG-Iter = "<<iter<<" ResNorm2 = "<<std::setprecision(15)<<resNorm<<std::endl;

      obj.columnOperator->apply(nullVec, pVec, MatOutVec, 1.0, 0.0);
      slaveToMasterCorrection(MatOutVec);
      setSlaveToZero(MatOutVec);
      setDirichletDofsToZero(MatOutVec);

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
      setDirichletDofsToZero(zVec);
      masterToSlaveCorrection(zVec);

      double resNewDotZ = columnResVec->dot(zVec);

      double beta = resNewDotZ/resOldDotZ;

      pVec->axpy(beta, pVec, zVec);
    }//end for iter
  }

  setDirichletVals(columnSolVec);

  if(useMLpc) {
    ML_Aggregate_Destroy(&agg_object);
    ML_Destroy(&ml_object);
  }

  freeMyData();

#ifdef USE_SILO
  manager->registerVectorAsData ( columnSolVec, "Displacement" );
  manager->writeFile<AMP::Mesh::SiloIO> ( exeName , 1 );
#endif

  ut->passes(exeName);
}

int main(int argc, char *argv[])
{
  AMP::AMPManager::startup(argc, argv);
  AMP::UnitTest ut;

  std::string exeName = "testContactV5withML";

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


