
#include "utils/InputManager.h"
#include "utils/AMPManager.h"
#include "utils/UnitTest.h"
#include "utils/Utilities.h"
#include "utils/WriteSolutionToFile.h"

#include <iostream>
#include <string>
#include <vector>
#include <algorithm>

#include "ampmesh/SiloIO.h"

#include "operators/ColumnOperator.h"
#include "operators/LinearBVPOperator.h"
#include "operators/OperatorBuilder.h"
#include "operators/boundary/DirichletVectorCorrection.h"

#include "ml_include.h"
#include "solvers/MLoptions.h"

#if 0
#include "TiedContactUtils.h"
#include "MPCtestUtils.h"

ML_Aggregate *agg_object;
ML* ml_object;

int getColCnt = 0;

struct MLdata {
  boost::shared_ptr<AMP::Operator::ColumnOperator> columnOperator;
  boost::shared_ptr<AMP::Operator::LinearBVPOperator> masterOperator;
  boost::shared_ptr<AMP::Operator::LinearBVPOperator> slaveOperator;
  std::vector<unsigned int> masterContactNodes;
  std::vector<unsigned int> slaveContactNodes;
  std::vector<unsigned int> masterVolumeNodes;
  std::vector<unsigned int> slaveVolumeNodes;
  AMP::LinearAlgebra::Variable::shared_ptr masterVar;
  AMP::LinearAlgebra::Variable::shared_ptr slaveVar;
  AMP::LinearAlgebra::Variable::shared_ptr columnVar;
  AMP::Mesh::MeshManager::Adapter::shared_ptr masterMeshAdapter;
  AMP::Mesh::MeshManager::Adapter::shared_ptr slaveMeshAdapter;
  AMP::Mesh::MeshManager::shared_ptr manager;
  AMP::LinearAlgebra::Vector::shared_ptr matVecInVec;
  AMP::LinearAlgebra::Vector::shared_ptr matVecOutVec;
  std::vector<double> getColInVec;
  std::vector<double> getColOutVec;
  std::vector<double> pcSolVec;
  std::vector<double> pcRhsVec;
  std::vector<unsigned int> ml2MasterOrSlave;
  std::vector<unsigned int> master2Ml;
  std::vector<unsigned int> slave2Ml;
  std::vector<unsigned int> masterContact2SlaveContact;
  size_t vecSize;
} mlData;

void freeMLdata() {
  mlData.columnOperator.reset();
  mlData.masterOperator.reset();
  mlData.slaveOperator.reset();
  mlData.masterVar.reset();
  mlData.slaveVar.reset();
  mlData.columnVar.reset();
  mlData.masterMeshAdapter.reset();
  mlData.slaveMeshAdapter.reset();
  mlData.manager.reset();
  mlData.matVecInVec.reset();
  mlData.matVecOutVec.reset();
}

int myMatVec(ML_Operator *data, int in_length, double in[], int out_length, double out[]) {
  AMP::LinearAlgebra::Vector::shared_ptr columnInVec = mlData.matVecInVec;

  AMP::LinearAlgebra::Vector::shared_ptr columnOutVec = mlData.matVecOutVec;

  AMP::LinearAlgebra::Vector::shared_ptr masterInVec =
    columnInVec->subsetVectorForVariable(mlData.masterVar);

  AMP::LinearAlgebra::Vector::shared_ptr slaveInVec =
    columnInVec->subsetVectorForVariable(mlData.slaveVar);

  AMP::LinearAlgebra::Vector::shared_ptr masterOutVec =
    columnOutVec->subsetVectorForVariable(mlData.masterVar);

  AMP::LinearAlgebra::Vector::shared_ptr slaveOutVec =
    columnOutVec->subsetVectorForVariable(mlData.slaveVar);

  AMP::Mesh::DOFMap::shared_ptr master_dof_map = 
    (mlData.masterMeshAdapter)->getDOFMap(mlData.masterVar);

  AMP::Mesh::DOFMap::shared_ptr slave_dof_map = 
    (mlData.slaveMeshAdapter)->getDOFMap(mlData.slaveVar);

  std::vector<unsigned int> dofs(3);
  dofs[0] = 0; dofs[1] = 1; dofs[2] = 2;

  int arrIdx = 0;
  for(size_t i = 0; i < ((mlData.masterVolumeNodes).size()); i++) {
    AMP::Mesh::LibMeshNode masterNd = (mlData.masterMeshAdapter)->getNode( (mlData.masterVolumeNodes)[i] );
    std::vector<unsigned int> masterGlobalIds;
    master_dof_map->getDOFs(masterNd, masterGlobalIds, dofs);
    for(size_t j = 0; j < dofs.size(); j++) {
      masterInVec->setLocalValueByGlobalID(masterGlobalIds[j], in[arrIdx]);
      arrIdx++;
    }
  }

  for(size_t i = 0; i < ((mlData.masterContactNodes).size()); i++) {
    AMP::Mesh::LibMeshNode masterNd = (mlData.masterMeshAdapter)->getNode( (mlData.masterContactNodes)[i] );
    std::vector<unsigned int> masterGlobalIds;
    master_dof_map->getDOFs(masterNd, masterGlobalIds, dofs);
    for(size_t j = 0; j < dofs.size(); j++) {
      masterInVec->setLocalValueByGlobalID(masterGlobalIds[j], in[arrIdx]);
      arrIdx++;
    }
  }

  for(size_t i = 0; i < ((mlData.slaveVolumeNodes).size()); i++) {
    AMP::Mesh::LibMeshNode slaveNd = (mlData.slaveMeshAdapter)->getNode( (mlData.slaveVolumeNodes)[i] );
    std::vector<unsigned int> slaveGlobalIds;
    slave_dof_map->getDOFs(slaveNd, slaveGlobalIds, dofs);
    for(size_t j = 0; j < dofs.size(); j++) {
      slaveInVec->setLocalValueByGlobalID(slaveGlobalIds[j], in[arrIdx]);
      arrIdx++;
    }
  }

  copyMasterToSlave(mlData.slaveVar, mlData.masterVar, columnInVec, mlData.slaveMeshAdapter,
      mlData.masterMeshAdapter, mlData.slaveContactNodes, mlData.masterContactNodes);

  AMP::LinearAlgebra::Vector::shared_ptr nullVec;

  (mlData.columnOperator)->apply(nullVec, columnInVec, columnOutVec, 1.0, 0.0);
  addSlaveToMaster(mlData.slaveVar, mlData.masterVar, columnOutVec, mlData.slaveMeshAdapter,
      mlData.masterMeshAdapter, mlData.slaveContactNodes, mlData.masterContactNodes);

  arrIdx = 0;
  for(size_t i = 0; i < ((mlData.masterVolumeNodes).size()); i++) {
    AMP::Mesh::LibMeshNode masterNd = (mlData.masterMeshAdapter)->getNode( (mlData.masterVolumeNodes)[i] );
    std::vector<unsigned int> masterGlobalIds;
    master_dof_map->getDOFs(masterNd, masterGlobalIds, dofs);
    for(size_t j = 0; j < dofs.size(); j++) {
      out[arrIdx] = masterOutVec->getLocalValueByGlobalID(masterGlobalIds[j]);
      arrIdx++;
    }
  }

  for(size_t i = 0; i < ((mlData.masterContactNodes).size()); i++) {
    AMP::Mesh::LibMeshNode masterNd = (mlData.masterMeshAdapter)->getNode( (mlData.masterContactNodes)[i] );
    std::vector<unsigned int> masterGlobalIds;
    master_dof_map->getDOFs(masterNd, masterGlobalIds, dofs);
    for(size_t j = 0; j < dofs.size(); j++) {
      out[arrIdx] = masterOutVec->getLocalValueByGlobalID(masterGlobalIds[j]);
      arrIdx++;
    }
  }

  for(size_t i = 0; i < ((mlData.slaveVolumeNodes).size()); i++) {
    AMP::Mesh::LibMeshNode slaveNd = (mlData.slaveMeshAdapter)->getNode( (mlData.slaveVolumeNodes)[i] );
    std::vector<unsigned int> slaveGlobalIds;
    slave_dof_map->getDOFs(slaveNd, slaveGlobalIds, dofs);
    for(size_t j = 0; j < dofs.size(); j++) {
      out[arrIdx] = slaveOutVec->getLocalValueByGlobalID(slaveGlobalIds[j]);
      arrIdx++;
    }
  }

  return 0;
}

void myGetColumn(int column, std::vector<unsigned int> &rows, std::vector<double> &values) {
  std::cout<<"GetColCnt = "<<(getColCnt++)<<std::endl;

  double * inVec = &((mlData.getColInVec)[0]); 
  double * outVec = &((mlData.getColOutVec)[0]);

  for(int i = 0; i < (mlData.vecSize); i++) {
    inVec[i] = 0.0;
  }
  inVec[column] = 1.0;

  myMatVec(NULL, (mlData.vecSize), inVec, (mlData.vecSize), outVec);

  rows.clear();
  values.clear();
  for(size_t i = 0; i < (mlData.vecSize); i++) {
    if(outVec[i]) {
      rows.push_back(i);
      values.push_back(outVec[i]);
    }
  }

}

void myGetSingleRow(int row, std::vector<unsigned int> &cols, std::vector<double> &values) {
  cols.clear();
  values.clear();
  if( row < (3*((mlData.masterVolumeNodes).size())) ) {
    //MasterVolume
    AMP::LinearAlgebra::Matrix::shared_ptr masterMat = (mlData.masterOperator)->getMatrix();
    std::vector<unsigned int> masterCols;
    masterMat->getRowByGlobalID(mlData.ml2MasterOrSlave[row], masterCols, values);
    for(size_t i = 0; i < (masterCols.size()); i++) {
      cols.push_back(mlData.master2Ml[masterCols[i]]);
    }
  } else if( row < (3*((mlData.masterVolumeNodes.size()) + (mlData.masterContactNodes.size()))) ) {
    //MasterContact
    AMP::LinearAlgebra::Matrix::shared_ptr masterMat = (mlData.masterOperator)->getMatrix();
    AMP::LinearAlgebra::Matrix::shared_ptr slaveMat = (mlData.slaveOperator)->getMatrix();
    std::vector<unsigned int> masterCols;
    std::vector<unsigned int> slaveCols;
    std::vector<double> masterVals;
    std::vector<double> slaveVals;
    masterMat->getRowByGlobalID(mlData.ml2MasterOrSlave[row], masterCols, masterVals);
    slaveMat->getRowByGlobalID(mlData.masterContact2SlaveContact[mlData.ml2MasterOrSlave[row]], slaveCols, slaveVals);

    for(size_t i = 0; i < (masterCols.size()); i++) {
      unsigned int key = mlData.master2Ml[masterCols[i]];
      std::vector<unsigned int>::iterator pos = std::lower_bound(cols.begin(), cols.end(), key);
      if(pos == cols.end()) {
        cols.push_back(key);
        values.push_back(masterVals[i]);
      } else if((*pos) == key) {
        AMP_ERROR("This should not happen.");
      } else {
        AMP_ASSERT((*pos) > key);
        values.insert(values.begin() + (pos - cols.begin()), masterVals[i]);
        cols.insert(pos, key);
      }
    }

    for(size_t i = 0; i < (slaveCols.size()); i++) {
      unsigned int key = mlData.slave2Ml[slaveCols[i]];
      std::vector<unsigned int>::iterator pos = std::lower_bound(cols.begin(), cols.end(), key);
      if(pos == cols.end()) {
        cols.push_back(key);
        values.push_back(slaveVals[i]);
      } else if((*pos) == key) {
        values[pos - cols.begin()] += slaveVals[i];
      } else {
        AMP_ASSERT((*pos) > key);
        values.insert(values.begin() + (pos - cols.begin()), slaveVals[i]);
        cols.insert(pos, key);
      }
    }

  } else {
    //SlaveVolume
    AMP::LinearAlgebra::Matrix::shared_ptr slaveMat = (mlData.slaveOperator)->getMatrix();
    std::vector<unsigned int> slaveCols;
    slaveMat->getRowByGlobalID(mlData.ml2MasterOrSlave[row], slaveCols, values);
    for(size_t i = 0; i < (slaveCols.size()); i++) {
      cols.push_back(mlData.slave2Ml[slaveCols[i]]);
    }
  }
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
  double * solArr = &((mlData.pcSolVec)[0]);
  double * rhsArr = &((mlData.pcRhsVec)[0]);

  AMP::LinearAlgebra::Vector::shared_ptr masterFvec = fVec->subsetVectorForVariable(mlData.masterVar);
  AMP::LinearAlgebra::Vector::shared_ptr slaveFvec = fVec->subsetVectorForVariable(mlData.slaveVar);

  AMP::LinearAlgebra::Vector::shared_ptr masterUvec = uVec->subsetVectorForVariable(mlData.masterVar);
  AMP::LinearAlgebra::Vector::shared_ptr slaveUvec = uVec->subsetVectorForVariable(mlData.slaveVar);

  AMP::Mesh::DOFMap::shared_ptr master_dof_map = (mlData.masterMeshAdapter)->getDOFMap(mlData.masterVar);
  AMP::Mesh::DOFMap::shared_ptr slave_dof_map = (mlData.slaveMeshAdapter)->getDOFMap(mlData.slaveVar);

  std::vector<unsigned int> dofs(3);
  dofs[0] = 0; dofs[1] = 1; dofs[2] = 2;

  int arrIdx = 0;
  for(size_t i = 0; i < ((mlData.masterVolumeNodes).size()); i++) {
    AMP::Mesh::LibMeshNode masterNd = (mlData.masterMeshAdapter)->getNode( (mlData.masterVolumeNodes)[i] );
    std::vector<unsigned int> masterGlobalIds;
    master_dof_map->getDOFs(masterNd, masterGlobalIds, dofs);
    for(size_t j = 0; j < dofs.size(); j++) {
      rhsArr[arrIdx] = masterFvec->getLocalValueByGlobalID(masterGlobalIds[j]);
      arrIdx++;
    }
  }

  for(size_t i = 0; i < ((mlData.masterContactNodes).size()); i++) {
    AMP::Mesh::LibMeshNode masterNd = (mlData.masterMeshAdapter)->getNode( (mlData.masterContactNodes)[i] );
    std::vector<unsigned int> masterGlobalIds;
    master_dof_map->getDOFs(masterNd, masterGlobalIds, dofs);
    for(size_t j = 0; j < dofs.size(); j++) {
      rhsArr[arrIdx] = masterFvec->getLocalValueByGlobalID(masterGlobalIds[j]);
      arrIdx++;
    }
  }

  for(size_t i = 0; i < ((mlData.slaveVolumeNodes).size()); i++) {
    AMP::Mesh::LibMeshNode slaveNd = (mlData.slaveMeshAdapter)->getNode( (mlData.slaveVolumeNodes)[i] );
    std::vector<unsigned int> slaveGlobalIds;
    slave_dof_map->getDOFs(slaveNd, slaveGlobalIds, dofs);
    for(size_t j = 0; j < dofs.size(); j++) {
      rhsArr[arrIdx] = slaveFvec->getLocalValueByGlobalID(slaveGlobalIds[j]);
      arrIdx++;
    }
  }

  for(size_t i = 0; i < mlData.vecSize; i++) {
    solArr[i] = 0;
  }

  ML_Iterate(ml_object, solArr, rhsArr);

  arrIdx = 0;
  for(size_t i = 0; i < ((mlData.masterVolumeNodes).size()); i++) {
    AMP::Mesh::LibMeshNode masterNd = (mlData.masterMeshAdapter)->getNode( (mlData.masterVolumeNodes)[i] );
    std::vector<unsigned int> masterGlobalIds;
    master_dof_map->getDOFs(masterNd, masterGlobalIds, dofs);
    for(size_t j = 0; j < dofs.size(); j++) {
      masterUvec->setLocalValueByGlobalID(masterGlobalIds[j], solArr[arrIdx]);
      arrIdx++;
    }
  }

  for(size_t i = 0; i < ((mlData.masterContactNodes).size()); i++) {
    AMP::Mesh::LibMeshNode masterNd = (mlData.masterMeshAdapter)->getNode( (mlData.masterContactNodes)[i] );
    std::vector<unsigned int> masterGlobalIds;
    master_dof_map->getDOFs(masterNd, masterGlobalIds, dofs);
    for(size_t j = 0; j < dofs.size(); j++) {
      masterUvec->setLocalValueByGlobalID(masterGlobalIds[j], solArr[arrIdx]);
      arrIdx++;
    }
  }

  for(size_t i = 0; i < ((mlData.slaveVolumeNodes).size()); i++) {
    AMP::Mesh::LibMeshNode slaveNd = (mlData.slaveMeshAdapter)->getNode( (mlData.slaveVolumeNodes)[i] );
    std::vector<unsigned int> slaveGlobalIds;
    slave_dof_map->getDOFs(slaveNd, slaveGlobalIds, dofs);
    for(size_t j = 0; j < dofs.size(); j++) {
      slaveUvec->setLocalValueByGlobalID(slaveGlobalIds[j], solArr[arrIdx]);
      arrIdx++;
    }
  }

  if ( uVec->isA<AMP::LinearAlgebra::DataChangeFirer>() )
  {
    uVec->castTo<AMP::LinearAlgebra::DataChangeFirer>().fireDataChange();
  }

}

void computeMl2MasterOrSlave() {
  AMP::Mesh::DOFMap::shared_ptr master_dof_map = (mlData.masterMeshAdapter)->getDOFMap(mlData.masterVar);
  AMP::Mesh::DOFMap::shared_ptr slave_dof_map = (mlData.slaveMeshAdapter)->getDOFMap(mlData.slaveVar);

  std::vector<unsigned int> dofs(3);
  dofs[0] = 0; dofs[1] = 1; dofs[2] = 2;

  int arrIdx = 0;
  for(size_t i = 0; i < (mlData.masterVolumeNodes.size()); i++) {
    AMP::Mesh::LibMeshNode masterNd = (mlData.masterMeshAdapter)->getNode( (mlData.masterVolumeNodes)[i] );
    std::vector<unsigned int> masterGlobalIds;
    master_dof_map->getDOFs(masterNd, masterGlobalIds, dofs);
    for(size_t j = 0; j < dofs.size(); j++) {
      (mlData.ml2MasterOrSlave)[arrIdx] = masterGlobalIds[j];
      arrIdx++;
    }
  }

  for(size_t i = 0; i < (mlData.masterContactNodes.size()); i++) {
    AMP::Mesh::LibMeshNode masterNd = (mlData.masterMeshAdapter)->getNode( (mlData.masterContactNodes)[i] );
    std::vector<unsigned int> masterGlobalIds;
    master_dof_map->getDOFs(masterNd, masterGlobalIds, dofs);
    for(size_t j = 0; j < dofs.size(); j++) {
      (mlData.ml2MasterOrSlave)[arrIdx] = masterGlobalIds[j];
      arrIdx++;
    }
  }

  for(size_t i = 0; i < (mlData.slaveVolumeNodes.size()); i++) {
    AMP::Mesh::LibMeshNode slaveNd = (mlData.slaveMeshAdapter)->getNode( (mlData.slaveVolumeNodes)[i] );
    std::vector<unsigned int> slaveGlobalIds;
    slave_dof_map->getDOFs(slaveNd, slaveGlobalIds, dofs);
    for(size_t j = 0; j < dofs.size(); j++) {
      (mlData.ml2MasterOrSlave)[arrIdx] = slaveGlobalIds[j];
      arrIdx++;
    }
  }

}

void computeMaster2Ml() {
  AMP::Mesh::DOFMap::shared_ptr master_dof_map = (mlData.masterMeshAdapter)->getDOFMap(mlData.masterVar);

  std::vector<unsigned int> dofs(3);
  dofs[0] = 0; dofs[1] = 1; dofs[2] = 2;

  int arrIdx = 0;
  for(size_t i = 0; i < (mlData.masterVolumeNodes.size()); i++) {
    AMP::Mesh::LibMeshNode masterNd = (mlData.masterMeshAdapter)->getNode( (mlData.masterVolumeNodes)[i] );
    std::vector<unsigned int> masterGlobalIds;
    master_dof_map->getDOFs(masterNd, masterGlobalIds, dofs);
    for(size_t j = 0; j < dofs.size(); j++) {
      (mlData.master2Ml)[masterGlobalIds[j]] = arrIdx;
      arrIdx++;
    }
  }

  for(size_t i = 0; i < (mlData.masterContactNodes.size()); i++) {
    AMP::Mesh::LibMeshNode masterNd = (mlData.masterMeshAdapter)->getNode( (mlData.masterContactNodes)[i] );
    std::vector<unsigned int> masterGlobalIds;
    master_dof_map->getDOFs(masterNd, masterGlobalIds, dofs);
    for(size_t j = 0; j < dofs.size(); j++) {
      (mlData.master2Ml)[masterGlobalIds[j]] = arrIdx;
      arrIdx++;
    }
  }

}

void computeSlave2Ml() {
  AMP::Mesh::DOFMap::shared_ptr master_dof_map = (mlData.masterMeshAdapter)->getDOFMap(mlData.masterVar);
  AMP::Mesh::DOFMap::shared_ptr slave_dof_map = (mlData.slaveMeshAdapter)->getDOFMap(mlData.slaveVar);

  std::vector<unsigned int> dofs(3);
  dofs[0] = 0; dofs[1] = 1; dofs[2] = 2;

  int arrIdx = (3*(mlData.masterVolumeNodes.size() + mlData.masterContactNodes.size()));
  for(size_t i = 0; i < (mlData.slaveVolumeNodes.size()); i++) {
    AMP::Mesh::LibMeshNode slaveNd = (mlData.slaveMeshAdapter)->getNode( (mlData.slaveVolumeNodes)[i] );
    std::vector<unsigned int> slaveGlobalIds;
    slave_dof_map->getDOFs(slaveNd, slaveGlobalIds, dofs);
    for(size_t j = 0; j < dofs.size(); j++) {
      (mlData.slave2Ml)[slaveGlobalIds[j]] = arrIdx;
      arrIdx++;
    }
  }

  for(size_t i = 0; i < (mlData.slaveContactNodes.size()); i++) {
    AMP::Mesh::LibMeshNode slaveNd = (mlData.slaveMeshAdapter)->getNode( (mlData.slaveContactNodes)[i] );
    AMP::Mesh::LibMeshNode masterNd = (mlData.masterMeshAdapter)->getNode( (mlData.masterContactNodes)[i] );
    std::vector<unsigned int> slaveGlobalIds;
    slave_dof_map->getDOFs(slaveNd, slaveGlobalIds, dofs);
    std::vector<unsigned int> masterGlobalIds;
    master_dof_map->getDOFs(masterNd, masterGlobalIds, dofs);
    for(size_t j = 0; j < dofs.size(); j++) {
      (mlData.slave2Ml)[slaveGlobalIds[j]] = (mlData.master2Ml)[masterGlobalIds[j]];
    }
  }

}

void computeMasterContact2SlaveContact() {
  AMP::Mesh::DOFMap::shared_ptr master_dof_map = (mlData.masterMeshAdapter)->getDOFMap(mlData.masterVar);
  AMP::Mesh::DOFMap::shared_ptr slave_dof_map = (mlData.slaveMeshAdapter)->getDOFMap(mlData.slaveVar);

  std::vector<unsigned int> dofs(3);
  dofs[0] = 0; dofs[1] = 1; dofs[2] = 2;

  for(size_t i = 0; i < (mlData.masterContactNodes.size()); i++) {
    AMP::Mesh::LibMeshNode slaveNd = (mlData.slaveMeshAdapter)->getNode( (mlData.slaveContactNodes)[i] );
    AMP::Mesh::LibMeshNode masterNd = (mlData.masterMeshAdapter)->getNode( (mlData.masterContactNodes)[i] );
    std::vector<unsigned int> slaveGlobalIds;
    slave_dof_map->getDOFs(slaveNd, slaveGlobalIds, dofs);
    std::vector<unsigned int> masterGlobalIds;
    master_dof_map->getDOFs(masterNd, masterGlobalIds, dofs);
    for(size_t j = 0; j < dofs.size(); j++) {
      (mlData.masterContact2SlaveContact)[masterGlobalIds[j]] = slaveGlobalIds[j];
    }
  }
}

void pcNoneApply(AMP::LinearAlgebra::Vector::shared_ptr fVec,
    AMP::LinearAlgebra::Vector::shared_ptr uVec) {
  //Identity PC
  uVec->copyVector(fVec);
}

#endif

void myTest(AMP::UnitTest *ut, std::string exeName) {
#if 0
  std::string input_file = "input_" + exeName;
  std::string log_file = "output_" + exeName;

  AMP::PIO::logOnlyNodeZero(log_file);

  boost::shared_ptr<AMP::InputDatabase> input_db(new AMP::InputDatabase("input_db"));
  AMP::InputManager::getManager()->parseInputFile(input_file, input_db);
  input_db->printClassData(AMP::plog);

  boost::shared_ptr<AMP::Database> tiedSurface_db = input_db->getDatabase("TiedSurface_1");
  const std::string masterMeshName = tiedSurface_db->getString("MasterMesh");
  const std::string slaveMeshName = tiedSurface_db->getString("SlaveMesh");
  const unsigned int masterId = tiedSurface_db->getInteger("MasterId");
  const unsigned int slaveId = tiedSurface_db->getInteger("SlaveId");

  AMP::Mesh::MeshManagerParameters::shared_ptr meshmgrParams ( new AMP::Mesh::MeshManagerParameters ( input_db ) );
  AMP::Mesh::MeshManager::shared_ptr manager ( new AMP::Mesh::MeshManager ( meshmgrParams ) );
  AMP::Mesh::MeshManager::Adapter::shared_ptr masterMeshAdapter = manager->getMesh ( masterMeshName );
  AMP::Mesh::MeshManager::Adapter::shared_ptr slaveMeshAdapter = manager->getMesh ( slaveMeshName );

  std::vector<unsigned int> masterContactNodes;
  std::vector<unsigned int> slaveContactNodes;
  std::vector<unsigned int> masterVolumeNodes;
  std::vector<unsigned int> slaveVolumeNodes;

  createMasterSlaveMap(masterMeshAdapter, slaveMeshAdapter, masterId, slaveId,
      masterContactNodes, slaveContactNodes, masterVolumeNodes, slaveVolumeNodes);

  boost::shared_ptr<AMP::Operator::ElementPhysicsModel> masterElementPhysicsModel;
  boost::shared_ptr<AMP::Operator::LinearBVPOperator> masterOperator = boost::dynamic_pointer_cast<
    AMP::Operator::LinearBVPOperator>(AMP::Operator::OperatorBuilder::createOperator(
          masterMeshAdapter, "MasterBVPOperator", input_db, masterElementPhysicsModel));

  boost::shared_ptr<AMP::Operator::ElementPhysicsModel> slaveElementPhysicsModel;
  boost::shared_ptr<AMP::Operator::LinearBVPOperator> slaveOperator = boost::dynamic_pointer_cast<
    AMP::Operator::LinearBVPOperator>(AMP::Operator::OperatorBuilder::createOperator(
          slaveMeshAdapter, "SlaveBVPOperator", input_db, slaveElementPhysicsModel));

  boost::shared_ptr<AMP::Operator::OperatorParameters> dummyParams;
  boost::shared_ptr<AMP::Operator::ColumnOperator> columnOperator(new AMP::Operator::ColumnOperator(dummyParams));
  columnOperator->append(masterOperator);
  columnOperator->append(slaveOperator);

  AMP::LinearAlgebra::Variable::shared_ptr masterVar = masterOperator->getOutputVariable();
  AMP::LinearAlgebra::Variable::shared_ptr slaveVar = slaveOperator->getOutputVariable();
  AMP::LinearAlgebra::Variable::shared_ptr columnVar = columnOperator->getOutputVariable();

  boost::shared_ptr<AMP::Operator::ElementPhysicsModel> dummyModel;
  boost::shared_ptr<AMP::Operator::DirichletVectorCorrection> loadOperator = boost::dynamic_pointer_cast<
    AMP::Operator::DirichletVectorCorrection>(AMP::Operator::OperatorBuilder::createOperator(
          slaveMeshAdapter, "LoadOperator", input_db, dummyModel));
  loadOperator->setVariable(slaveVar);

  AMP::LinearAlgebra::Vector::shared_ptr nullVec;
  AMP::LinearAlgebra::Vector::shared_ptr columnSolVec = manager->createVector(columnVar);
  AMP::LinearAlgebra::Vector::shared_ptr columnRhsVec = columnSolVec->cloneVector();
  AMP::LinearAlgebra::Vector::shared_ptr columnResVec = columnSolVec->cloneVector();

  columnSolVec->zero();
  columnRhsVec->zero();
  columnResVec->zero();

  loadOperator->apply(nullVec, nullVec, columnRhsVec, 1.0, 0.0);

  AMP::LinearAlgebra::Vector::shared_ptr masterRhsVec = columnRhsVec->subsetVectorForVariable(masterVar);
  AMP::LinearAlgebra::Vector::shared_ptr slaveRhsVec = columnRhsVec->subsetVectorForVariable(slaveVar);

  std::cout<<"before correction: master RHS norm = "<<(masterRhsVec->L2Norm())<<std::endl;
  std::cout<<"before correction: slave RHS norm = "<<(slaveRhsVec->L2Norm())<<std::endl;
  std::cout<<"before correction: column RHS norm = "<<(columnRhsVec->L2Norm())<<std::endl;

  int maxCGiters = input_db->getInteger("maxCGiters");

  bool useMLpc = input_db->getBool("USE_ML_PC");

  if(useMLpc) {
    int maxMLiters = 1;
    int printInfoLevel = 10;
    int numGrids = 10;
    int numPDEs = 3;
    int maxCoarseSize = 128;
    size_t vecSize = (3*(masterVolumeNodes.size() + 
          slaveVolumeNodes.size() + masterContactNodes.size()));

    mlData.columnOperator = columnOperator;
    mlData.masterOperator = masterOperator;
    mlData.slaveOperator = slaveOperator;
    mlData.masterContactNodes = masterContactNodes;
    mlData.slaveContactNodes = slaveContactNodes;
    mlData.masterVolumeNodes = masterVolumeNodes;
    mlData.slaveVolumeNodes = slaveVolumeNodes;
    mlData.masterVar = masterVar;
    mlData.slaveVar = slaveVar;
    mlData.columnVar = columnVar;
    mlData.masterMeshAdapter = masterMeshAdapter;
    mlData.slaveMeshAdapter = slaveMeshAdapter;
    mlData.manager = manager;
    mlData.matVecInVec = manager->createVector(columnVar);
    mlData.matVecOutVec = manager->createVector(columnVar);
    (mlData.getColInVec).resize(vecSize);
    (mlData.getColOutVec).resize(vecSize);
    (mlData.pcSolVec).resize(vecSize);
    (mlData.pcRhsVec).resize(vecSize);
    (mlData.ml2MasterOrSlave).resize(vecSize);
    (mlData.master2Ml).resize(masterRhsVec->getLocalSize());
    (mlData.slave2Ml).resize(slaveRhsVec->getLocalSize());
    (mlData.masterContact2SlaveContact).resize(masterRhsVec->getLocalSize());
    mlData.vecSize = vecSize;

    computeMl2MasterOrSlave();
    computeMaster2Ml();
    computeSlave2Ml();
    computeMasterContact2SlaveContact();

    ML_Create (&ml_object, numGrids);
    ML_Init_Amatrix(ml_object, 0, vecSize, vecSize, NULL);
    ML_Set_Amatrix_Getrow(ml_object, 0, &myGetMultipleRows, NULL, vecSize);
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

  {
    std::cout<<"MPC-MLCG: "<<std::endl;
    AMP::LinearAlgebra::Vector::shared_ptr MatOutVec = columnSolVec->cloneVector();
    AMP::LinearAlgebra::Vector::shared_ptr pVec = columnSolVec->cloneVector();
    AMP::LinearAlgebra::Vector::shared_ptr zVec = columnSolVec->cloneVector();

    addSlaveToMaster(slaveVar, masterVar, columnRhsVec, slaveMeshAdapter, 
        masterMeshAdapter, slaveContactNodes, masterContactNodes);
    setSlaveToZero(slaveVar, columnRhsVec, slaveMeshAdapter, slaveContactNodes);

    std::cout<<"after correction: master RHS norm = "<<(masterRhsVec->L2Norm())<<std::endl;
    std::cout<<"after correction: slave RHS norm = "<<(slaveRhsVec->L2Norm())<<std::endl;
    std::cout<<"after correction: column RHS norm = "<<(columnRhsVec->L2Norm())<<std::endl;

    copyMasterToSlave(slaveVar, masterVar, columnSolVec, slaveMeshAdapter,
        masterMeshAdapter, slaveContactNodes, masterContactNodes);

    columnOperator->apply(nullVec, columnSolVec, MatOutVec, 1.0, 0.0);
    addSlaveToMaster(slaveVar, masterVar, MatOutVec, slaveMeshAdapter, 
        masterMeshAdapter, slaveContactNodes, masterContactNodes);
    setSlaveToZero(slaveVar, MatOutVec, slaveMeshAdapter, slaveContactNodes);

    columnResVec->subtract(columnRhsVec, MatOutVec);

    if(useMLpc) {
      pcMLapply(columnResVec, zVec);
    } else {
      pcNoneApply(columnResVec, zVec);
    }
    copyMasterToSlave(slaveVar, masterVar, zVec, slaveMeshAdapter, 
        masterMeshAdapter, slaveContactNodes, masterContactNodes);

    pVec->copyVector(zVec);

    for(int iter = 0; iter <= maxCGiters; iter++) {
      double resNorm = columnResVec->L2Norm();

      std::cout<<"CG-Iter = "<<iter<<" ResNorm2 = "<<std::setprecision(15)<<resNorm<<std::endl;

      columnOperator->apply(nullVec, pVec, MatOutVec, 1.0, 0.0);
      addSlaveToMaster(slaveVar, masterVar, MatOutVec, slaveMeshAdapter, 
          masterMeshAdapter, slaveContactNodes, masterContactNodes);
      setSlaveToZero(slaveVar, MatOutVec, slaveMeshAdapter, slaveContactNodes);

      double matOutNorm = MatOutVec->L2Norm();
      std::cout<<"CG-Iter = "<<iter<<" MatOutNorm2 = "<<std::setprecision(15)<<matOutNorm<<std::endl;

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
      copyMasterToSlave(slaveVar, masterVar, zVec, slaveMeshAdapter, 
          masterMeshAdapter, slaveContactNodes, masterContactNodes);

      double resNewDotZ = columnResVec->dot(zVec);

      double beta = resNewDotZ/resOldDotZ;

      std::cout<<"CG-Iter = "<<iter
        <<" resOldDotZ = "<<std::setprecision(15)<<resOldDotZ
        <<" alphaDenom = "<<std::setprecision(15)<<alphaDenom
        <<" alpha = "<<std::setprecision(15)<<alpha
        <<" resNewDotZ = "<<std::setprecision(15)<<resNewDotZ
        <<" beta = "<<std::setprecision(15)<<beta
        <<std::endl<<std::endl;

      pVec->axpy(beta, pVec, zVec);
    }
  }

  if(useMLpc)
  {
    ML_Aggregate_Destroy(&agg_object);
    ML_Destroy(&ml_object);
  }

  freeMLdata();

#ifdef USE_SILO
  manager->registerVectorAsData ( columnSolVec, "Displacement" );
  manager->writeFile<AMP::Mesh::SiloIO> ( exeName , 1 );
#endif

#endif
  ut->passes(exeName);
}

int main(int argc, char *argv[])
{
  AMP::AMPManager::startup(argc, argv);
  AMP::UnitTest ut;

  AMP_ERROR("Not yet converted!");

  std::string exeName = "testTiedContactV2";

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


