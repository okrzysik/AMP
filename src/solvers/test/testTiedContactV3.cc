
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

#include "boost/shared_ptr.hpp"

#include "operators/ColumnOperator.h"
#include "operators/LinearBVPOperator.h"
#include "operators/OperatorBuilder.h"
#include "operators/boundary/DirichletVectorCorrection.h"

#include "ml_include.h"
#include "solvers/MLoptions.h"

#include "TiedContactUtils.h"
#include "MPCtestUtils.h"
#include "indexHolder.h"

ML_Aggregate *agg_object;
ML* ml_object;

struct MyData {
  boost::shared_ptr<AMP::Operator::ColumnOperator> columnOperator; 
  std::vector<AMP::Mesh::MeshManager::Adapter::shared_ptr> meshAdapters;
  AMP::Mesh::MeshManager::shared_ptr manager;
  std::vector<unsigned int> masterContactNodeIds;
  std::vector<unsigned int> masterContactMeshIds;
  std::vector<unsigned int> slaveContactNodeIds;
  std::vector<unsigned int> slaveContactMeshIds;
  std::vector<unsigned int> slave2MasterMap;
  std::vector<std::vector<unsigned int> > nonContactNodeList;
  size_t vecSize;
} myData;

void freeMyData() {
  myData.columnOperator.reset();
  myData.manager.reset();
  myData.meshAdapters.clear();
}

void addSlaveToMaster(AMP::LinearAlgebra::Vector::shared_ptr vec) {
  std::vector<AMP::LinearAlgebra::Variable::shared_ptr> vars;
  vars.resize(myData.meshAdapters.size());
  for(size_t i = 0; i < vars.size(); i++) {
    vars[i] = ((myData.columnOperator)->getOperator(i))->getOutputVariable();
  }

  std::vector<AMP::Mesh::DOFMap::shared_ptr> dof_maps;
  dof_maps.resize(vars.size());
  for(size_t i = 0; i < vars.size(); i++) {
    dof_maps[i] = (myData.meshAdapters[i])->getDOFMap(vars[i]);
  }

  std::vector<AMP::LinearAlgebra::Vector::shared_ptr> subVecs;
  subVecs.resize(vars.size());
  for(size_t i = 0; i < vars.size(); i++) {
    subVecs[i] = vec->subsetVectorForVariable(vars[i]);
  }

  std::vector<unsigned int> dofs(3);
  dofs[0] = 0; dofs[1] = 1; dofs[2] = 2;

  for(size_t i = 0; i < myData.slaveContactNodeIds.size(); i++) {
    unsigned int slaveMeshId = myData.slaveContactMeshIds[i];
    unsigned int slaveNodeId = myData.slaveContactNodeIds[i];
    unsigned int masterIndex = myData.slave2MasterMap[i];
    unsigned int masterMeshId = myData.masterContactMeshIds[masterIndex];
    unsigned int masterNodeId = myData.masterContactNodeIds[masterIndex];
    AMP::Mesh::MeshManager::Adapter::shared_ptr slaveMeshAdapter = myData.meshAdapters[slaveMeshId];
    AMP::Mesh::MeshManager::Adapter::shared_ptr masterMeshAdapter = myData.meshAdapters[masterMeshId];
    AMP::LinearAlgebra::Vector::shared_ptr slaveVec = subVecs[slaveMeshId];
    AMP::LinearAlgebra::Vector::shared_ptr masterVec = subVecs[masterMeshId];
    AMP::Mesh::DOFMap::shared_ptr slave_dof_map = dof_maps[slaveMeshId];
    AMP::Mesh::DOFMap::shared_ptr master_dof_map = dof_maps[masterMeshId];
    AMP::Mesh::LibMeshNode slaveNd = slaveMeshAdapter->getNode(slaveNodeId);
    AMP::Mesh::LibMeshNode masterNd = masterMeshAdapter->getNode(masterNodeId);
    std::vector<unsigned int> slaveGlobalIds;
    std::vector<unsigned int> masterGlobalIds;
    slave_dof_map->getDOFs(slaveNd, slaveGlobalIds, dofs);
    master_dof_map->getDOFs(masterNd, masterGlobalIds, dofs);
    for(size_t j = 0; j < dofs.size(); j++) {
      double slaveVal = slaveVec->getLocalValueByGlobalID(slaveGlobalIds[j]);
      masterVec->addLocalValueByGlobalID(masterGlobalIds[j], slaveVal);
    }//end for j
  }//end for i
  vec->makeConsistent( AMP::LinearAlgebra::Vector::CONSISTENT_ADD );
}

void copyMasterToSlave(AMP::LinearAlgebra::Vector::shared_ptr vec) {
  std::vector<AMP::LinearAlgebra::Variable::shared_ptr> vars;
  vars.resize(myData.meshAdapters.size());
  for(size_t i = 0; i < vars.size(); i++) {
    vars[i] = ((myData.columnOperator)->getOperator(i))->getOutputVariable();
  }

  std::vector<AMP::Mesh::DOFMap::shared_ptr> dof_maps;
  dof_maps.resize(vars.size());
  for(size_t i = 0; i < vars.size(); i++) {
    dof_maps[i] = (myData.meshAdapters[i])->getDOFMap(vars[i]);
  }

  std::vector<AMP::LinearAlgebra::Vector::shared_ptr> subVecs;
  subVecs.resize(vars.size());
  for(size_t i = 0; i < vars.size(); i++) {
    subVecs[i] = vec->subsetVectorForVariable(vars[i]);
  }

  std::vector<unsigned int> dofs(3);
  dofs[0] = 0; dofs[1] = 1; dofs[2] = 2;

  for(size_t i = 0; i < myData.slaveContactNodeIds.size(); i++) {
    unsigned int slaveMeshId = myData.slaveContactMeshIds[i];
    unsigned int slaveNodeId = myData.slaveContactNodeIds[i];
    unsigned int masterIndex = myData.slave2MasterMap[i];
    unsigned int masterMeshId = myData.masterContactMeshIds[masterIndex];
    unsigned int masterNodeId = myData.masterContactNodeIds[masterIndex];
    AMP::Mesh::MeshManager::Adapter::shared_ptr slaveMeshAdapter = myData.meshAdapters[slaveMeshId];
    AMP::Mesh::MeshManager::Adapter::shared_ptr masterMeshAdapter = myData.meshAdapters[masterMeshId];
    AMP::Mesh::DOFMap::shared_ptr slave_dof_map = dof_maps[slaveMeshId];
    AMP::Mesh::DOFMap::shared_ptr master_dof_map = dof_maps[masterMeshId];
    AMP::LinearAlgebra::Vector::shared_ptr slaveVec = subVecs[slaveMeshId];
    AMP::LinearAlgebra::Vector::shared_ptr masterVec = subVecs[masterMeshId];
    AMP::Mesh::LibMeshNode slaveNd = slaveMeshAdapter->getNode(slaveNodeId);
    AMP::Mesh::LibMeshNode masterNd = masterMeshAdapter->getNode(masterNodeId);
    std::vector<unsigned int> slaveGlobalIds;
    std::vector<unsigned int> masterGlobalIds;
    slave_dof_map->getDOFs(slaveNd, slaveGlobalIds, dofs);
    master_dof_map->getDOFs(masterNd, masterGlobalIds, dofs);
    for(size_t j = 0; j < dofs.size(); j++) {
      double masterVal = masterVec->getLocalValueByGlobalID(masterGlobalIds[j]);
      slaveVec->setLocalValueByGlobalID(slaveGlobalIds[j], masterVal);
    }//end for j
  }//end for i
  vec->makeConsistent( AMP::LinearAlgebra::Vector::CONSISTENT_SET );
}

void setSlaveToZero(AMP::LinearAlgebra::Vector::shared_ptr vec) {
  std::vector<AMP::LinearAlgebra::Variable::shared_ptr> vars;
  vars.resize(myData.meshAdapters.size());
  for(size_t i = 0; i < vars.size(); i++) {
    vars[i] = ((myData.columnOperator)->getOperator(i))->getOutputVariable();
  }

  std::vector<AMP::Mesh::DOFMap::shared_ptr> dof_maps;
  dof_maps.resize(vars.size());
  for(size_t i = 0; i < vars.size(); i++) {
    dof_maps[i] = (myData.meshAdapters[i])->getDOFMap(vars[i]);
  }

  std::vector<AMP::LinearAlgebra::Vector::shared_ptr> subVecs;
  subVecs.resize(vars.size());
  for(size_t i = 0; i < vars.size(); i++) {
    subVecs[i] = vec->subsetVectorForVariable(vars[i]);
  }

  std::vector<unsigned int> dofs(3);
  dofs[0] = 0; dofs[1] = 1; dofs[2] = 2;

  for(size_t i = 0; i < myData.slaveContactNodeIds.size(); i++) {
    unsigned int slaveMeshId = myData.slaveContactMeshIds[i];
    unsigned int slaveNodeId = myData.slaveContactNodeIds[i];
    AMP::Mesh::MeshManager::Adapter::shared_ptr slaveMeshAdapter = myData.meshAdapters[slaveMeshId];
    AMP::Mesh::DOFMap::shared_ptr slave_dof_map = dof_maps[slaveMeshId];
    AMP::LinearAlgebra::Vector::shared_ptr slaveVec = subVecs[slaveMeshId];
    AMP::Mesh::LibMeshNode slaveNd =  slaveMeshAdapter->getNode(slaveNodeId);
    std::vector<unsigned int> slaveGlobalIds;
    slave_dof_map->getDOFs(slaveNd, slaveGlobalIds, dofs);
    for(size_t j = 0; j < dofs.size(); j++) {
      slaveVec->setLocalValueByGlobalID(slaveGlobalIds[j], 0.0);
    }//end for j
  }//end for i
  vec->makeConsistent( AMP::LinearAlgebra::Vector::CONSISTENT_SET );
}

void pcNoneApply(AMP::LinearAlgebra::Vector::shared_ptr fVec,
    AMP::LinearAlgebra::Vector::shared_ptr uVec) {
  //Identity PC
  uVec->copyVector(fVec);
}

int myMatVec(ML_Operator *data, int in_length, double in[], int out_length, double out[]) {
  //ToDo

  return 0;
}

void myGetSingleRow(int row, std::vector<unsigned int> &cols, std::vector<double> &values) {
  cols.clear();
  values.clear();
  //ToDo

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
}

void createMasterAndSlaveContactNodes(std::vector<AMP::Mesh::MeshManager::Adapter::shared_ptr> meshAdapters,
    std::vector<unsigned int> tiedMeshIds, std::vector<unsigned int> tiedBndIds,  
    std::vector<unsigned int> & masterContactNodeIds, std::vector<unsigned int> & masterContactMeshIds, 
    std::vector<unsigned int> & slaveContactNodeIds, std::vector<unsigned int> & slaveContactMeshIds,
    std::vector<unsigned int> & slave2MasterMap) {

  std::vector<PointAndId> tiedPointSet;

  std::vector<unsigned int> contactNodeIds;
  std::vector<unsigned int> contactMeshIds;

  for(size_t i = 0; i < tiedBndIds.size(); i++) {
    AMP::Mesh::MeshManager::Adapter::OwnedBoundaryNodeIterator bnd = 
      meshAdapters[tiedMeshIds[i]]->beginOwnedBoundary( tiedBndIds[i] );
    AMP::Mesh::MeshManager::Adapter::OwnedBoundaryNodeIterator end_bnd = 
      meshAdapters[tiedMeshIds[i]]->endOwnedBoundary( tiedBndIds[i] );
    for( ; bnd != end_bnd; ++bnd) {
      PointAndId tmp;
      tmp.d_id = contactNodeIds.size(); 
      contactNodeIds.push_back(bnd->globalID());
      contactMeshIds.push_back(tiedMeshIds[i]);
      tmp.d_pt[0] = bnd->x();
      tmp.d_pt[1] = bnd->y();
      tmp.d_pt[2] = bnd->z();
      tiedPointSet.push_back(tmp);
    }//end for bnd
  }

  std::sort(tiedPointSet.begin(), tiedPointSet.end());

  std::cout<<"Number of contact nodes (before removing duplicates) = " <<(tiedPointSet.size())<<std::endl;

  std::vector<PointAndId> tmpPointSet;
  tmpPointSet.push_back(tiedPointSet[0]);
  for(size_t i = 1; i < tiedPointSet.size(); i++) {
    unsigned int tmpMeshId = contactMeshIds[tmpPointSet[tmpPointSet.size() - 1].d_id];
    std::vector<IndexHolder<unsigned int> > tmpVec;   
    while((i < tiedPointSet.size()) &&
        (tiedPointSet[i].pointEquals(tiedPointSet[i - 1]))) {
      unsigned int currMeshId = contactMeshIds[tiedPointSet[i].d_id];
      if(tmpMeshId != currMeshId) {
        IndexHolder<unsigned int> tmpElem;
        tmpElem.value = &(contactMeshIds[tiedPointSet[i].d_id]);
        tmpElem.index = i;
        tmpVec.push_back(tmpElem);
      }
      i++;
    }
    std::sort(tmpVec.begin(), tmpVec.end());
    std::vector<unsigned int> tmpIndices;
    tmpIndices.push_back(tmpVec[0].index);
    for(size_t j = 1; j < tmpVec.size(); j++) {
      if(tmpVec[j] != tmpVec[j - 1]) {
        tmpIndices.push_back(tmpVec[j].index);
      }
    }
    std::sort(tmpIndices.begin(), tmpIndices.end());
    for(size_t j = 0; j < tmpIndices.size(); j++) {
      tmpPointSet.push_back(tiedPointSet[tmpIndices[j]]);
    }
    if(i < tiedPointSet.size()) {
      tmpPointSet.push_back(tiedPointSet[i]);
    }
  }
  tiedPointSet = tmpPointSet;
  tmpPointSet.clear();

  masterContactNodeIds.clear();
  masterContactMeshIds.clear();
  slaveContactNodeIds.clear();
  slaveContactMeshIds.clear();
  slave2MasterMap.clear();

  masterContactNodeIds.push_back(contactNodeIds[tiedPointSet[0].d_id]);
  masterContactMeshIds.push_back(contactMeshIds[tiedPointSet[0].d_id]);
  for(size_t i = 1; i < tiedPointSet.size(); i++) {
    if(tiedPointSet[i].pointEquals(tiedPointSet[i - 1])) {
      slaveContactNodeIds.push_back(contactNodeIds[tiedPointSet[i].d_id]);
      slaveContactMeshIds.push_back(contactMeshIds[tiedPointSet[i].d_id]);
      slave2MasterMap.push_back(masterContactNodeIds.size() - 1);
    } else {
      masterContactNodeIds.push_back(contactNodeIds[tiedPointSet[i].d_id]);
      masterContactMeshIds.push_back(contactMeshIds[tiedPointSet[i].d_id]);
    }
  }

  std::cout<<"Number of contact nodes (final) = " <<(tiedPointSet.size())<<std::endl;
  std::cout<<"Number of Master nodes = " <<(masterContactNodeIds.size())<<std::endl;
  std::cout<<"Number of Slave nodes = " <<(slaveContactNodeIds.size())<<std::endl;

}

void createNonContactNodeList(std::vector<AMP::Mesh::MeshManager::Adapter::shared_ptr> meshAdapters, 
    std::vector<unsigned int> masterContactNodeIds, std::vector<unsigned int> masterContactMeshIds, 
    std::vector<unsigned int> slaveContactNodeIds, std::vector<unsigned int> slaveContactMeshIds, 
    std::vector<std::vector<unsigned int> > & nonContactNodeList) {
  nonContactNodeList.resize(meshAdapters.size());
  for(size_t i = 0; i < nonContactNodeList.size(); i++) {
    nonContactNodeList[i].resize((meshAdapters[i])->numLocalNodes());
    for(size_t j = 0; j < nonContactNodeList[i].size(); j++) {
      nonContactNodeList[i][j] = j;
    }
  }
  for(size_t i = 0; i < masterContactNodeIds.size(); i++) {
    nonContactNodeList[masterContactMeshIds[i]][masterContactNodeIds[i]] = 
      static_cast<unsigned int>(-1);
  }
  for(size_t i = 0; i < slaveContactNodeIds.size(); i++) {
    nonContactNodeList[slaveContactMeshIds[i]][slaveContactNodeIds[i]] = 
      static_cast<unsigned int>(-1);
  }
  for(size_t i = 0; i < nonContactNodeList.size(); i++) {
    std::vector<unsigned int> tmpVec;
    for(size_t j = 0; j < nonContactNodeList[i].size(); j++) {
      if(nonContactNodeList[i][j] != static_cast<unsigned int>(-1)) {
        tmpVec.push_back(nonContactNodeList[i][j]);
      }
    }
    nonContactNodeList[i] = tmpVec;
  }
}


void myTest(AMP::UnitTest *ut, std::string exeName) {
  std::string input_file = "input_" + exeName;
  std::string log_file = "output_" + exeName;

  AMP::PIO::logOnlyNodeZero(log_file);

  boost::shared_ptr<AMP::InputDatabase> input_db(new AMP::InputDatabase("input_db"));
  AMP::InputManager::getManager()->parseInputFile(input_file, input_db);
  input_db->printClassData(AMP::plog);

  AMP::Mesh::MeshManagerParameters::shared_ptr meshmgrParams ( new AMP::Mesh::MeshManagerParameters ( input_db ) );
  AMP::Mesh::MeshManager::shared_ptr manager ( new AMP::Mesh::MeshManager ( meshmgrParams ) );

  unsigned int numMeshes = input_db->getInteger("NumberOfMeshes");

  std::vector<std::string> meshNames;
  for(unsigned int i = 1; i <= numMeshes; i++) {
    char tmp[200];
    sprintf(tmp, "Mesh_%u", i);
    meshNames.push_back((input_db->getDatabase(tmp))->getString("MeshName"));
  }

  std::vector<AMP::Mesh::MeshManager::Adapter::shared_ptr> meshAdapters;
  for(size_t i = 0; i < meshNames.size(); i++) {
    meshAdapters.push_back(manager->getMesh( meshNames[i] ));
  }

  unsigned int numContactSurfaces = input_db->getInteger("NumberOfContactSurfaces");

  std::vector<unsigned int> tiedMeshIds;
  std::vector<unsigned int> tiedBndIds;

  for(unsigned int i = 1; i <= numContactSurfaces; i++) {
    char tmp1[200];
    sprintf(tmp1, "Contact_Mesh_%u", i);

    char tmp2[200];
    sprintf(tmp2, "Contact_BoundaryId_%u", i);

    tiedMeshIds.push_back( (input_db->getInteger(tmp1)) - 1 );
    tiedBndIds.push_back( input_db->getInteger(tmp2) );
  }

  std::vector<unsigned int> masterContactNodeIds;
  std::vector<unsigned int> masterContactMeshIds;
  std::vector<unsigned int> slaveContactNodeIds;
  std::vector<unsigned int> slaveContactMeshIds;
  std::vector<unsigned int> slave2MasterMap;
  std::vector<std::vector<unsigned int> > nonContactNodeList;

  createMasterAndSlaveContactNodes(meshAdapters, tiedMeshIds, tiedBndIds, masterContactNodeIds,
      masterContactMeshIds, slaveContactNodeIds, slaveContactMeshIds, slave2MasterMap);

  createNonContactNodeList(meshAdapters, masterContactNodeIds, masterContactMeshIds, 
      slaveContactNodeIds, slaveContactMeshIds, nonContactNodeList);

  boost::shared_ptr<AMP::Operator::OperatorParameters> dummyParams;
  boost::shared_ptr<AMP::Operator::ColumnOperator> columnOperator(new 
      AMP::Operator::ColumnOperator(dummyParams));
  for(unsigned int i = 1; i <= numMeshes; i++) {
    char key[200];
    sprintf(key, "BVPoperator_%u", i);
    boost::shared_ptr<AMP::Operator::ElementPhysicsModel> currModel;
    boost::shared_ptr<AMP::Operator::LinearBVPOperator> currOperator = boost::dynamic_pointer_cast<
      AMP::Operator::LinearBVPOperator>(AMP::Operator::OperatorBuilder::createOperator(
            meshAdapters[i - 1], key, input_db, currModel));
    columnOperator->append(currOperator);
  }

  unsigned int loadMeshId = (input_db->getInteger("LoadMesh")) - 1;

  boost::shared_ptr<AMP::Operator::ElementPhysicsModel> dummyModel;
  boost::shared_ptr<AMP::Operator::DirichletVectorCorrection> loadOperator = boost::dynamic_pointer_cast<
    AMP::Operator::DirichletVectorCorrection>(AMP::Operator::OperatorBuilder::createOperator(
          meshAdapters[loadMeshId], "LoadOperator", input_db, dummyModel));
  loadOperator->setVariable((columnOperator->getOperator(loadMeshId))->getOutputVariable());


  AMP::LinearAlgebra::Variable::shared_ptr columnVar = columnOperator->getOutputVariable();

  AMP::LinearAlgebra::Vector::shared_ptr nullVec;
  AMP::LinearAlgebra::Vector::shared_ptr columnSolVec = manager->createVector(columnVar);
  AMP::LinearAlgebra::Vector::shared_ptr columnRhsVec = columnSolVec->cloneVector();
  AMP::LinearAlgebra::Vector::shared_ptr columnResVec = columnSolVec->cloneVector();

  columnSolVec->zero();
  columnRhsVec->zero();
  columnResVec->zero();

  loadOperator->apply(nullVec, nullVec, columnRhsVec, 1.0, 0.0);

  int maxCGiters = input_db->getInteger("maxCGiters");

  bool useMLpc = input_db->getBool("USE_ML_PC");

  size_t vecSize = 0;
  for(size_t i = 0; i < nonContactNodeList.size(); i++) {
    vecSize += nonContactNodeList[i].size();
  }
  vecSize += masterContactNodeIds.size();
  vecSize *= 3;

  myData.columnOperator = columnOperator;
  myData.manager = manager;
  myData.meshAdapters = meshAdapters;
  myData.masterContactNodeIds = masterContactNodeIds;
  myData.masterContactMeshIds = masterContactMeshIds;
  myData.slaveContactNodeIds = slaveContactNodeIds;
  myData.slaveContactMeshIds = slaveContactMeshIds;
  myData.slave2MasterMap = slave2MasterMap;
  myData.nonContactNodeList = nonContactNodeList;
  myData.vecSize = vecSize;

  if(useMLpc) {
    int maxMLiters = 1;
    int printInfoLevel = 10;
    int numGrids = 10;
    int numPDEs = 3;
    int maxCoarseSize = 128;

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

    columnOperator->apply(nullVec, columnSolVec, MatOutVec, 1.0, 0.0);

    columnResVec->subtract(columnRhsVec, MatOutVec);

    addSlaveToMaster(columnResVec);
    setSlaveToZero(columnResVec);

    if(useMLpc) {
      pcMLapply(columnResVec, zVec);
    } else {
      pcNoneApply(columnResVec, zVec);
    }

    copyMasterToSlave(zVec);

    pVec->copyVector(zVec);

    for(int iter = 0; iter <= maxCGiters; iter++) {
      double resNorm = columnResVec->L2Norm();

      std::cout<<"CG-Iter = "<<iter<<" ResNorm2 = "<<std::setprecision(15)<<resNorm<<std::endl;

      columnOperator->apply(nullVec, pVec, MatOutVec, 1.0, 0.0);
      addSlaveToMaster(MatOutVec);
      setSlaveToZero(MatOutVec);

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
      copyMasterToSlave(zVec);

      double resNewDotZ = columnResVec->dot(zVec);

      double beta = resNewDotZ/resOldDotZ;

      pVec->axpy(beta, pVec, zVec);
    }
  }

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

  std::string exeName = "testTiedContactV3";

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



