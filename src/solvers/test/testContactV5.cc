
#include "utils/InputManager.h"
#include "utils/AMPManager.h"
#include "utils/UnitTest.h"
#include "utils/Utilities.h"

#include <iostream>
#include <string>
#include <vector>
#include <algorithm>
#include <cmath>

#include "ampmesh/SiloIO.h"

#include "operators/ColumnOperator.h"
#include "operators/OperatorBuilder.h"
#include "operators/mechanics/MechanicsLinearFEOperator.h"
#include "operators/boundary/DirichletVectorCorrection.h"

#include "fe_interface.h"
#include "cell_hex8.h"

//#include "ContactSearchUtils.h"

#if 0

void pcNoneApply(AMP::LinearAlgebra::Vector::shared_ptr fVec,
    AMP::LinearAlgebra::Vector::shared_ptr uVec) {
  //Identity PC
  uVec->copyVector(fVec);
}

void setDirichletVals(AMP::Mesh::MeshManager::Adapter::shared_ptr masterMeshAdapter, 
    std::vector<unsigned int> const & bndIds,
    std::vector<std::vector<unsigned int> > const & bndDofIds,
    std::vector<std::vector<double> > const & bndDofVals,
    AMP::LinearAlgebra::Variable::shared_ptr masterVar,
    AMP::LinearAlgebra::Vector::shared_ptr vec) {
  AMP::Mesh::DOFMap::shared_ptr master_dof_map = masterMeshAdapter->getDOFMap(masterVar);

  AMP::LinearAlgebra::Vector::shared_ptr masterVec =  vec->subsetVectorForVariable(masterVar);

  for(size_t i = 0; i < bndIds.size(); i++) {
    AMP::Mesh::MeshManager::Adapter::OwnedBoundaryNodeIterator bnd = masterMeshAdapter->beginOwnedBoundary( bndIds[i] );
    AMP::Mesh::MeshManager::Adapter::OwnedBoundaryNodeIterator end_bnd = masterMeshAdapter->endOwnedBoundary( bndIds[i] );

    for( ; bnd != end_bnd; ++bnd) {
      std::vector<unsigned int> bndGlobalIds;
      master_dof_map->getDOFs(*bnd, bndGlobalIds, bndDofIds[i]);

      for(size_t j = 0; j < bndGlobalIds.size(); j++) {
        masterVec->setLocalValueByGlobalID(bndGlobalIds[j], bndDofVals[i][j]);
      }//end for j
    }//end for bnd
  }//end for i

  vec->makeConsistent( AMP::LinearAlgebra::Vector::CONSISTENT_SET );
}

void setDirichletDofsToZero(AMP::Mesh::MeshManager::Adapter::shared_ptr masterMeshAdapter, 
    std::vector<unsigned int> const & bndIds,
    std::vector<std::vector<unsigned int> > const & bndDofIds,
    AMP::LinearAlgebra::Variable::shared_ptr masterVar,
    AMP::LinearAlgebra::Vector::shared_ptr vec) {
  AMP::Mesh::DOFMap::shared_ptr master_dof_map = masterMeshAdapter->getDOFMap(masterVar);

  AMP::LinearAlgebra::Vector::shared_ptr masterVec =  vec->subsetVectorForVariable(masterVar);

  for(size_t i = 0; i < bndIds.size(); i++) {
    AMP::Mesh::MeshManager::Adapter::OwnedBoundaryNodeIterator bnd = masterMeshAdapter->beginOwnedBoundary( bndIds[i] );
    AMP::Mesh::MeshManager::Adapter::OwnedBoundaryNodeIterator end_bnd = masterMeshAdapter->endOwnedBoundary( bndIds[i] );

    for( ; bnd != end_bnd; ++bnd) {
      std::vector<unsigned int> bndGlobalIds;
      master_dof_map->getDOFs(*bnd, bndGlobalIds, bndDofIds[i]);

      for(size_t j = 0; j < bndGlobalIds.size(); j++) {
        masterVec->setLocalValueByGlobalID(bndGlobalIds[j], 0.0);
      }//end for j
    }//end for bnd
  }//end for i

  vec->makeConsistent( AMP::LinearAlgebra::Vector::CONSISTENT_SET );
}

void setSlaveToZero(AMP::Mesh::MeshManager::Adapter::shared_ptr slaveMeshAdapter,
    AMP::LinearAlgebra::Variable::shared_ptr slaveVar,
    std::vector<size_t> const & slaveNodes,
    AMP::LinearAlgebra::Vector::shared_ptr vec) {
  AMP::Mesh::DOFMap::shared_ptr slave_dof_map = slaveMeshAdapter->getDOFMap(slaveVar);

  AMP::LinearAlgebra::Vector::shared_ptr slaveVec = vec->subsetVectorForVariable(slaveVar);

  std::vector<unsigned int> dofs(3);
  dofs[0] = 0; dofs[1] = 1; dofs[2] = 2;

  for(size_t i = 0; i < slaveNodes.size(); i++) {
    AMP::Mesh::LibMeshNode slaveNd =  slaveMeshAdapter->getNode( slaveNodes[i] );
    std::vector<unsigned int> slaveGlobalIds;
    slave_dof_map->getDOFs(slaveNd, slaveGlobalIds, dofs);
    for(size_t j = 0; j < dofs.size(); j++) {
      slaveVec->setLocalValueByGlobalID(slaveGlobalIds[j], 0.0);
    }//end for j
  }//end for i

  vec->makeConsistent( AMP::LinearAlgebra::Vector::CONSISTENT_SET );
}

void masterToSlaveCorrection(AMP::Mesh::MeshManager::Adapter::shared_ptr slaveMeshAdapter,
    AMP::Mesh::MeshManager::Adapter::shared_ptr masterMeshAdapter,
    AMP::LinearAlgebra::Variable::shared_ptr slaveVar,
    AMP::LinearAlgebra::Variable::shared_ptr masterVar,
    std::vector<size_t> const & slaveNodes,
    std::vector<std::vector<size_t> > const & slave2MasterNodes,
    std::vector<std::vector<double> > const & slave2MasterFactors,
    AMP::LinearAlgebra::Vector::shared_ptr vec) {
  AMP::Mesh::DOFMap::shared_ptr slave_dof_map = slaveMeshAdapter->getDOFMap(slaveVar);

  AMP::Mesh::DOFMap::shared_ptr master_dof_map = masterMeshAdapter->getDOFMap(masterVar);

  AMP::LinearAlgebra::Vector::shared_ptr slaveVec = vec->subsetVectorForVariable(slaveVar);

  AMP::LinearAlgebra::Vector::shared_ptr masterVec = vec->subsetVectorForVariable(masterVar);

  std::vector<unsigned int> dofs(3);
  dofs[0] = 0; dofs[1] = 1; dofs[2] = 2;

  for(size_t i = 0; i < slaveNodes.size(); i++) {
    AMP::Mesh::LibMeshNode slaveNd =  slaveMeshAdapter->getNode( slaveNodes[i] );
    std::vector<unsigned int> slaveGlobalIds;
    slave_dof_map->getDOFs(slaveNd, slaveGlobalIds, dofs);
    double slaveVal[] = {0.0, 0.0, 0.0};
    for(size_t j = 0; j < slave2MasterNodes[i].size(); j++) {
      AMP::Mesh::LibMeshNode masterNd =  masterMeshAdapter->getNode( slave2MasterNodes[i][j] );
      std::vector<unsigned int> masterGlobalIds;
      master_dof_map->getDOFs(masterNd, masterGlobalIds, dofs);
      for(size_t k = 0; k < dofs.size(); k++) {
        double masterVal = masterVec->getValueByGlobalID(masterGlobalIds[k]);
        slaveVal[k] += (slave2MasterFactors[i][j]*masterVal);
      }//end for k
    }//end for j
    for(size_t k = 0; k < dofs.size(); k++) {
      slaveVec->setLocalValueByGlobalID(slaveGlobalIds[k], slaveVal[k]);
    }//end for k
  }//end for i

  vec->makeConsistent( AMP::LinearAlgebra::Vector::CONSISTENT_SET );
}

void slaveToMasterCorrection(AMP::Mesh::MeshManager::Adapter::shared_ptr slaveMeshAdapter,
    AMP::Mesh::MeshManager::Adapter::shared_ptr masterMeshAdapter,
    AMP::LinearAlgebra::Variable::shared_ptr slaveVar,
    AMP::LinearAlgebra::Variable::shared_ptr masterVar,
    std::vector<size_t> const & slaveNodes,
    std::vector<std::vector<size_t> > const & slave2MasterNodes,
    std::vector<std::vector<double> > const & slave2MasterFactors,
    AMP::LinearAlgebra::Vector::shared_ptr vec) {
  AMP::Mesh::DOFMap::shared_ptr slave_dof_map = slaveMeshAdapter->getDOFMap(slaveVar);

  AMP::Mesh::DOFMap::shared_ptr master_dof_map = masterMeshAdapter->getDOFMap(masterVar);

  AMP::LinearAlgebra::Vector::shared_ptr slaveVec = vec->subsetVectorForVariable(slaveVar);

  AMP::LinearAlgebra::Vector::shared_ptr masterVec = vec->subsetVectorForVariable(masterVar);

  std::vector<unsigned int> dofs(3);
  dofs[0] = 0; dofs[1] = 1; dofs[2] = 2;

  for(size_t i = 0; i < slaveNodes.size(); i++) {
    AMP::Mesh::LibMeshNode slaveNd =  slaveMeshAdapter->getNode( slaveNodes[i] );
    std::vector<unsigned int> slaveGlobalIds;
    slave_dof_map->getDOFs(slaveNd, slaveGlobalIds, dofs);
    for(size_t j = 0; j < slave2MasterNodes[i].size(); j++) {
      AMP::Mesh::LibMeshNode masterNd =  masterMeshAdapter->getNode( slave2MasterNodes[i][j] );
      std::vector<unsigned int> masterGlobalIds;
      master_dof_map->getDOFs(masterNd, masterGlobalIds, dofs);
      for(size_t k = 0; k < dofs.size(); k++) {
        double slaveVal = slaveVec->getValueByGlobalID(slaveGlobalIds[k]);
        double masterVal = (slave2MasterFactors[i][j]*slaveVal);
        masterVec->addValueByGlobalID(masterGlobalIds[k], masterVal);
      }//end for k
    }//end for j
  }//end for i

  vec->makeConsistent( AMP::LinearAlgebra::Vector::CONSISTENT_ADD );
}

#endif

void myTest(AMP::UnitTest *ut, std::string exeName) {
#if 0
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

  AMP::Mesh::MeshManagerParameters::shared_ptr meshmgrParams ( new AMP::Mesh::MeshManagerParameters ( input_db ) );
  AMP::Mesh::MeshManager::shared_ptr manager ( new AMP::Mesh::MeshManager ( meshmgrParams ) );
  AMP::Mesh::MeshManager::Adapter::shared_ptr masterMeshAdapter = manager->getMesh ( masterMeshName );
  AMP::Mesh::MeshManager::Adapter::shared_ptr slaveMeshAdapter = manager->getMesh ( slaveMeshName );

  std::cout<<"Master volume has "<<(masterMeshAdapter->numLocalNodes())<<" nodes."<<std::endl;
  std::cout<<"Slave volume has "<<(slaveMeshAdapter->numLocalNodes())<<" nodes."<<std::endl;

  double minXYZ[3];
  double maxXYZ[3];

  computeRGboundingBox(precision, masterMeshAdapter, minXYZ, maxXYZ);

  double rgH[3];
  std::vector<std::vector<size_t> > rg2ElemMap;

  computeRG2ElemMap(precision, rgDim, masterMeshAdapter, minXYZ, maxXYZ, rg2ElemMap, rgH);

  std::vector<size_t> slaveNodes;
  std::vector<size_t> slave2MasterElem;

  computeSlave2MasterElem(slaveId, masterId, precision, rgH, rgDim, 
      slaveMeshAdapter, masterMeshAdapter, minXYZ, maxXYZ, rg2ElemMap, slaveNodes, slave2MasterElem);

  std::cout<<"# slave contact nodes = "<<(slaveNodes.size())<<std::endl;

  std::vector<std::vector<size_t> > slave2MasterNodes;
  std::vector<std::vector<double> > slave2MasterFactors;

  computeSlave2MasterNodes(precision, slaveId, masterId, slaveMeshAdapter, masterMeshAdapter,
      slaveNodes, slave2MasterElem, slave2MasterNodes, slave2MasterFactors);

  boost::shared_ptr<AMP::Operator::ElementPhysicsModel> masterElementPhysicsModel;
  boost::shared_ptr<AMP::Operator::MechanicsLinearFEOperator> masterOperator = boost::dynamic_pointer_cast<
    AMP::Operator::MechanicsLinearFEOperator>(AMP::Operator::OperatorBuilder::createOperator(
          masterMeshAdapter, "MasterVolumeOperator", input_db, masterElementPhysicsModel));

  boost::shared_ptr<AMP::Operator::ElementPhysicsModel> slaveElementPhysicsModel;
  boost::shared_ptr<AMP::Operator::MechanicsLinearFEOperator> slaveOperator = boost::dynamic_pointer_cast<
    AMP::Operator::MechanicsLinearFEOperator>(AMP::Operator::OperatorBuilder::createOperator(
          slaveMeshAdapter, "SlaveVolumeOperator", input_db, slaveElementPhysicsModel));

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
          slaveMeshAdapter, "SlaveLoadOperator", input_db, dummyModel));
  loadOperator->setVariable(slaveVar);

  AMP::LinearAlgebra::Vector::shared_ptr nullVec;
  AMP::LinearAlgebra::Vector::shared_ptr columnSolVec = manager->createVector(columnVar);
  AMP::LinearAlgebra::Vector::shared_ptr columnRhsVec = columnSolVec->cloneVector();
  AMP::LinearAlgebra::Vector::shared_ptr columnResVec = columnSolVec->cloneVector();

  columnRhsVec->zero();

  loadOperator->apply(nullVec, nullVec, columnRhsVec, 1.0, 0.0);

  boost::shared_ptr<AMP::Database> bndDatabase = input_db->getDatabase("MasterDirichletBoundary");
  std::vector<unsigned int> bndIds(bndDatabase->getInteger("number_of_ids"));
  std::vector<std::vector<unsigned int> > bndDofIds(bndIds.size());
  std::vector<std::vector<double> > bndDofVals(bndIds.size());
  for(unsigned int i = 0; i < bndIds.size(); i++) {
    char tmp[200];
    sprintf(tmp, "id_%u", i);
    bndIds[i] = bndDatabase->getInteger(tmp);
    sprintf(tmp, "number_of_dofs_%u", i); 
    bndDofIds[i].resize(bndDatabase->getInteger(tmp));
    bndDofVals[i].resize(bndDofIds[i].size());
    for(unsigned int j = 0; j < bndDofIds[i].size(); j++) {
      char key[200];
      sprintf(key, "dof_%u_%u", i, j);
      bndDofIds[i][j] = bndDatabase->getInteger(key);
      sprintf(key, "value_%u_%u", i, j);
      bndDofVals[i][j] = bndDatabase->getDouble(key);
    }//end for j
  }//end for i

  columnSolVec->zero();

  setDirichletVals(masterMeshAdapter, bndIds, bndDofIds, 
      bndDofVals, masterVar, columnSolVec);

  masterToSlaveCorrection(slaveMeshAdapter, masterMeshAdapter, slaveVar, masterVar,
      slaveNodes, slave2MasterNodes, slave2MasterFactors, columnSolVec);

  columnOperator->apply(nullVec, columnSolVec, columnResVec, 1.0, 0.0);

  columnRhsVec->subtract(columnRhsVec, columnResVec);

  int maxCGiters = input_db->getInteger("maxCGiters");

  {
    std::cout<<"MPC-MLCG: "<<std::endl;
    AMP::LinearAlgebra::Vector::shared_ptr MatOutVec = columnSolVec->cloneVector();
    AMP::LinearAlgebra::Vector::shared_ptr pVec = columnSolVec->cloneVector();
    AMP::LinearAlgebra::Vector::shared_ptr zVec = columnSolVec->cloneVector();

    columnSolVec->zero();

    columnOperator->apply(nullVec, columnSolVec, MatOutVec, 1.0, 0.0);
    columnResVec->subtract(columnRhsVec, MatOutVec);
    slaveToMasterCorrection(slaveMeshAdapter, masterMeshAdapter, slaveVar, masterVar,
        slaveNodes, slave2MasterNodes, slave2MasterFactors, columnResVec);
    setSlaveToZero(slaveMeshAdapter, slaveVar, slaveNodes, columnResVec);
    setDirichletDofsToZero(masterMeshAdapter, bndIds, bndDofIds, masterVar, columnResVec);

    pcNoneApply(columnResVec, zVec);
    setDirichletDofsToZero(masterMeshAdapter, bndIds, bndDofIds, masterVar, zVec);
    masterToSlaveCorrection(slaveMeshAdapter, masterMeshAdapter, slaveVar, masterVar,
        slaveNodes, slave2MasterNodes, slave2MasterFactors, zVec);

    pVec->copyVector(zVec);

    for(int iter = 0; iter <= maxCGiters; iter++) {
      double resNorm = columnResVec->L2Norm();

      std::cout<<"CG-Iter = "<<iter<<" ResNorm2 = "<<std::setprecision(15)<<resNorm<<std::endl;

      columnOperator->apply(nullVec, pVec, MatOutVec, 1.0, 0.0);
      slaveToMasterCorrection(slaveMeshAdapter, masterMeshAdapter, slaveVar, masterVar,
          slaveNodes, slave2MasterNodes, slave2MasterFactors, MatOutVec);
      setSlaveToZero(slaveMeshAdapter, slaveVar, slaveNodes, MatOutVec);
      setDirichletDofsToZero(masterMeshAdapter, bndIds, bndDofIds, masterVar, MatOutVec);

      double resOldDotZ = columnResVec->dot(zVec);

      double alphaDenom = MatOutVec->dot(pVec);

      double alpha = resOldDotZ/alphaDenom;

      columnSolVec->axpy(alpha, pVec, columnSolVec);

      columnResVec->axpy(-alpha, MatOutVec, columnResVec);

      pcNoneApply(columnResVec, zVec);
      setDirichletDofsToZero(masterMeshAdapter, bndIds, bndDofIds, masterVar, zVec);
      masterToSlaveCorrection(slaveMeshAdapter, masterMeshAdapter, slaveVar, masterVar,
          slaveNodes, slave2MasterNodes, slave2MasterFactors, zVec);

      double resNewDotZ = columnResVec->dot(zVec);

      double beta = resNewDotZ/resOldDotZ;

      pVec->axpy(beta, pVec, zVec);
    }//end for iter
  }

  setDirichletVals(masterMeshAdapter, bndIds, bndDofIds, 
      bndDofVals, masterVar, columnSolVec);

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

  std::string exeName = "testContactV5";

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


