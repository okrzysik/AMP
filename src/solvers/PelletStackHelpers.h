
#ifndef included_AMP_PelletStackHelpers
#define included_AMP_PelletStackHelpers

#include "operators/OperatorBuilder.h"
#include "operators/LinearBVPOperator.h"
#include "operators/NonlinearBVPOperator.h"
#include "operators/boundary/DirichletVectorCorrection.h"
#include "operators/PelletStackOperator.h"
#include "operators/CoupledOperator.h"
#include "operators/map/NodeToNodeMap.h"
#include "operators/map/AsyncMapColumnOperator.h"
#include "operators/mechanics/MechanicsNonlinearFEOperator.h"

#include "ampmesh/MeshManager.h"
#include "ampmesh/MeshAdapter.h"
#include "vectors/CommCollectVector.h"

#include "solvers/PetscKrylovSolver.h"
#include "solvers/ColumnSolver.h"
#include "solvers/PelletStackMechanicsSolver.h"
#include "solvers/TrilinosMLSolver.h"

void helperCreateStackOperatorForPelletMechanics(AMP::Mesh::MeshManager::shared_ptr manager,
    boost::shared_ptr<AMP::Operator::AsyncMapColumnOperator> n2nmaps, AMP::AMP_MPI globalComm,
    boost::shared_ptr<AMP::InputDatabase> global_input_db, 
    boost::shared_ptr<AMP::Operator::PelletStackOperator> & pelletStackOp)
{
  boost::shared_ptr<AMP::Database> pelletStackOp_db = global_input_db->getDatabase("PelletStackOperator");
  pelletStackOp_db->putInteger("TOTAL_NUMBER_OF_PELLETS", global_input_db->getInteger("NumberOfPelletMeshes"));
  boost::shared_ptr<AMP::Operator::PelletStackOperatorParameters> pelletStackOpParams(new 
      AMP::Operator::PelletStackOperatorParameters(pelletStackOp_db));
  pelletStackOpParams->d_pelletStackComm = globalComm;
  pelletStackOpParams->d_n2nMaps = n2nmaps;
  pelletStackOpParams->d_meshManager = manager;
  pelletStackOp.reset(new AMP::Operator::PelletStackOperator(pelletStackOpParams));
}

void helperCreateColumnOperatorsForPelletMechanics(std::vector<unsigned int> localPelletIds, 
    std::vector<AMP::Mesh::MeshManager::Adapter::shared_ptr> localMeshes,
    boost::shared_ptr<AMP::InputDatabase> global_input_db,
    boost::shared_ptr<AMP::Operator::ColumnOperator> & nonlinearColumnOperator,
    boost::shared_ptr<AMP::Operator::ColumnOperator> & linearColumnOperator) 
{
  boost::shared_ptr<AMP::Operator::OperatorParameters> emptyParams;
  nonlinearColumnOperator.reset(new AMP::Operator::ColumnOperator(emptyParams));
  linearColumnOperator.reset(new AMP::Operator::ColumnOperator(emptyParams));
  for(unsigned int id = 0; id < localPelletIds.size(); id++) {
    std::string prefix = "";
    if(localPelletIds[id] == 0)
    {
      prefix = "Bottom";
    }

    AMP::Mesh::MeshManager::Adapter::shared_ptr meshAdapter = localMeshes[id];

    boost::shared_ptr<AMP::Operator::ElementPhysicsModel> mechModel;
    boost::shared_ptr<AMP::Operator::NonlinearBVPOperator> nonlinOperator =
      boost::dynamic_pointer_cast<AMP::Operator::NonlinearBVPOperator>(AMP::Operator::OperatorBuilder::createOperator(meshAdapter,
            prefix+"PelletMechanicsNonlinearBVPOperator", global_input_db, mechModel));
    nonlinearColumnOperator->append(nonlinOperator);

    boost::shared_ptr<AMP::Operator::LinearBVPOperator> linOperator =
      boost::dynamic_pointer_cast<AMP::Operator::LinearBVPOperator>(AMP::Operator::OperatorBuilder::createOperator(meshAdapter,
            prefix+"PelletMechanicsLinearBVPOperator", global_input_db, mechModel));
    linearColumnOperator->append(linOperator);
  }//end for id
} 

void helperCreateCoupledOperatorForPelletMechanics(boost::shared_ptr<AMP::Operator::AsyncMapColumnOperator> n2nmaps, 
						   boost::shared_ptr<AMP::Operator::ColumnOperator> nonlinearColumnOperator,
						   boost::shared_ptr<AMP::Operator::CoupledOperator> & coupledOp) 
{
  boost::shared_ptr<AMP::InputDatabase> emptyDb;
  boost::shared_ptr<AMP::Operator::CoupledOperatorParameters> coupledOpParams(new
      AMP::Operator::CoupledOperatorParameters(emptyDb));
  coupledOpParams->d_MapOperator = n2nmaps;
  coupledOpParams->d_BVPOperator = nonlinearColumnOperator;
  coupledOp.reset(new AMP::Operator::CoupledOperator(coupledOpParams));
}

void helperSetFrozenVectorForMapsForPelletMechanics(AMP::Mesh::MeshManager::shared_ptr manager, 
    boost::shared_ptr<AMP::Operator::CoupledOperator> coupledOp) 
{
  boost::shared_ptr<AMP::Operator::AsyncMapColumnOperator> n2nmaps = 
    boost::dynamic_pointer_cast<AMP::Operator::AsyncMapColumnOperator>(coupledOp->getOperator(1)); 
  boost::shared_ptr<AMP::Operator::ColumnOperator> nonlinearColumnOperator = 
    boost::dynamic_pointer_cast<AMP::Operator::ColumnOperator>(coupledOp->getOperator(2));
  AMP::LinearAlgebra::Variable::shared_ptr dispVar = nonlinearColumnOperator->getOutputVariable();
  AMP::LinearAlgebra::Vector::shared_ptr dirichletValues =  manager->createVector ( dispVar );
  n2nmaps->setVector(dirichletValues);
  for(int id = 0; id < nonlinearColumnOperator->getNumberOfOperators(); id++) {
    boost::shared_ptr<AMP::Operator::DirichletVectorCorrection> dirichletOp =
      boost::dynamic_pointer_cast<AMP::Operator::DirichletVectorCorrection>(
          boost::dynamic_pointer_cast<AMP::Operator::NonlinearBVPOperator>(
            nonlinearColumnOperator->getOperator(id))->getBoundaryOperator());
    dirichletOp->setDirichletValues(dirichletValues);
  }//end for id
}

void helperCreateAllOperatorsForPelletMechanics(AMP::Mesh::MeshManager::shared_ptr manager,
    AMP::AMP_MPI globalComm, boost::shared_ptr<AMP::InputDatabase> global_input_db,
    boost::shared_ptr<AMP::Operator::CoupledOperator> & coupledOp,
    boost::shared_ptr<AMP::Operator::ColumnOperator> & linearColumnOperator, 
    boost::shared_ptr<AMP::Operator::PelletStackOperator> & pelletStackOp)
{
  boost::shared_ptr<AMP::Operator::AsyncMapColumnOperator> n2nmaps =
    AMP::Operator::AsyncMapColumnOperator::build<AMP::Operator::NodeToNodeMap> ( manager , global_input_db );

  helperCreateStackOperatorForPelletMechanics(manager, n2nmaps, globalComm, global_input_db, pelletStackOp);

  std::vector<unsigned int> localPelletIds = pelletStackOp->getLocalPelletIds();
  std::vector<AMP::Mesh::MeshManager::Adapter::shared_ptr> localMeshes = pelletStackOp->getLocalMeshes();

  boost::shared_ptr<AMP::Operator::ColumnOperator> nonlinearColumnOperator;
  helperCreateColumnOperatorsForPelletMechanics(localPelletIds, localMeshes, global_input_db, nonlinearColumnOperator, linearColumnOperator);
  helperCreateCoupledOperatorForPelletMechanics(n2nmaps, nonlinearColumnOperator, coupledOp);
  helperSetFrozenVectorForMapsForPelletMechanics(manager, coupledOp);
}

void helperCreateVectorsForPelletMechanics(AMP::Mesh::MeshManager::shared_ptr manager,
    AMP::Operator::Operator::shared_ptr coupledOp, AMP::AMP_MPI globalComm, 
    AMP::LinearAlgebra::Vector::shared_ptr & solVec,
    AMP::LinearAlgebra::Vector::shared_ptr & rhsVec,
    AMP::LinearAlgebra::Vector::shared_ptr & scaledRhsVec) 
{
  AMP::LinearAlgebra::Variable::shared_ptr dispVar = coupledOp->getOutputVariable();
  solVec = AMP::LinearAlgebra::CommCollectVector::view (
      manager->createVector ( dispVar ) , globalComm );
  rhsVec = AMP::LinearAlgebra::CommCollectVector::view (
      manager->createVector ( dispVar ) , globalComm );
  scaledRhsVec = AMP::LinearAlgebra::CommCollectVector::view (
      manager->createVector ( dispVar ) , globalComm );
}

void helperBuildPointLoadRHSForPelletMechanics(boost::shared_ptr<AMP::InputDatabase> global_input_db, 
    boost::shared_ptr<AMP::Operator::CoupledOperator> coupledOp,
    AMP::LinearAlgebra::Vector::shared_ptr rhsVec) 
{
  AMP::LinearAlgebra::Vector::shared_ptr nullVec;
  boost::shared_ptr<AMP::Operator::ColumnOperator> nonlinearColumnOperator = 
    boost::dynamic_pointer_cast<AMP::Operator::ColumnOperator>(coupledOp->getOperator(2));
  rhsVec->zero();
  for(int id = 0; id < nonlinearColumnOperator->getNumberOfOperators(); id++) {
    AMP::Operator::Operator::shared_ptr currOp = nonlinearColumnOperator->getOperator(id);
    AMP::Mesh::MeshManager::Adapter::shared_ptr meshAdapter = currOp->getMeshAdapter();
    boost::shared_ptr<AMP::Operator::ElementPhysicsModel> dummyModel;
    boost::shared_ptr<AMP::Operator::DirichletVectorCorrection> loadOp = 
      boost::dynamic_pointer_cast<AMP::Operator::DirichletVectorCorrection>(
          AMP::Operator::OperatorBuilder::createOperator(meshAdapter,
            "PointLoad", global_input_db, dummyModel));
    loadOp->setVariable(currOp->getOutputVariable());
    loadOp->apply(nullVec, nullVec, rhsVec, 1.0, 0.0);
  }//end for id
}

void helperApplyBoundaryCorrectionsForPelletMechanics(boost::shared_ptr<AMP::Operator::CoupledOperator> coupledOp,
    AMP::LinearAlgebra::Vector::shared_ptr solVec, 
    AMP::LinearAlgebra::Vector::shared_ptr rhsVec) 
{
  boost::shared_ptr<AMP::Operator::ColumnOperator> nonlinearColumnOperator = 
    boost::dynamic_pointer_cast<AMP::Operator::ColumnOperator>(coupledOp->getOperator(2));
  for(int id = 0; id < nonlinearColumnOperator->getNumberOfOperators(); id++) {
    boost::shared_ptr<AMP::Operator::NonlinearBVPOperator> nonlinOperator =
      boost::dynamic_pointer_cast<AMP::Operator::NonlinearBVPOperator>(
          nonlinearColumnOperator->getOperator(id));
    nonlinOperator->modifyInitialSolutionVector(solVec);
    nonlinOperator->modifyRHSvector(rhsVec);
  }//end for id
}

void helperCreateTemperatureVectorsForPelletMechanics(AMP::Mesh::MeshManager::shared_ptr manager, 
    boost::shared_ptr<AMP::Operator::CoupledOperator> coupledOp,
    AMP::LinearAlgebra::Vector::shared_ptr & initialTemperatureVec,
    AMP::LinearAlgebra::Vector::shared_ptr & finalTemperatureVec) 
{
  boost::shared_ptr<AMP::Operator::ColumnOperator> nonlinearColumnOperator = 
    boost::dynamic_pointer_cast<AMP::Operator::ColumnOperator>(coupledOp->getOperator(2));
  boost::shared_ptr<AMP::LinearAlgebra::MultiVariable> tempVar( new AMP::LinearAlgebra::MultiVariable("ColumnVariable"));
  for(int id = 0; id < nonlinearColumnOperator->getNumberOfOperators(); id++) {
    AMP::LinearAlgebra::Variable::shared_ptr opVar = 
      (nonlinearColumnOperator->getOperator(id))->getInputVariable(AMP::Operator::Mechanics::TEMPERATURE);
    if(opVar.get() != NULL) {
      tempVar->add(opVar);
    }
  }//end for id
  initialTemperatureVec = manager->createVector(tempVar);
  finalTemperatureVec = initialTemperatureVec->cloneVector();
}

void helperSetReferenceTemperatureForPelletMechanics(boost::shared_ptr<AMP::Operator::CoupledOperator> coupledOp,
    AMP::LinearAlgebra::Vector::shared_ptr initialTemperatureVec) 
{
  boost::shared_ptr<AMP::Operator::ColumnOperator> nonlinearColumnOperator = 
    boost::dynamic_pointer_cast<AMP::Operator::ColumnOperator>(coupledOp->getOperator(2));
  for(int id = 0; id < nonlinearColumnOperator->getNumberOfOperators(); id++) {
    boost::shared_ptr<AMP::Operator::NonlinearBVPOperator> bvpOp = boost::dynamic_pointer_cast<
      AMP::Operator::NonlinearBVPOperator>(nonlinearColumnOperator->getOperator(id));
    boost::shared_ptr<AMP::Operator::MechanicsNonlinearFEOperator> mechOp = boost::dynamic_pointer_cast<
      AMP::Operator::MechanicsNonlinearFEOperator>(bvpOp->getVolumeOperator());
    mechOp->setReferenceTemperature(initialTemperatureVec);
  }//end for id
}

void helperSetFinalTemperatureForPelletMechanics(boost::shared_ptr<AMP::Operator::CoupledOperator> coupledOp,
    AMP::LinearAlgebra::Vector::shared_ptr finalTemperatureVec) 
{
  boost::shared_ptr<AMP::Operator::ColumnOperator> nonlinearColumnOperator = 
    boost::dynamic_pointer_cast<AMP::Operator::ColumnOperator>(coupledOp->getOperator(2));
  for(int id = 0; id < nonlinearColumnOperator->getNumberOfOperators(); id++) {
    boost::shared_ptr<AMP::Operator::NonlinearBVPOperator> bvpOp = boost::dynamic_pointer_cast<
      AMP::Operator::NonlinearBVPOperator>(nonlinearColumnOperator->getOperator(id));
    boost::shared_ptr<AMP::Operator::MechanicsNonlinearFEOperator> mechOp = boost::dynamic_pointer_cast<
      AMP::Operator::MechanicsNonlinearFEOperator>(bvpOp->getVolumeOperator());
    mechOp->setVector(AMP::Operator::Mechanics::TEMPERATURE, finalTemperatureVec); 
  }//end for id
}

void helperBuildColumnSolverForPelletMechanics(boost::shared_ptr<AMP::Database> columnSolver_db,
    boost::shared_ptr<AMP::Operator::ColumnOperator> linearColumnOperator,
    boost::shared_ptr<AMP::Solver::ColumnSolver> & columnSolver) 
{
  boost::shared_ptr<AMP::Database> ikspSolver_db = columnSolver_db->getDatabase("KrylovSolver");
  boost::shared_ptr<AMP::Database> mlSolver_db = ikspSolver_db->getDatabase("MLSolver");
  boost::shared_ptr<AMP::Solver::SolverStrategyParameters> columnSolverParams(new
      AMP::Solver::SolverStrategyParameters(columnSolver_db));
  columnSolverParams->d_pOperator = linearColumnOperator;
  columnSolver.reset(new AMP::Solver::ColumnSolver(columnSolverParams));
  for(int id = 0; id < linearColumnOperator->getNumberOfOperators(); id++) {
    AMP::Operator::Operator::shared_ptr currOp = linearColumnOperator->getOperator(id);
    boost::shared_ptr<AMP::Solver::SolverStrategyParameters> mlSolverParams(new
        AMP::Solver::SolverStrategyParameters(mlSolver_db));
    mlSolverParams->d_pOperator = currOp;
    boost::shared_ptr<AMP::Solver::TrilinosMLSolver> mlSolver(new
        AMP::Solver::TrilinosMLSolver(mlSolverParams));

    boost::shared_ptr<AMP::Solver::PetscKrylovSolverParameters> ikspSolverParams(new
        AMP::Solver::PetscKrylovSolverParameters(ikspSolver_db));
    ikspSolverParams->d_pOperator = currOp;
    ikspSolverParams->d_comm = (currOp->getMeshAdapter())->getComm();
    ikspSolverParams->d_pPreconditioner = mlSolver;
    boost::shared_ptr<AMP::Solver::PetscKrylovSolver> ikspSolver(new
        AMP::Solver::PetscKrylovSolver(ikspSolverParams));

    columnSolver->append(ikspSolver);
  }//end for id
}

void helperBuildStackSolverForPelletMechanics(boost::shared_ptr<AMP::Database> pelletStackSolver_db, 
    boost::shared_ptr<AMP::Operator::PelletStackOperator> pelletStackOp, 
    boost::shared_ptr<AMP::Operator::ColumnOperator> linearColumnOperator, 
    boost::shared_ptr<AMP::Solver::SolverStrategy> & pelletStackSolver) 
{
  boost::shared_ptr<AMP::Database> columnSolver_db = pelletStackSolver_db->getDatabase("ColumnSolver");
  boost::shared_ptr<AMP::Solver::ColumnSolver> columnSolver;
  helperBuildColumnSolverForPelletMechanics(columnSolver_db, linearColumnOperator, columnSolver);
  boost::shared_ptr<AMP::Solver::PelletStackMechanicsSolverParameters> pelletStackSolverParams(new
      AMP::Solver::PelletStackMechanicsSolverParameters(pelletStackSolver_db));
  pelletStackSolverParams->d_columnSolver = columnSolver;
  pelletStackSolverParams->d_pOperator = pelletStackOp;
  pelletStackSolver.reset(new AMP::Solver::PelletStackMechanicsSolver(pelletStackSolverParams));
}

void helperResetNonlinearOperatorForPelletMechanics(boost::shared_ptr<AMP::Operator::CoupledOperator> coupledOp) 
{
  boost::shared_ptr<AMP::Operator::ColumnOperator> nonlinearColumnOperator =
    boost::dynamic_pointer_cast<AMP::Operator::ColumnOperator>(coupledOp->getOperator(2));
  boost::shared_ptr<AMP::InputDatabase> tmp_db (new AMP::InputDatabase("Dummy"));
  boost::shared_ptr<AMP::Operator::MechanicsNonlinearFEOperatorParameters> tmpParams(new
      AMP::Operator::MechanicsNonlinearFEOperatorParameters(tmp_db));

  for(int id = 0; id < nonlinearColumnOperator->getNumberOfOperators(); id++) {
    boost::shared_ptr<AMP::Operator::NonlinearBVPOperator> nonlinOperator =
      boost::dynamic_pointer_cast<AMP::Operator::NonlinearBVPOperator>(
          nonlinearColumnOperator->getOperator(id));
    (nonlinOperator->getVolumeOperator())->reset(tmpParams);
  }//end for id
}

#endif




