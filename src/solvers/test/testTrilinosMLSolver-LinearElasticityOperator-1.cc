
#include <iostream>
#include <string>
#include <cassert>

#include "utils/AMPManager.h"
#include "utils/UnitTest.h"
#include "utils/Utilities.h"

#include "discretization/simpleDOF_Manager.h"
#include "vectors/VectorBuilder.h"
#include "ampmesh/SiloIO.h"

/* AMP files */
#include "utils/Database.h"
#include "utils/InputDatabase.h"
#include "utils/InputManager.h"
#include "utils/AMP_MPI.h"
#include "utils/AMPManager.h"
#include "utils/PIO.h"

#include "operators/boundary/DirichletVectorCorrection.h"
#include "operators/LinearBVPOperator.h"
#include "operators/OperatorBuilder.h"

#include "solvers/TrilinosMLSolver.h"

void linearElasticTest(AMP::UnitTest *ut )
{
  std::string exeName("testTrilinosMLSolver-LinearElasticityOperator-1");
  std::string input_file = "input_" + exeName;
  std::string log_file = "output_" + exeName;

  AMP::PIO::logOnlyNodeZero(log_file);

#ifdef USES_SILO
  // Create the silo writer and register the data
  AMP::Mesh::SiloIO::shared_ptr siloWriter( new AMP::Mesh::SiloIO);
#endif

  boost::shared_ptr<AMP::InputDatabase> input_db(new AMP::InputDatabase("input_db"));
  AMP::InputManager::getManager()->parseInputFile(input_file, input_db);
  input_db->printClassData(AMP::plog);

  AMP_INSIST(input_db->keyExists("Mesh"), "Key ''Mesh'' is missing!");
  boost::shared_ptr<AMP::Database> mesh_db = input_db->getDatabase("Mesh");
  boost::shared_ptr<AMP::Mesh::MeshParameters> meshParams(new AMP::Mesh::MeshParameters(mesh_db));
  meshParams->setComm(AMP::AMP_MPI(AMP_COMM_WORLD));
  AMP::Mesh::Mesh::shared_ptr meshAdapter = AMP::Mesh::Mesh::buildMesh(meshParams);

  boost::shared_ptr<AMP::Operator::ElementPhysicsModel> elementPhysicsModel;
  boost::shared_ptr<AMP::Operator::LinearBVPOperator> bvpOperator =
    boost::dynamic_pointer_cast<AMP::Operator::LinearBVPOperator>(AMP::Operator::OperatorBuilder::createOperator(meshAdapter,
          "MechanicsBVPOperator", input_db, elementPhysicsModel));

  AMP::LinearAlgebra::Variable::shared_ptr var = bvpOperator->getOutputVariable();

  boost::shared_ptr<AMP::Operator::ElementPhysicsModel> dummyModel;
  boost::shared_ptr<AMP::Operator::DirichletVectorCorrection> dirichletVecOp =
    boost::dynamic_pointer_cast<AMP::Operator::DirichletVectorCorrection>(AMP::Operator::OperatorBuilder::createOperator(meshAdapter,
          "Load_Boundary", input_db, dummyModel));
  //This has an in-place apply. So, it has an empty input variable and
  //the output variable is the same as what it is operating on. 
  dirichletVecOp->setVariable(var);

  AMP::Discretization::DOFManager::shared_ptr dofMap = AMP::Discretization::simpleDOFManager::create(
      meshAdapter, AMP::Mesh::Vertex, 1, 3, true); 

  AMP::LinearAlgebra::Vector::shared_ptr nullVec;
  AMP::LinearAlgebra::Vector::shared_ptr mechSolVec = AMP::LinearAlgebra::createVector(dofMap, var, true);
  AMP::LinearAlgebra::Vector::shared_ptr mechRhsVec = mechSolVec->cloneVector();
  AMP::LinearAlgebra::Vector::shared_ptr mechResVec = mechSolVec->cloneVector();

  mechSolVec->setToScalar(0.5);
  mechRhsVec->setToScalar(0.0);
  mechResVec->setToScalar(0.0);

  dirichletVecOp->apply(nullVec, nullVec, mechRhsVec, 1.0, 0.0);

  double rhsNorm = mechRhsVec->L2Norm();

  std::cout<<"RHS Norm: "<<rhsNorm<<std::endl;

  double initSolNorm = mechSolVec->L2Norm();

  std::cout<<"Initial Solution Norm: "<<initSolNorm<<std::endl;

  bvpOperator->apply(mechRhsVec, mechSolVec, mechResVec, 1.0, -1.0);

  double initResidualNorm = mechResVec->L2Norm();

  std::cout<<"Initial Residual Norm: "<<initResidualNorm<<std::endl;

  boost::shared_ptr<AMP::Database> mlSolver_db = input_db->getDatabase("LinearSolver"); 

  boost::shared_ptr<AMP::Solver::SolverStrategyParameters> mlSolverParams(new AMP::Solver::SolverStrategyParameters(mlSolver_db));

  mlSolverParams->d_pOperator = bvpOperator;

  // create the ML solver interface
  boost::shared_ptr<AMP::Solver::TrilinosMLSolver> mlSolver(new AMP::Solver::TrilinosMLSolver(mlSolverParams));

  mlSolver->setZeroInitialGuess(false);

  mlSolver->solve(mechRhsVec, mechSolVec);

#ifdef USES_SILO
  siloWriter->registerVector(mechSolVec, meshAdapter, AMP::Mesh::Vertex, "Solution" );
  siloWriter->writeFile(exeName, 0);
#endif

  bvpOperator->apply(mechRhsVec, mechSolVec, mechResVec);

  double finalResidualNorm = mechResVec->L2Norm();

  std::cout<<"Final Residual Norm: "<<finalResidualNorm<<std::endl;

  if(finalResidualNorm > (1e-10*initResidualNorm)) {
    ut->failure("TrilinosMLSolver successfully solves a linear elasticity problem");
  } else {
    ut->passes("TrilinosMLSolver successfully solves a linear elasticity problem");
  }

  input_db.reset();

}

int main(int argc, char *argv[])
{
  AMP::AMPManager::startup(argc, argv);
  AMP::UnitTest ut;

  try {
    linearElasticTest(&ut);
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



