#include "utils/AMPManager.h"
#include "utils/UnitTest.h"
#include "utils/Utilities.h"
#include <iostream>
#include <string>

#include "utils/shared_ptr.h"

#include "operators/libmesh/VolumeIntegralOperator.h"
#include "operators/NeutronicsRhs.h"

#include "utils/Database.h"
#include "utils/InputDatabase.h"
#include "utils/InputManager.h"
#include "utils/AMP_MPI.h"
#include "utils/AMPManager.h"
#include "utils/PIO.h"
#include "materials/Material.h"

#include "ampmesh/Mesh.h"
#include "utils/Writer.h"

#include "vectors/Variable.h"
#include "vectors/MultiVariable.h"
#include "vectors/Vector.h"
#include "vectors/VectorBuilder.h"

#include "discretization/DOF_Manager.h"
#include "discretization/simpleDOF_Manager.h"

#include "operators/mechanics/MechanicsLinearFEOperator.h"
#include "operators/mechanics/MechanicsNonlinearFEOperator.h"

#include "operators/diffusion/DiffusionLinearFEOperator.h"
#include "operators/diffusion/DiffusionNonlinearFEOperator.h"

#include "operators/boundary/DirichletVectorCorrection.h"

#include "operators/BVPOperatorParameters.h"
#include "operators/LinearBVPOperator.h"
#include "operators/NonlinearBVPOperator.h"
#include "operators/ColumnOperator.h"
#include "operators/OperatorBuilder.h"

#include "solvers/ColumnSolver.h"
#include "solvers/petsc/PetscKrylovSolverParameters.h"
#include "solvers/petsc/PetscKrylovSolver.h"
#include "solvers/petsc/PetscSNESSolverParameters.h"
#include "solvers/petsc/PetscSNESSolver.h"

#include "solvers/trilinos/TrilinosMLSolver.h"


void myTest(AMP::UnitTest *ut, std::string exeName)
{
  std::string input_file = "input_" + exeName;
  std::string log_file = "output_" + exeName;

  AMP::PIO::logOnlyNodeZero(log_file);
  AMP::AMP_MPI globalComm(AMP_COMM_WORLD);

  // Read the input file
  AMP::shared_ptr<AMP::InputDatabase>  input_db ( new AMP::InputDatabase ( "input_db" ) );
  AMP::InputManager::getManager()->parseInputFile ( input_file , input_db );
  input_db->printClassData (AMP::plog);

  // Get the Mesh database and create the mesh parameters
  AMP::shared_ptr<AMP::Database> database = input_db->getDatabase( "Mesh" );
  AMP::shared_ptr<AMP::Mesh::MeshParameters> params(new AMP::Mesh::MeshParameters(database));
  params->setComm(globalComm);

  // Create the meshes from the input database
  AMP::Mesh::Mesh::shared_ptr meshAdapter = AMP::Mesh::Mesh::buildMesh(params);

  AMP::Discretization::DOFManager::shared_ptr DOF_scalar = AMP::Discretization::simpleDOFManager::create(meshAdapter,AMP::Mesh::Vertex,1,1,true);


  //----------------------------------------------------------------------------------------------------------------------------------------------//
  // create a nonlinear BVP operator for nonlinear fick diffusion
  AMP_INSIST( input_db->keyExists("testNonlinearFickOperator"), "key missing!" );

  AMP::shared_ptr<AMP::Operator::ElementPhysicsModel> fickTransportModel;
  AMP::shared_ptr<AMP::Operator::NonlinearBVPOperator> nonlinearFickOperator = AMP::dynamic_pointer_cast<
    AMP::Operator::NonlinearBVPOperator>(AMP::Operator::OperatorBuilder::createOperator(meshAdapter,
											"testNonlinearFickOperator",
											input_db,
											fickTransportModel));

  //----------------------------------------------------------------------------------------------------------------------------------------------//
  // initialize the input variable
  AMP::shared_ptr<AMP::Operator::DiffusionNonlinearFEOperator> fickVolumeOperator =
		  AMP::dynamic_pointer_cast<AMP::Operator::DiffusionNonlinearFEOperator>(nonlinearFickOperator->getVolumeOperator());

  AMP::shared_ptr<AMP::LinearAlgebra::Variable> fickVariable = fickVolumeOperator->getOutputVariable();

  // create solution, rhs, and residual vectors
  AMP::LinearAlgebra::Vector::shared_ptr solVec = AMP::LinearAlgebra::createVector( DOF_scalar, fickVariable );
  AMP::LinearAlgebra::Vector::shared_ptr rhsVec = AMP::LinearAlgebra::createVector( DOF_scalar, fickVariable );
  AMP::LinearAlgebra::Vector::shared_ptr resVec = AMP::LinearAlgebra::createVector( DOF_scalar, fickVariable );

  //----------------------------------------------------------------------------------------------------------------------------------------------//
  // register some variables for plotting
#ifdef USE_EXT_SILO
  AMP::Utilities::Writer::shared_ptr siloWriter = AMP::Utilities::Writer::buildWriter("Silo");
  siloWriter->registerVector( solVec, meshAdapter, AMP::Mesh::Vertex, "Solution" );
  siloWriter->registerVector( resVec, meshAdapter, AMP::Mesh::Vertex, "Residual" );
#endif

  //----------------------------------------------------------------------------------------------------------------------------------------------//
  // now construct the linear BVP operator for fick
  AMP_INSIST( input_db->keyExists("testLinearFickOperator"), "key missing!" );
  AMP::shared_ptr<AMP::Operator::LinearBVPOperator> linearFickOperator = AMP::dynamic_pointer_cast<AMP::Operator::LinearBVPOperator>(
																	 AMP::Operator::OperatorBuilder::createOperator(meshAdapter,
																							"testLinearFickOperator",
																							input_db,
																							fickTransportModel));

  //----------------------------------------------------------------------------------------------------------------------------------------------//
  //Initial guess

  solVec->setToScalar(.05);
  double initialGuessNorm  = solVec->L2Norm();
  std::cout << "initial guess norm = " << initialGuessNorm <<"\n";

  nonlinearFickOperator->modifyInitialSolutionVector(solVec);

  initialGuessNorm  = solVec->L2Norm();
  std::cout << "initial guess norm  after apply = " << initialGuessNorm <<"\n";

  rhsVec->setToScalar(0.0);
  nonlinearFickOperator->modifyRHSvector(rhsVec);

  //----------------------------------------------------------------------------------------------------------------------------------------------/

  AMP::shared_ptr<AMP::Database> nonlinearSolver_db = input_db->getDatabase("NonlinearSolver"); 
  AMP::shared_ptr<AMP::Database> linearSolver_db = nonlinearSolver_db->getDatabase("LinearSolver"); 

  //----------------------------------------------------------------------------------------------------------------------------------------------//
  // initialize the nonlinear solver
  AMP::shared_ptr<AMP::Solver::PetscSNESSolverParameters> nonlinearSolverParams(new
      AMP::Solver::PetscSNESSolverParameters(nonlinearSolver_db));

  // change the next line to get the correct communicator out
  nonlinearSolverParams->d_comm = globalComm;
  nonlinearSolverParams->d_pOperator = nonlinearFickOperator;
  nonlinearSolverParams->d_pInitialGuess = solVec;

  AMP::shared_ptr<AMP::Solver::PetscSNESSolver> nonlinearSolver(new AMP::Solver::PetscSNESSolver(nonlinearSolverParams));

  //----------------------------------------------------------------------------------------------------------------------------------------------//
  AMP::shared_ptr<AMP::Database> fickPreconditioner_db = linearSolver_db->getDatabase("Preconditioner");
  AMP::shared_ptr<AMP::Solver::SolverStrategyParameters> fickPreconditionerParams(new AMP::Solver::SolverStrategyParameters(fickPreconditioner_db));
  fickPreconditionerParams->d_pOperator = linearFickOperator;
  AMP::shared_ptr<AMP::Solver::TrilinosMLSolver> linearFickPreconditioner(new AMP::Solver::TrilinosMLSolver(fickPreconditionerParams));

  //----------------------------------------------------------------------------------------------------------------------------------------------//
  // register the preconditioner with the Jacobian free Krylov solver
  AMP::shared_ptr<AMP::Solver::PetscKrylovSolver> linearSolver = nonlinearSolver->getKrylovSolver();

  linearSolver->setPreconditioner(linearFickPreconditioner);

  nonlinearFickOperator->apply(rhsVec, solVec, resVec, 1.0, -1.0);
  double initialResidualNorm  = resVec->L2Norm();

  AMP::pout<<"Initial Residual Norm: "<<initialResidualNorm<<std::endl;

  nonlinearSolver->setZeroInitialGuess(false);

  nonlinearSolver->solve(rhsVec, solVec);

  nonlinearFickOperator->apply(rhsVec, solVec, resVec, 1.0, -1.0);

  double finalResidualNorm  = resVec->L2Norm();

  std::cout<<"Final Residual Norm: "<<finalResidualNorm<<std::endl;

  solVec->makeConsistent ( AMP::LinearAlgebra::Vector::CONSISTENT_SET );
  resVec->makeConsistent ( AMP::LinearAlgebra::Vector::CONSISTENT_SET );

#ifdef USE_EXT_SILO
  siloWriter->writeFile( exeName, 0 );
#endif

  if(finalResidualNorm>1.0e-08) {
    ut->failure(exeName);
  } else {
	ut->passes(exeName);
  }

}

int main(int argc, char *argv[])
{
    AMP::AMPManager::startup(argc, argv);
    AMP::UnitTest ut;

    std::vector<std::string> exeNames;
    exeNames.push_back("testPetscSNESSolver-NonlinearFick-cylinder-1a");
    exeNames.push_back("testPetscSNESSolver-NonlinearFick-cylinder-1b");
    exeNames.push_back("testPetscSNESSolver-NonlinearFick-cylinder-1c");
    exeNames.push_back("testPetscSNESSolver-NonlinearFick-cylinder-1d");

    for(unsigned int i = 0; i < exeNames.size(); i++) {
        try {      
            myTest(&ut, exeNames[i]);
        } catch (std::exception &err) {
            std::cout << "ERROR: While testing "<<argv[0] << err.what() << std::endl;
            ut.failure("ERROR: While testing");
        } catch( ... ) {
            std::cout << "ERROR: While testing "<<argv[0] << "An unknown exception was thrown." << std::endl;
            ut.failure("ERROR: While testing");
        }
    }
   
    ut.report();

    int num_failed = ut.NumFailGlobal();
    AMP::AMPManager::shutdown();
    return num_failed;
}   


