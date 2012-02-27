#include "utils/AMPManager.h"
#include "utils/UnitTest.h"
#include "utils/Utilities.h"

#include <iostream>
#include <string>
#include <cmath>

#include "boost/shared_ptr.hpp"

#include "operators/VolumeIntegralOperator.h"
#include "operators/NeutronicsRhs.h"

#include "utils/Database.h"
#include "utils/InputDatabase.h"
#include "utils/InputManager.h"
#include "utils/AMP_MPI.h"
#include "utils/AMPManager.h"
#include "utils/PIO.h"
#include "materials/Material.h"

#include "ampmesh/Mesh.h"
#include "ampmesh/SiloIO.h"

#include "discretization/DOF_Manager.h"
#include "discretization/simpleDOF_Manager.h"
#include "vectors/Variable.h"
#include "vectors/VectorBuilder.h"
#include "vectors/Variable.h"
#include "vectors/Vector.h"
#include "vectors/VectorSelector.h"

#include "operators/boundary/DirichletVectorCorrection.h"

#include "operators/diffusion/DiffusionLinearFEOperator.h"
#include "operators/diffusion/DiffusionNonlinearFEOperator.h"
#include "operators/diffusion/FickSoretNonlinearFEOperator.h"
#include "operators/BVPOperatorParameters.h"
#include "operators/LinearBVPOperator.h"
#include "operators/NonlinearBVPOperator.h"
#include "operators/ColumnOperator.h"
#include "operators/OperatorBuilder.h"

#include "../ColumnSolver.h"
#include "../PetscKrylovSolverParameters.h"
#include "../PetscKrylovSolver.h"
#include "../PetscSNESSolverParameters.h"
#include "../PetscSNESSolver.h"

#include "../TrilinosMLSolver.h"


void fickTest(AMP::UnitTest *ut, std::string exeName, std::vector<double> &results)
{
  std::string input_file = "input_" + exeName;
  std::string log_file = "output_" + exeName;

  AMP::PIO::logOnlyNodeZero(log_file);
  AMP::AMP_MPI globalComm(AMP_COMM_WORLD);

  boost::shared_ptr<AMP::InputDatabase> input_db(new AMP::InputDatabase("input_db"));
  AMP::InputManager::getManager()->parseInputFile(input_file, input_db);
  input_db->printClassData(AMP::plog);

  // Get the Mesh database and create the mesh parameters
  boost::shared_ptr<AMP::Database> database = input_db->getDatabase( "Mesh" );
  boost::shared_ptr<AMP::Mesh::MeshParameters> params(new AMP::Mesh::MeshParameters(database));
  params->setComm(globalComm);

  // Create the meshes from the input database
  AMP::Mesh::Mesh::shared_ptr manager = AMP::Mesh::Mesh::buildMesh(params);
  AMP::Mesh::Mesh::shared_ptr meshAdapter = manager->Subset( "cylinder" );

  //----------------------------------------------------------------------------------------------------------------------------------------------//
  // create a nonlinear BVP operator for nonlinear fick diffusion
  AMP_INSIST( input_db->keyExists("testNonlinearFickOperator"), "key missing!" );

  boost::shared_ptr<AMP::Operator::ElementPhysicsModel> fickTransportModel;
  boost::shared_ptr<AMP::Operator::NonlinearBVPOperator> nonlinearFickOperator = boost::dynamic_pointer_cast<
  AMP::Operator::NonlinearBVPOperator>(AMP::Operator::OperatorBuilder::createOperator(meshAdapter,
										      "testNonlinearFickOperator",
										      input_db,
										      fickTransportModel));

  //----------------------------------------------------------------------------------------------------------------------------------------------//
  // initialize the input variable
  boost::shared_ptr<AMP::Operator::DiffusionNonlinearFEOperator> fickVolumeOperator =
	        boost::dynamic_pointer_cast<AMP::Operator::DiffusionNonlinearFEOperator>(nonlinearFickOperator->getVolumeOperator());

  boost::shared_ptr<AMP::LinearAlgebra::Variable> fickVariable = fickVolumeOperator->getOutputVariable();

  // create solution, rhs, and residual vectors
  AMP::Discretization::DOFManager::shared_ptr nodalScalarDOF = AMP::Discretization::simpleDOFManager::create(meshAdapter,AMP::Mesh::Vertex,1,1,true);
  AMP::LinearAlgebra::Vector::shared_ptr solVec = AMP::LinearAlgebra::createVector( nodalScalarDOF, fickVariable, true );
  AMP::LinearAlgebra::Vector::shared_ptr rhsVec = AMP::LinearAlgebra::createVector( nodalScalarDOF, fickVariable, true );
  AMP::LinearAlgebra::Vector::shared_ptr resVec = AMP::LinearAlgebra::createVector( nodalScalarDOF, fickVariable, true );

  //----------------------------------------------------------------------------------------------------------------------------------------------//
  // register some variables for plotting
  AMP::Mesh::SiloIO::shared_ptr  siloWriter( new AMP::Mesh::SiloIO);
  siloWriter->registerVector( solVec, meshAdapter, AMP::Mesh::Vertex, "Solution" );
  siloWriter->registerVector( resVec, meshAdapter, AMP::Mesh::Vertex, "Residual" );

  //----------------------------------------------------------------------------------------------------------------------------------------------//
  // now construct the linear BVP operator for fick
  AMP_INSIST( input_db->keyExists("testLinearFickOperator"), "key missing!" );
  boost::shared_ptr<AMP::Operator::LinearBVPOperator> linearFickOperator = boost::dynamic_pointer_cast<AMP::Operator::LinearBVPOperator>(
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

  boost::shared_ptr<AMP::Database> nonlinearSolver_db = input_db->getDatabase("NonlinearSolver");
  boost::shared_ptr<AMP::Database> linearSolver_db = nonlinearSolver_db->getDatabase("LinearSolver");

  //----------------------------------------------------------------------------------------------------------------------------------------------//
  // initialize the nonlinear solver
  boost::shared_ptr<AMP::Solver::PetscSNESSolverParameters> nonlinearSolverParams(new
      AMP::Solver::PetscSNESSolverParameters(nonlinearSolver_db));

  // change the next line to get the correct communicator out
  nonlinearSolverParams->d_comm = globalComm;
  nonlinearSolverParams->d_pOperator = nonlinearFickOperator;
  nonlinearSolverParams->d_pInitialGuess = solVec;

  boost::shared_ptr<AMP::Solver::PetscSNESSolver> nonlinearSolver(new AMP::Solver::PetscSNESSolver(nonlinearSolverParams));

  //----------------------------------------------------------------------------------------------------------------------------------------------//
  boost::shared_ptr<AMP::Database> fickPreconditioner_db = linearSolver_db->getDatabase("Preconditioner");
  boost::shared_ptr<AMP::Solver::SolverStrategyParameters> fickPreconditionerParams(new AMP::Solver::SolverStrategyParameters(fickPreconditioner_db));
  fickPreconditionerParams->d_pOperator = linearFickOperator;
  boost::shared_ptr<AMP::Solver::TrilinosMLSolver> linearFickPreconditioner(new AMP::Solver::TrilinosMLSolver(fickPreconditionerParams));

  //----------------------------------------------------------------------------------------------------------------------------------------------//
  // register the preconditioner with the Jacobian free Krylov solver
  boost::shared_ptr<AMP::Solver::PetscKrylovSolver> linearSolver = nonlinearSolver->getKrylovSolver();

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

#ifdef USE_SILO
  siloWriter->writeFile( exeName , 0 );
#endif

  // store result
  {
      AMP::Mesh::MeshIterator iterator = meshAdapter->getIterator(AMP::Mesh::Vertex,0);
      size_t numNodes = iterator.size();
      results.resize(numNodes);
      std::vector<size_t> dofs;
      for(size_t iNode=0; iNode<numNodes; iNode++ ) {
          nodalScalarDOF->getDOFs(iterator->globalID(),dofs);
          size_t gid = dofs[0];
		  results[iNode] = solVec->getValueByGlobalID(gid);
          ++iterator;
      }
  }

  ut->passes(exeName);

}



void fickSoretTest(AMP::UnitTest *ut, std::string exeName, std::vector<double> &results)
{
  std::string input_file = "input_" + exeName;
  std::string log_file = "output_" + exeName;

  AMP::PIO::logOnlyNodeZero(log_file);
  AMP::AMP_MPI globalComm(AMP_COMM_WORLD);

  boost::shared_ptr<AMP::InputDatabase> input_db(new AMP::InputDatabase("input_db"));
  AMP::InputManager::getManager()->parseInputFile(input_file, input_db);
  input_db->printClassData(AMP::plog);

  // Get the Mesh database and create the mesh parameters
  boost::shared_ptr<AMP::Database> database = input_db->getDatabase( "Mesh" );
  boost::shared_ptr<AMP::Mesh::MeshParameters> params(new AMP::Mesh::MeshParameters(database));
  params->setComm(globalComm);

  // Create the meshes from the input database
  AMP::Mesh::Mesh::shared_ptr manager = AMP::Mesh::Mesh::buildMesh(params);
  AMP::Mesh::Mesh::shared_ptr meshAdapter = manager->Subset( "cylinder" );

  //----------------------------------------------------------------------------------------------------------------------------------------------//
  // create a nonlinear BVP operator for nonlinear Fick-Soret diffusion
  AMP_INSIST( input_db->keyExists("testNonlinearFickSoretBVPOperator"), "key missing!" );

  // Create nonlinear FickSoret BVP operator and access volume nonlinear FickSoret operator
  boost::shared_ptr<AMP::Operator::ElementPhysicsModel> elementPhysicsModel;
  boost::shared_ptr<AMP::Operator::Operator> nlinBVPOperator =
    AMP::Operator::OperatorBuilder::createOperator(meshAdapter,
						   "testNonlinearFickSoretBVPOperator",
						   input_db,
						   elementPhysicsModel);
  boost::shared_ptr<AMP::Operator::NonlinearBVPOperator> nlinBVPOp =
	        boost::dynamic_pointer_cast<AMP::Operator::NonlinearBVPOperator>(nlinBVPOperator);
  boost::shared_ptr<AMP::Operator::FickSoretNonlinearFEOperator> nlinOp =
		 boost::dynamic_pointer_cast<AMP::Operator::FickSoretNonlinearFEOperator>(nlinBVPOp->getVolumeOperator());
  boost::shared_ptr<AMP::Operator::DiffusionNonlinearFEOperator> fickOp =
		 boost::dynamic_pointer_cast<AMP::Operator::DiffusionNonlinearFEOperator>(nlinOp->getFickOperator());
  boost::shared_ptr<AMP::Operator::DiffusionNonlinearFEOperator> soretOp =
		 boost::dynamic_pointer_cast<AMP::Operator::DiffusionNonlinearFEOperator>(nlinOp->getSoretOperator());

  //----------------------------------------------------------------------------------------------------------------------------------------------//
  // use the linear BVP operator to create a Fick linear operator with bc's
  AMP_INSIST( input_db->keyExists("testLinearFickBVPOperator"), "key missing!" );

  boost::shared_ptr<AMP::Operator::Operator> linBVPOperator =
    AMP::Operator::OperatorBuilder::createOperator(meshAdapter,
						   "testLinearFickBVPOperator",
						   input_db,
						   elementPhysicsModel);
  boost::shared_ptr<AMP::Operator::LinearBVPOperator> linBVPOp =
	        boost::dynamic_pointer_cast<AMP::Operator::LinearBVPOperator>(linBVPOperator);
  //boost::shared_ptr<AMP::Operator::DiffusionLinearFEOperator> linOp =
 //		 boost::dynamic_pointer_cast<AMP::Operator::DiffusionLinearFEOperator>(linBVPOp->getVolumeOperator());

  //----------------------------------------------------------------------------------------------------------------------------------------------//
  // Set up input and output variables
  //AMP::LinearAlgebra::Variable::shared_ptr tVar(fickOp->getInputVariable(AMP::Operator::Diffusion::TEMPERATURE));
  //AMP::LinearAlgebra::Variable::shared_ptr cVar(fickOp->getInputVariable(AMP::Operator::Diffusion::CONCENTRATION));
  AMP::LinearAlgebra::Variable::shared_ptr tVar(fickOp->getInputVariable());
  AMP::LinearAlgebra::Variable::shared_ptr cVar(fickOp->getInputVariable());
  boost::shared_ptr<AMP::LinearAlgebra::Variable> fsOutVar(nlinBVPOp->getOutputVariable());

  //----------------------------------------------------------------------------------------------------------------------------------------------//
  // create solution, rhs, and residual vectors
  AMP::Discretization::DOFManager::shared_ptr nodalScalarDOF = AMP::Discretization::simpleDOFManager::create(meshAdapter,AMP::Mesh::Vertex,1,1,true);
  AMP::LinearAlgebra::Vector::shared_ptr solVec = AMP::LinearAlgebra::createVector( nodalScalarDOF, cVar, true );
  AMP::LinearAlgebra::Vector::shared_ptr rhsVec = AMP::LinearAlgebra::createVector( nodalScalarDOF, fsOutVar, true );
  AMP::LinearAlgebra::Vector::shared_ptr resVec = AMP::LinearAlgebra::createVector( nodalScalarDOF, fsOutVar, true );

  //----------------------------------------------------------------------------------------------------------------------------------------------//
  // create parameters for reset test and reset fick and soret operators
  AMP::LinearAlgebra::Vector::shared_ptr tVec = AMP::LinearAlgebra::createVector( nodalScalarDOF, tVar, true );
  tVec->setToScalar(300.);

  fickOp->setVector(0, tVec);
  soretOp->setVector(0, tVec);

  //----------------------------------------------------------------------------------------------------------------------------------------------//
  // register some variables for plotting
  AMP::Mesh::SiloIO::shared_ptr  siloWriter( new AMP::Mesh::SiloIO);
  siloWriter->registerVector( solVec, meshAdapter, AMP::Mesh::Vertex, "Solution" );
  siloWriter->registerVector( resVec, meshAdapter, AMP::Mesh::Vertex, "Residual" );


  //----------------------------------------------------------------------------------------------------------------------------------------------//
  //Initial guess

  solVec->setToScalar(.05);
  double initialGuessNorm  = solVec->L2Norm();
  std::cout << "initial guess norm = " << initialGuessNorm <<"\n";

  nlinBVPOp->modifyInitialSolutionVector(solVec);

  initialGuessNorm  = solVec->L2Norm();
  std::cout << "initial guess norm  after apply = " << initialGuessNorm <<"\n";

  rhsVec->setToScalar(0.0);
  nlinBVPOp->modifyRHSvector(rhsVec);

  //----------------------------------------------------------------------------------------------------------------------------------------------/

  boost::shared_ptr<AMP::Database> nonlinearSolver_db = input_db->getDatabase("NonlinearSolver"); 
  boost::shared_ptr<AMP::Database> linearSolver_db = nonlinearSolver_db->getDatabase("LinearSolver"); 

  //----------------------------------------------------------------------------------------------------------------------------------------------//
  // initialize the nonlinear solver
  boost::shared_ptr<AMP::Solver::PetscSNESSolverParameters> nonlinearSolverParams(new
      AMP::Solver::PetscSNESSolverParameters(nonlinearSolver_db));

  // change the next line to get the correct communicator out
  nonlinearSolverParams->d_comm = globalComm;
  nonlinearSolverParams->d_pOperator = nlinBVPOp;
  nonlinearSolverParams->d_pInitialGuess = solVec;

  boost::shared_ptr<AMP::Solver::PetscSNESSolver> nonlinearSolver(new AMP::Solver::PetscSNESSolver(nonlinearSolverParams));

  //----------------------------------------------------------------------------------------------------------------------------------------------//
  boost::shared_ptr<AMP::Database> fickPreconditioner_db = linearSolver_db->getDatabase("Preconditioner");
  boost::shared_ptr<AMP::Solver::SolverStrategyParameters> fickPreconditionerParams(new AMP::Solver::SolverStrategyParameters(fickPreconditioner_db));
  fickPreconditionerParams->d_pOperator = linBVPOp;
  boost::shared_ptr<AMP::Solver::TrilinosMLSolver> linearFickPreconditioner(new AMP::Solver::TrilinosMLSolver(fickPreconditionerParams));

  //----------------------------------------------------------------------------------------------------------------------------------------------//
  // register the preconditioner with the Jacobian free Krylov solver
  boost::shared_ptr<AMP::Solver::PetscKrylovSolver> linearSolver = nonlinearSolver->getKrylovSolver();

  linearSolver->setPreconditioner(linearFickPreconditioner);

  nlinBVPOp->apply(rhsVec, solVec, resVec, 1.0, -1.0);
  double initialResidualNorm  = resVec->L2Norm();

  AMP::pout<<"Initial Residual Norm: "<<initialResidualNorm<<std::endl;

  nonlinearSolver->setZeroInitialGuess(false);

  nonlinearSolver->solve(rhsVec, solVec);

  nlinBVPOp->apply(rhsVec, solVec, resVec, 1.0, -1.0);

  double finalResidualNorm  = resVec->L2Norm();

  std::cout<<"Final Residual Norm: "<<finalResidualNorm<<std::endl;

  solVec->makeConsistent ( AMP::LinearAlgebra::Vector::CONSISTENT_SET );
  resVec->makeConsistent ( AMP::LinearAlgebra::Vector::CONSISTENT_SET );

#ifdef USE_SILO
  siloWriter->writeFile( exeName , 0 );
#endif

  // store result
  {
      AMP::Mesh::MeshIterator iterator = meshAdapter->getIterator(AMP::Mesh::Vertex,0);
      size_t numNodes = iterator.size();
      results.resize(numNodes);
      std::vector<size_t> dofs;
      for(size_t iNode=0; iNode<numNodes; iNode++ ) {
          nodalScalarDOF->getDOFs(iterator->globalID(),dofs);
          size_t gid = dofs[0];
		  results[iNode] = solVec->getValueByGlobalID(gid);
          ++iterator;
      }
  }

  ut->passes(exeName);
}



int main(int argc, char *argv[])
{
    AMP::AMPManager::startup(argc, argv);
    AMP::UnitTest ut;

    try {
        std::vector<double>  fickOnly, fickSoretOff, fickSoretZero, fickOnlyReal, fickSoretOffReal;

        fickTest(&ut, "testPetscSNESSolver-NonlinearFick-cylinder-TUI-1", fickOnly);
        fickSoretTest(&ut, "testPetscSNESSolver-NonlinearFickSoret-cylinder-TUI-1", fickSoretOff);
        fickSoretTest(&ut, "testPetscSNESSolver-NonlinearFickSoret-cylinder-TUI-2", fickSoretZero);
        fickTest(&ut, "testPetscSNESSolver-NonlinearFick-cylinder-TUI-2", fickOnlyReal);
        fickSoretTest(&ut, "testPetscSNESSolver-NonlinearFickSoret-cylinder-TUI-3", fickSoretOffReal);
        AMP_INSIST(fickOnly.size()==fickSoretOff.size() and fickSoretOff.size()==fickSoretZero.size()
		        and fickOnlyReal.size()==fickSoretOffReal.size(),
		        "sizes of results do not match");

        double l2err1 = 0., l2err2 = 0.;
        for (size_t i=0; i<fickOnly.size(); i++) {
	        double err = fickOnly[i] - fickSoretOff[i];
	        l2err1 += err*err;
	        err = fickSoretOff[i] - fickSoretZero[i];
	        l2err2 += err*err;
        }
        l2err1 = sqrt(l2err1); l2err2 = sqrt(l2err2);

        std::cout << "fick/soretOff err = " << l2err1 << "  soretOff/soretZero err = " << l2err2 << std::endl;

        double l2err3 = 0.;
        for (size_t i=0; i<fickOnlyReal.size(); i++) {
	        double err = fickOnlyReal[i] - fickSoretOffReal[i];
	        l2err3 += err*err;
        }
        l2err3 = sqrt(l2err3);

        std::cout << "fick/soretOff real err = " << l2err3  << std::endl;

        if (l2err1 < 1.e-6 and l2err2 < 1.e-6 and l2err3 < 1.e-6) {
	        ut.passes("fick, fick-soret/off, and fick-soret/zero all agree");
        }
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


