#include "utils/AMPManager.h"
#include "utils/UnitTest.h"
#include "utils/Utilities.h"
#include <iostream>
#include <string>

#include "utils/shared_ptr.h"

#include "utils/Database.h"
#include "utils/InputDatabase.h"
#include "utils/InputManager.h"
#include "utils/AMP_MPI.h"
#include "utils/AMPManager.h"
#include "utils/PIO.h"
#include "materials/Material.h"


#include "ampmesh/MeshVariable.h"
#include "utils/Writer.h"


#include "operators/diffusion/DiffusionLinearFEOperator.h"
#include "operators/diffusion/DiffusionNonlinearFEOperator.h"
#include "operators/DirichletVectorCorrection.h"
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


void myTest(AMP::UnitTest *ut, std::string exeName)
{
  std::string input_file = "input_" + exeName;
  std::string log_file = "output_" + exeName;

  AMP::PIO::logOnlyNodeZero(log_file);

  AMP::shared_ptr<AMP::InputDatabase> input_db(new AMP::InputDatabase("input_db"));
  AMP::AMP_MPI globalComm(AMP_COMM_WORLD);
  AMP::InputManager::getManager()->parseInputFile(input_file, input_db);
  input_db->printClassData(AMP::plog);

  AMP::Mesh::MeshManagerParameters::shared_ptr  meshmgrParams ( new AMP::Mesh::MeshManagerParameters ( input_db ) );
  AMP::Mesh::MeshManager::shared_ptr  manager ( new AMP::Mesh::MeshManager ( meshmgrParams ) );
  AMP::Mesh::MeshManager::Adapter::shared_ptr meshAdapter = manager->getMesh ( "cylinder" );

  //----------------------------------------------------------------------------------------------------------------------------------------------//
  // create a nonlinear BVP operator for nonlinear oxygen diffusion
  AMP_INSIST( input_db->keyExists("testNonlinearOxygenOperator"), "key missing!" );

  AMP::shared_ptr<AMP::Operator::ElementPhysicsModel> oxygenTransportModel;
  AMP::shared_ptr<AMP::Database> nonlinearOxygenDatabase = input_db->getDatabase("testNonlinearOxygenOperator");
  AMP::shared_ptr<AMP::Operator::NonlinearBVPOperator> nonlinearOxygenOperator = AMP::dynamic_pointer_cast<
    AMP::Operator::NonlinearBVPOperator>(AMP::Operator::OperatorBuilder::createOperator(meshAdapter, nonlinearOxygenDatabase, oxygenTransportModel));

  //----------------------------------------------------------------------------------------------------------------------------------------------//
  // create a nonlinear BVP operator for nonlinear thermal diffusion
  AMP_INSIST( input_db->keyExists("testNonlinearThermalOperator"), "key missing!" );

  AMP::shared_ptr<AMP::Operator::ElementPhysicsModel> thermalTransportModel;
  AMP::shared_ptr<AMP::Database> nonlinearThermalDatabase = input_db->getDatabase("testNonlinearThermalOperator");
  AMP::shared_ptr<AMP::Operator::NonlinearBVPOperator> nonlinearThermalOperator = AMP::dynamic_pointer_cast<
    AMP::Operator::NonlinearBVPOperator>(AMP::Operator::OperatorBuilder::createOperator(meshAdapter, nonlinearThermalDatabase, thermalTransportModel));

  //----------------------------------------------------------------------------------------------------------------------------------------------//
  // create a column operator object for nonlinear thermal-oxygen diffusion
  AMP::shared_ptr<AMP::Operator::OperatorParameters> params;
  AMP::shared_ptr<AMP::Operator::ColumnOperator> nonlinearThermalOxygenOperator(new AMP::Operator::ColumnOperator(params));
  nonlinearThermalOxygenOperator->append(nonlinearOxygenOperator);
  nonlinearThermalOxygenOperator->append(nonlinearThermalOperator);

  //----------------------------------------------------------------------------------------------------------------------------------------------//
  AMP::shared_ptr<AMP::Operator::DiffusionNonlinearFEOperator> oxygenVolumeOperator = AMP::dynamic_pointer_cast<
    AMP::Operator::DiffusionNonlinearFEOperator>(nonlinearOxygenOperator->getVolumeOperator());
  AMP::shared_ptr<AMP::Operator::DiffusionNonlinearFEOperator> thermalVolumeOperator = AMP::dynamic_pointer_cast<
    AMP::Operator::DiffusionNonlinearFEOperator>(nonlinearThermalOperator->getVolumeOperator());

  // initialize the output multi-variable
  AMP::LinearAlgebra::Variable::shared_ptr columnVariable = nonlinearThermalOxygenOperator->getOutputVariable();

  // create solution, rhs, and residual vectors
  AMP::LinearAlgebra::Vector::shared_ptr solVec = meshAdapter->createVector( columnVariable );
  AMP::LinearAlgebra::Vector::shared_ptr rhsVec = meshAdapter->createVector( columnVariable );
  AMP::LinearAlgebra::Vector::shared_ptr resVec = meshAdapter->createVector( columnVariable );

  // just making sure
  rhsVec->zero();

  // create the following shared pointers for ease of use
  AMP::LinearAlgebra::Vector::shared_ptr nullVec;

  AMP::LinearAlgebra::Variable::shared_ptr oxygenVariable = oxygenVolumeOperator->getOutputVariable();
  AMP::LinearAlgebra::Variable::shared_ptr thermalVariable = thermalVolumeOperator->getOutputVariable(); 

  AMP::LinearAlgebra::Vector::shared_ptr oxygenNlSolVec = solVec->subsetVectorForVariable( oxygenVariable );
  AMP::LinearAlgebra::Vector::shared_ptr oxygenNlRhsVec = rhsVec->subsetVectorForVariable( oxygenVariable );
  AMP::LinearAlgebra::Vector::shared_ptr oxygenNlResVec = resVec->subsetVectorForVariable( oxygenVariable );

  AMP::LinearAlgebra::Vector::shared_ptr thermalNlSolVec = solVec->subsetVectorForVariable( thermalVariable );
  AMP::LinearAlgebra::Vector::shared_ptr thermalNlRhsVec = rhsVec->subsetVectorForVariable( thermalVariable );
  AMP::LinearAlgebra::Vector::shared_ptr thermalNlResVec = resVec->subsetVectorForVariable( thermalVariable );

  //----------------------------------------------------------------------------------------------------------------------------------------------//
  // register some variables for plotting
  meshAdapter->registerVectorAsData ( oxygenNlSolVec , "OxygenSolution");
  meshAdapter->registerVectorAsData ( thermalNlSolVec, "ThermalSolution" );
  meshAdapter->registerVectorAsData ( oxygenNlResVec, "OxygenResidual" );
  meshAdapter->registerVectorAsData ( thermalNlResVec, "ThermalResidual" );

  //----------------------------------------------------------------------------------------------------------------------------------------------//
  // now construct the linear BVP operator for oxygen
  AMP_INSIST( input_db->keyExists("testLinearOxygenOperator"), "key missing!" );
  AMP::shared_ptr<AMP::Database> linearOxygenDatabase = input_db->getDatabase("testLinearOxygenOperator");
  AMP::shared_ptr<AMP::Operator::LinearBVPOperator> linearOxygenOperator = AMP::dynamic_pointer_cast<AMP::Operator::LinearBVPOperator>(
      AMP::Operator::OperatorBuilder::createOperator(meshAdapter, linearOxygenDatabase, oxygenTransportModel));

  //----------------------------------------------------------------------------------------------------------------------------------------------//
  // now construct the linear BVP operator for thermal
  AMP_INSIST( input_db->keyExists("testLinearThermalOperator"), "key missing!" );
  AMP::shared_ptr<AMP::Database> linearThermalDatabase = input_db->getDatabase("testLinearThermalOperator");
  AMP::shared_ptr<AMP::Operator::LinearBVPOperator> linearThermalOperator = AMP::dynamic_pointer_cast<AMP::Operator::LinearBVPOperator>(
      AMP::Operator::OperatorBuilder::createOperator(meshAdapter, linearThermalDatabase, thermalTransportModel));

  //----------------------------------------------------------------------------------------------------------------------------------------------//
  // create a column operator object for linear thermal-oxygen
  AMP::shared_ptr<AMP::Operator::ColumnOperator> linearThermalOxygenOperator(new AMP::Operator::ColumnOperator(params));
  linearThermalOxygenOperator->append(linearOxygenOperator);
  linearThermalOxygenOperator->append(linearThermalOperator);

  //----------------------------------------------------------------------------------------------------------------------------------------------//
  //Random initial guess
  solVec->setRandomValues();
  const double referenceTemperature=750.0;
  thermalNlSolVec->addScalar(thermalNlSolVec, referenceTemperature);
  oxygenNlSolVec->scale(0.01, oxygenNlSolVec);

  nonlinearOxygenOperator->modifyInitialSolutionVector(solVec);
  nonlinearThermalOperator->modifyInitialSolutionVector(solVec);

  //----------------------------------------------------------------------------------------------------------------------------------------------//

  //We need to reset the linear operator before the solve since TrilinosML does
  //the factorization of the matrix during construction and so the matrix must
  //be correct before constructing the TrilinosML object.
  //The thermal operator does not expect an apply to be called before calling
  //getJacobianParams and so it need not be called. So, any of the following
  //apply calls will work:
  //oxygenVolumeOperator->apply(nullVec, solVec, resVec, 1.0, 0.0);
  //nonlinearOxygenOperator->apply(nullVec, solVec, resVec, 1.0, 0.0);
  nonlinearThermalOxygenOperator->apply(solVec, resVec);
  linearThermalOxygenOperator->reset(nonlinearThermalOxygenOperator->getParameters("Jacobian", solVec));
  //----------------------------------------------------------------------------------------------------------------------------------------------/

  //----------------------------------------------------------------------------------------------------------------------------------------------//
  AMP::shared_ptr<AMP::Database> nonlinearSolver_db = input_db->getDatabase("NonlinearSolver"); 
  AMP::shared_ptr<AMP::Database> linearSolver_db = nonlinearSolver_db->getDatabase("LinearSolver"); 

  //----------------------------------------------------------------------------------------------------------------------------------------------//
  // initialize the nonlinear solver
  AMP::shared_ptr<AMP::Solver::PetscSNESSolverParameters> nonlinearSolverParams(new
      AMP::Solver::PetscSNESSolverParameters(nonlinearSolver_db));

  // change the next line to get the correct communicator out
  nonlinearSolverParams->d_comm = globalComm;
  nonlinearSolverParams->d_pOperator = nonlinearThermalOxygenOperator;
  nonlinearSolverParams->d_pInitialGuess = solVec;

  AMP::shared_ptr<AMP::Solver::PetscSNESSolver> nonlinearSolver(new AMP::Solver::PetscSNESSolver(nonlinearSolverParams));

  //----------------------------------------------------------------------------------------------------------------------------------------------//
  // initialize the column preconditioner which is a diagonal block preconditioner
  AMP::shared_ptr<AMP::Database> columnPreconditioner_db = linearSolver_db->getDatabase("Preconditioner");
  AMP::shared_ptr<AMP::Solver::SolverStrategyParameters> columnPreconditionerParams(new AMP::Solver::SolverStrategyParameters(columnPreconditioner_db));
  columnPreconditionerParams->d_pOperator = linearThermalOxygenOperator;
  AMP::shared_ptr<AMP::Solver::ColumnSolver> columnPreconditioner(new AMP::Solver::ColumnSolver(columnPreconditionerParams));

  AMP::shared_ptr<AMP::Database> oxygenPreconditioner_db = columnPreconditioner_db->getDatabase("oxygenPreconditioner"); 
  AMP::shared_ptr<AMP::Solver::SolverStrategyParameters> oxygenPreconditionerParams(new AMP::Solver::SolverStrategyParameters(oxygenPreconditioner_db));
  oxygenPreconditionerParams->d_pOperator = linearOxygenOperator;
  AMP::shared_ptr<AMP::Solver::TrilinosMLSolver> linearOxygenPreconditioner(new AMP::Solver::TrilinosMLSolver(oxygenPreconditionerParams));

  AMP::shared_ptr<AMP::Database> thermalPreconditioner_db = columnPreconditioner_db->getDatabase("thermalPreconditioner"); 
  AMP::shared_ptr<AMP::Solver::SolverStrategyParameters> thermalPreconditionerParams(new AMP::Solver::SolverStrategyParameters(thermalPreconditioner_db));
  thermalPreconditionerParams->d_pOperator = linearThermalOperator;
  AMP::shared_ptr<AMP::Solver::TrilinosMLSolver> linearThermalPreconditioner(new AMP::Solver::TrilinosMLSolver(thermalPreconditionerParams));

  columnPreconditioner->append(linearOxygenPreconditioner);
  columnPreconditioner->append(linearThermalPreconditioner);

  //----------------------------------------------------------------------------------------------------------------------------------------------//
  // register the preconditioner with the Jacobian free Krylov solver
  AMP::shared_ptr<AMP::Solver::PetscKrylovSolver> linearSolver = nonlinearSolver->getKrylovSolver();

  linearSolver->setPreconditioner(columnPreconditioner);

  nonlinearThermalOxygenOperator->residual(rhsVec, solVec, resVec);
  double initialResidualNorm  = resVec->L2Norm();

  AMP::pout<<"Initial Residual Norm: "<<initialResidualNorm<<std::endl;

  nonlinearSolver->setZeroInitialGuess(false);

  nonlinearSolver->solve(rhsVec, solVec);

  nonlinearThermalOxygenOperator->residual(rhsVec, solVec, resVec);

  double finalResidualNorm  = resVec->L2Norm();

  AMP::pout<<"Final Residual Norm: "<<finalResidualNorm<<std::endl;

#ifdef USE_EXT_SILO
  manager->writeFile<AMP::Mesh::SiloIO> ( exeName , 1 );
#endif

  if(finalResidualNorm>initialResidualNorm*1.0e-10+1.0e-05)
  {
    ITFAILS;
  }
  else
  {
    ut.passes("PetscSNES Solver successfully solves a nonlinear thermal-oxygen diffusion equation with JFNK, FGMRES for Krylov, block diagonal preconditioning with ML solvers");
  }
  ut.passes(exeName);

}

int main(int argc, char *argv[])
{
    AMP::AMPManager::startup(argc, argv);
    AMP::UnitTest ut;

  AMP::AMP_MPI::initialize();  std::vector<std::string> exeNames;
//  exeNames.push_back("testPetscSNESSolver-NonlinearThermalOxygenDiffusion-1");
  exeNames.push_back("testPetscSNESSolver-NonlinearThermalOxygenDiffusion-2");

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


