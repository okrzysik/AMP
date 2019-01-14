#include "AMP/ampmesh/Mesh.h"
#include "AMP/materials/Material.h"
#include "AMP/operators/BVPOperatorParameters.h"
#include "AMP/operators/ColumnOperator.h"
#include "AMP/operators/DirichletVectorCorrection.h"
#include "AMP/operators/LinearBVPOperator.h"
#include "AMP/operators/NonlinearBVPOperator.h"
#include "AMP/operators/OperatorBuilder.h"
#include "AMP/operators/diffusion/DiffusionLinearFEOperator.h"
#include "AMP/operators/diffusion/DiffusionNonlinearFEOperator.h"
#include "AMP/solvers/ColumnSolver.h"
#include "AMP/solvers/PetscKrylovSolver.h"
#include "AMP/solvers/PetscKrylovSolverParameters.h"
#include "AMP/solvers/PetscSNESSolver.h"
#include "AMP/solvers/PetscSNESSolverParameters.h"
#include "AMP/solvers/TrilinosMLSolver.h"
#include "AMP/utils/AMPManager.h"
#include "AMP/utils/AMP_MPI.h"
#include "AMP/utils/Database.h"
#include "AMP/utils/InputDatabase.h"
#include "AMP/utils/InputManager.h"
#include "AMP/utils/PIO.h"
#include "AMP/utils/UnitTest.h"
#include "AMP/utils/Utilities.h"
#include "AMP/utils/Writer.h"
#include "AMP/utils/shared_ptr.h"

#include <iostream>
#include <string>


void myTest( AMP::UnitTest *ut, const std::string &exeName )
{
    std::string input_file = "input_" + exeName;
    std::string log_file   = "output_" + exeName;

    AMP::PIO::logOnlyNodeZero( log_file );

    auto input_db = AMP::make_shared<AMP::InputDatabase>( "input_db" );
    AMP::AMP_MPI globalComm( AMP_COMM_WORLD );
    AMP::InputManager::getManager()->parseInputFile( input_file, input_db );
    input_db->printClassData( AMP::plog );

    auto meshmgrParams = AMP::make_shared<AMP::Mesh::MeshManagerParameters>( input_db );
    auto manager       = AMP::make_shared<AMP::Mesh::MeshManager>( meshmgrParams );
    auto meshAdapter   = manager->getMesh( "cylinder" );

    //----------------------------------------------------------------------------------------------------------------------------------------------//
    // create a nonlinear BVP operator for nonlinear oxygen diffusion
    AMP_INSIST( input_db->keyExists( "testNonlinearOxygenOperator" ), "key missing!" );

    AMP::shared_ptr<AMP::Operator::ElementPhysicsModel> oxygenTransportModel;
    auto nonlinearOxygenDatabase = input_db->getDatabase( "testNonlinearOxygenOperator" );
    auto nonlinearOxygenOperator = AMP::dynamic_pointer_cast<AMP::Operator::NonlinearBVPOperator>(
        AMP::Operator::OperatorBuilder::createOperator(
            meshAdapter, nonlinearOxygenDatabase, oxygenTransportModel ) );

    //----------------------------------------------------------------------------------------------------------------------------------------------//
    // create a nonlinear BVP operator for nonlinear thermal diffusion
    AMP_INSIST( input_db->keyExists( "testNonlinearThermalOperator" ), "key missing!" );

    AMP::shared_ptr<AMP::Operator::ElementPhysicsModel> thermalTransportModel;
    auto nonlinearThermalDatabase = input_db->getDatabase( "testNonlinearThermalOperator" );
    auto nonlinearThermalOperator = AMP::dynamic_pointer_cast<AMP::Operator::NonlinearBVPOperator>(
        AMP::Operator::OperatorBuilder::createOperator(
            meshAdapter, nonlinearThermalDatabase, thermalTransportModel ) );

    //----------------------------------------------------------------------------------------------------------------------------------------------//
    // create a column operator object for nonlinear thermal-oxygen diffusion
    AMP::shared_ptr<AMP::Operator::OperatorParameters> params;
    auto nonlinearThermalOxygenOperator = AMP::make_shared<AMP::Operator::ColumnOperator>( params );
    nonlinearThermalOxygenOperator->append( nonlinearOxygenOperator );
    nonlinearThermalOxygenOperator->append( nonlinearThermalOperator );

    //----------------------------------------------------------------------------------------------------------------------------------------------//
    auto oxygenVolumeOperator =
        AMP::dynamic_pointer_cast<AMP::Operator::DiffusionNonlinearFEOperator>(
            nonlinearOxygenOperator->getVolumeOperator() );
    auto thermalVolumeOperator =
        AMP::dynamic_pointer_cast<AMP::Operator::DiffusionNonlinearFEOperator>(
            nonlinearThermalOperator->getVolumeOperator() );

    // initialize the output multi-variable
    auto columnVariable = nonlinearThermalOxygenOperator->getOutputVariable();

    // create solution, rhs, and residual vectors
    auto solVec = meshAdapter->createVector( columnVariable );
    auto rhsVec = meshAdapter->createVector( columnVariable );
    auto resVec = meshAdapter->createVector( columnVariable );

    // just making sure
    rhsVec->zero();

    // create the following shared pointers for ease of use
    AMP::LinearAlgebra::Vector::shared_ptr nullVec;

    auto oxygenVariable  = oxygenVolumeOperator->getOutputVariable();
    auto thermalVariable = thermalVolumeOperator->getOutputVariable();

    auto oxygenNlSolVec = solVec->subsetVectorForVariable( oxygenVariable );
    auto oxygenNlRhsVec = rhsVec->subsetVectorForVariable( oxygenVariable );
    auto oxygenNlResVec = resVec->subsetVectorForVariable( oxygenVariable );

    auto thermalNlSolVec = solVec->subsetVectorForVariable( thermalVariable );
    auto thermalNlRhsVec = rhsVec->subsetVectorForVariable( thermalVariable );
    auto thermalNlResVec = resVec->subsetVectorForVariable( thermalVariable );

    //----------------------------------------------------------------------------------------------------------------------------------------------//
    // register some variables for plotting
    meshAdapter->registerVectorAsData( oxygenNlSolVec, "OxygenSolution" );
    meshAdapter->registerVectorAsData( thermalNlSolVec, "ThermalSolution" );
    meshAdapter->registerVectorAsData( oxygenNlResVec, "OxygenResidual" );
    meshAdapter->registerVectorAsData( thermalNlResVec, "ThermalResidual" );

    //----------------------------------------------------------------------------------------------------------------------------------------------//
    // now construct the linear BVP operator for oxygen
    AMP_INSIST( input_db->keyExists( "testLinearOxygenOperator" ), "key missing!" );
    auto linearOxygenDatabase = input_db->getDatabase( "testLinearOxygenOperator" );
    auto linearOxygenOperator = AMP::dynamic_pointer_cast<AMP::Operator::LinearBVPOperator>(
        AMP::Operator::OperatorBuilder::createOperator(
            meshAdapter, linearOxygenDatabase, oxygenTransportModel ) );

    //----------------------------------------------------------------------------------------------------------------------------------------------//
    // now construct the linear BVP operator for thermal
    AMP_INSIST( input_db->keyExists( "testLinearThermalOperator" ), "key missing!" );
    auto linearThermalDatabase = input_db->getDatabase( "testLinearThermalOperator" );
    auto linearThermalOperator = AMP::dynamic_pointer_cast<AMP::Operator::LinearBVPOperator>(
        AMP::Operator::OperatorBuilder::createOperator(
            meshAdapter, linearThermalDatabase, thermalTransportModel ) );

    //----------------------------------------------------------------------------------------------------------------------------------------------//
    // create a column operator object for linear thermal-oxygen
    auto linearThermalOxygenOperator = AMP::make_shared<AMP::Operator::ColumnOperator>( params );
    linearThermalOxygenOperator->append( linearOxygenOperator );
    linearThermalOxygenOperator->append( linearThermalOperator );

    //----------------------------------------------------------------------------------------------------------------------------------------------//
    // Random initial guess
    solVec->setRandomValues();
    const double referenceTemperature = 750.0;
    thermalNlSolVec->addScalar( thermalNlSolVec, referenceTemperature );
    oxygenNlSolVec->scale( 0.01, oxygenNlSolVec );

    nonlinearOxygenOperator->modifyInitialSolutionVector( solVec );
    nonlinearThermalOperator->modifyInitialSolutionVector( solVec );

    //----------------------------------------------------------------------------------------------------------------------------------------------//

    // We need to reset the linear operator before the solve since TrilinosML does
    // the factorization of the matrix during construction and so the matrix must
    // be correct before constructing the TrilinosML object.
    // The thermal operator does not expect an apply to be called before calling
    // getJacobianParams and so it need not be called. So, any of the following
    // apply calls will work:
    // oxygenVolumeOperator->apply(nullVec, solVec, resVec, 1.0, 0.0);
    // nonlinearOxygenOperator->apply(nullVec, solVec, resVec, 1.0, 0.0);
    nonlinearThermalOxygenOperator->apply( solVec, resVec );
    linearThermalOxygenOperator->reset(
        nonlinearThermalOxygenOperator->getParameters( "Jacobian", solVec ) );
    //----------------------------------------------------------------------------------------------------------------------------------------------/

    //----------------------------------------------------------------------------------------------------------------------------------------------//
    auto nonlinearSolver_db = input_db->getDatabase( "NonlinearSolver" );
    auto linearSolver_db    = nonlinearSolver_db->getDatabase( "LinearSolver" );

    //----------------------------------------------------------------------------------------------------------------------------------------------//
    // initialize the nonlinear solver
    auto nonlinearSolverParams =
        AMP::make_shared<AMP::Solver::PetscSNESSolverParameters>( nonlinearSolver_db );

    // change the next line to get the correct communicator out
    nonlinearSolverParams->d_comm          = globalComm;
    nonlinearSolverParams->d_pOperator     = nonlinearThermalOxygenOperator;
    nonlinearSolverParams->d_pInitialGuess = solVec;

    auto nonlinearSolver = AMP::make_shared<AMP::Solver::PetscSNESSolver>( nonlinearSolverParams );

    //----------------------------------------------------------------------------------------------------------------------------------------------//
    // initialize the column preconditioner which is a diagonal block preconditioner
    auto columnPreconditioner_db = linearSolver_db->getDatabase( "Preconditioner" );
    auto columnPreconditionerParams =
        AMP::make_shared<AMP::Solver::SolverStrategyParameters>( columnPreconditioner_db );
    columnPreconditionerParams->d_pOperator = linearThermalOxygenOperator;
    auto columnPreconditioner =
        AMP::make_shared<AMP::Solver::ColumnSolver>( columnPreconditionerParams );

    auto oxygenPreconditioner_db = columnPreconditioner_db->getDatabase( "oxygenPreconditioner" );
    auto oxygenPreconditionerParams =
        AMP::make_shared<AMP::Solver::SolverStrategyParameters>( oxygenPreconditioner_db );
    oxygenPreconditionerParams->d_pOperator = linearOxygenOperator;
    auto linearOxygenPreconditioner =
        AMP::make_shared<AMP::Solver::TrilinosMLSolver>( oxygenPreconditionerParams );

    auto thermalPreconditioner_db = columnPreconditioner_db->getDatabase( "thermalPreconditioner" );
    auto thermalPreconditionerParams =
        AMP::make_shared<AMP::Solver::SolverStrategyParameters>( thermalPreconditioner_db );
    thermalPreconditionerParams->d_pOperator = linearThermalOperator;
    auto linearThermalPreconditioner =
        AMP::make_shared<AMP::Solver::TrilinosMLSolver>( thermalPreconditionerParams );

    columnPreconditioner->append( linearOxygenPreconditioner );
    columnPreconditioner->append( linearThermalPreconditioner );

    //----------------------------------------------------------------------------------------------------------------------------------------------//
    // register the preconditioner with the Jacobian free Krylov solver
    auto linearSolver = nonlinearSolver->getKrylovSolver();

    linearSolver->setPreconditioner( columnPreconditioner );

    nonlinearThermalOxygenOperator->residual( rhsVec, solVec, resVec );
    double initialResidualNorm = resVec->L2Norm();

    AMP::pout << "Initial Residual Norm: " << initialResidualNorm << std::endl;

    nonlinearSolver->setZeroInitialGuess( false );

    nonlinearSolver->solve( rhsVec, solVec );

    nonlinearThermalOxygenOperator->residual( rhsVec, solVec, resVec );

    double finalResidualNorm = resVec->L2Norm();

    AMP::pout << "Final Residual Norm: " << finalResidualNorm << std::endl;

#ifdef USE_EXT_SILO
    manager->writeFile<AMP::Mesh::SiloIO>( exeName, 1 );
#endif

    if ( finalResidualNorm > initialResidualNorm * 1.0e-10 + 1.0e-05 ) {
        ITFAILS;
    } else {
        ut.passes( "PetscSNES Solver successfully solves a nonlinear thermal-oxygen diffusion "
                   "equation with JFNK, "
                   "FGMRES for Krylov, block diagonal preconditioning with ML solvers" );
    }
    ut.passes( exeName );
}

int main( int argc, char *argv[] )
{
    AMP::AMPManager::startup( argc, argv );
    AMP::UnitTest ut;

    AMP::AMP_MPI::initialize();
    std::vector<std::string> exeNames;
    //  exeNames.push_back("testPetscSNESSolver-NonlinearThermalOxygenDiffusion-1");
    exeNames.push_back( "testPetscSNESSolver-NonlinearThermalOxygenDiffusion-2" );

    for ( unsigned int i = 0; i < exeNames.size(); i++ )
        myTest( &ut, exeNames[i] );

    ut.report();

    int num_failed = ut.NumFailGlobal();
    AMP::AMPManager::shutdown();
    return num_failed;
}
