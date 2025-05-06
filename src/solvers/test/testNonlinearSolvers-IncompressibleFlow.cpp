#include "AMP/IO/PIO.h"
#include "AMP/discretization/DOF_Manager.h"
#include "AMP/discretization/simpleDOF_Manager.h"
#include "AMP/mesh/Mesh.h"
#include "AMP/mesh/MeshFactory.h"
#include "AMP/mesh/MeshParameters.h"
#include "AMP/operators/BVPOperatorParameters.h"
#include "AMP/operators/ColumnOperator.h"
#include "AMP/operators/NonlinearBVPOperator.h"
#include "AMP/operators/OperatorBuilder.h"
#include "AMP/operators/flow/NavierStokesLSWFFEOperator.h"
#include "AMP/operators/flow/NavierStokesLSWFLinearFEOperator.h"
#include "AMP/operators/libmesh/VolumeIntegralOperator.h"
#include "AMP/solvers/SolverFactory.h"
#include "AMP/solvers/SolverStrategy.h"
#include "AMP/solvers/SolverStrategyParameters.h"
#include "AMP/solvers/testHelpers/SolverTestParameters.h"
#include "AMP/utils/AMPManager.h"
#include "AMP/utils/AMP_MPI.h"
#include "AMP/utils/Database.h"
#include "AMP/utils/UnitTest.h"
#include "AMP/utils/Utilities.h"
#include "AMP/vectors/VectorBuilder.h"

#include <iostream>
#include <memory>
#include <string>

#include "AMP/solvers/testHelpers/testSolverHelpers.h"

static void myTest( AMP::UnitTest *ut, const std::string &inputName )
{
    std::string input_file = inputName;

    AMP::AMP_MPI globalComm( AMP_COMM_WORLD );
    size_t N_error0 = ut->NumFailLocal();
    auto input_db   = AMP::Database::parseInputFile( input_file );
    input_db->print( AMP::plog );

    // create the Mesh
    const auto meshAdapter = createMesh( input_db );

    // Create a DOF manager for a nodal vector
    int DOFsPerNode     = 10;
    int nodalGhostWidth = 1;
    bool split          = true;
    auto nodalDofMap    = AMP::Discretization::simpleDOFManager::create(
        meshAdapter, AMP::Mesh::GeomType::Vertex, nodalGhostWidth, DOFsPerNode, split );

    // create a nonlinear BVP operator for nonlinear flow
    AMP_INSIST( input_db->keyExists( "NonlinearFlowOperator" ), "key missing!" );
    std::shared_ptr<AMP::Operator::ElementPhysicsModel> flowTransportModel;
    auto nonlinearFlowOperator = std::dynamic_pointer_cast<AMP::Operator::NonlinearBVPOperator>(
        AMP::Operator::OperatorBuilder::createOperator(
            meshAdapter, "NonlinearFlowOperator", input_db, flowTransportModel ) );

    // initialize the input variable
    auto flowVolumeOperator = std::dynamic_pointer_cast<AMP::Operator::NavierStokesLSWFFEOperator>(
        nonlinearFlowOperator->getVolumeOperator() );

    auto flowVariable = flowVolumeOperator->getOutputVariable();

    // create solution, rhs, and residual vectors
    auto solVec = AMP::LinearAlgebra::createVector( nodalDofMap, flowVariable );
    auto rhsVec = AMP::LinearAlgebra::createVector( nodalDofMap, flowVariable );
    auto resVec = AMP::LinearAlgebra::createVector( nodalDofMap, flowVariable );

    // create the following shared pointers for ease of use

    // Initial guess
    solVec->setToScalar( 0. );
    std::cout << "initial guess norm = " << solVec->L2Norm() << "\n";
    nonlinearFlowOperator->modifyInitialSolutionVector( solVec );
    std::cout << "initial guess norm  after apply = " << solVec->L2Norm() << "\n";

    rhsVec->zero();

    nonlinearFlowOperator->modifyRHSvector( rhsVec );

    double initialRhsNorm = static_cast<double>( rhsVec->L2Norm() );
    std::cout << "rhs norm  after modifyRHSvector = " << initialRhsNorm << "\n";
    double expectedVal = 0.;
    if ( !AMP::Utilities::approx_equal( expectedVal, initialRhsNorm, 1e-5 ) )
        ut->failure( "the rhs norm after modifyRHSvector has changed." );

    // Create the solver
    auto nonlinearSolver = AMP::Solver::Test::buildSolver(
        "NonlinearSolver", input_db, globalComm, solVec, nonlinearFlowOperator );

#if 0
    // now construct the linear BVP operator for flow
    AMP_INSIST( input_db->keyExists( "LinearFlowOperator" ), "key missing!" );
    auto linearFlowOperator = std::dynamic_pointer_cast<AMP::Operator::LinearBVPOperator>(
        AMP::Operator::OperatorBuilder::createOperator(
            meshAdapter, "LinearFlowOperator", input_db, flowTransportModel ) );

    // Get the solver databases
    auto nonlinearSolver_db = input_db->getDatabase( "NonlinearSolver" );



    auto linearSolver_db    = nonlinearSolver_db->getDatabase( "LinearSolver" );

    // Create the preconditioner
    auto flowPreconditioner_db = linearSolver_db->getDatabase( "Preconditioner" );
    auto flowPreconditionerParams =
        std::make_shared<AMP::Solver::SolverStrategyParameters>( flowPreconditioner_db );
    flowPreconditionerParams->d_pOperator = linearFlowOperator;
    auto linearFlowPreconditioner =
        std::make_shared<AMP::Solver::TrilinosMLSolver>( flowPreconditionerParams );

    // initialize the linear solver
    auto linearSolverParams =
        std::make_shared<AMP::Solver::SolverStrategyParameters>( linearSolver_db );
    linearSolverParams->d_pOperator     = linearFlowOperator;
    linearSolverParams->d_comm          = globalComm;
    linearSolverParams->d_pNestedSolver = linearFlowPreconditioner;
    //    std::shared_ptr<AMP::Solver::PetscKrylovSolver> linearSolver(new
    //    AMP::Solver::PetscKrylovSolver(linearSolverParams));

    // Crete the solvers
    auto nonlinearSolverParams =
        std::make_shared<AMP::Solver::SolverStrategyParameters>( nonlinearSolver_db );
    nonlinearSolverParams->d_comm      = globalComm;
    nonlinearSolverParams->d_pOperator = nonlinearFlowOperator;
    //    nonlinearSolverParams->d_pKrylovSolver = linearSolver;
    nonlinearSolverParams->d_pInitialGuess = solVec;
    auto nonlinearSolver = std::make_shared<AMP::Solver::PetscSNESSolver>( nonlinearSolverParams );

    auto linearSolver = nonlinearSolver->getKrylovSolver();
    linearSolver->setNestedSolver( linearFlowPreconditioner );
#endif

    nonlinearFlowOperator->residual( rhsVec, solVec, resVec );
    auto initialResidualNorm = static_cast<double>( resVec->L2Norm() );

    AMP::pout << "Initial Residual Norm: " << initialResidualNorm << std::endl;
    expectedVal = 3625.84;
    if ( !AMP::Utilities::approx_equal( expectedVal, initialResidualNorm, 1e-5 ) ) {
        ut->failure( "the Initial Residual Norm has changed." );
    }

    nonlinearSolver->setZeroInitialGuess( false );

    nonlinearSolver->apply( rhsVec, solVec );

    nonlinearFlowOperator->residual( rhsVec, solVec, resVec );

    auto finalResidualNorm = static_cast<double>( resVec->L2Norm() );
    auto finalSolutionNorm = static_cast<double>( solVec->L2Norm() );
    auto finalRhsNorm      = static_cast<double>( rhsVec->L2Norm() );

    std::cout << "Final Residual Norm: " << finalResidualNorm << std::endl;
    std::cout << "Final Solution Norm: " << finalSolutionNorm << std::endl;
    std::cout << "Final Rhs Norm: " << finalRhsNorm << std::endl;

    if ( fabs( finalResidualNorm ) > 1e-9 )
        ut->failure( "the Final Residual is larger than the tolerance" );
    if ( !AMP::Utilities::approx_equal( 45431.3, finalSolutionNorm, 1e-5 ) )
        ut->failure( "the Final Solution Norm has changed." );
    if ( !AMP::Utilities::approx_equal( initialRhsNorm, finalRhsNorm, 1e-12 ) )
        ut->failure( "the Final Rhs Norm has changed." );

    if ( N_error0 == ut->NumFailLocal() )
        ut->passes( inputName );
    else
        ut->failure( inputName );
}

int main( int argc, char *argv[] )
{
    AMP::AMPManager::startup( argc, argv );
    AMP::UnitTest ut;

    std::vector<std::string> inputNames;

    if ( argc > 1 ) {
        inputNames.push_back( argv[1] );
    } else {
#ifdef AMP_USE_PETSC
    #ifdef AMP_USE_TRILINOS_ML
        inputNames.emplace_back( "input_testPetscSNESSolver-IncompressibleFlow-1" );
    #endif
#endif
    }
    for ( auto &inputName : inputNames )
        myTest( &ut, inputName );

    ut.report();

    int num_failed = ut.NumFailGlobal();
    AMP::AMPManager::shutdown();
    return num_failed;
}
