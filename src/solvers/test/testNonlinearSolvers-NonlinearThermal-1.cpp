#include "AMP/IO/PIO.h"
#include "AMP/discretization/DOF_Manager.h"
#include "AMP/discretization/simpleDOF_Manager.h"
#include "AMP/mesh/Mesh.h"
#include "AMP/mesh/MeshFactory.h"
#include "AMP/mesh/MeshParameters.h"
#include "AMP/operators/BVPOperatorParameters.h"
#include "AMP/operators/ColumnOperator.h"
#include "AMP/operators/LinearBVPOperator.h"
#include "AMP/operators/NeutronicsRhs.h"
#include "AMP/operators/NonlinearBVPOperator.h"
#include "AMP/operators/OperatorBuilder.h"
#include "AMP/operators/boundary/DirichletVectorCorrection.h"
#include "AMP/operators/diffusion/DiffusionLinearFEOperator.h"
#include "AMP/operators/diffusion/DiffusionNonlinearFEOperator.h"
#include "AMP/operators/libmesh/VolumeIntegralOperator.h"
#include "AMP/operators/mechanics/MechanicsLinearFEOperator.h"
#include "AMP/operators/mechanics/MechanicsNonlinearFEOperator.h"
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

#include "ProfilerApp.h"

#include <memory>

#include "testSolverHelpers.h"

void myTest( AMP::UnitTest *ut, const std::string &inputName )
{
    PROFILE2( inputName );

    std::string input_file = inputName;

    AMP::AMP_MPI globalComm( AMP_COMM_WORLD );
    size_t N_error0 = ut->NumFailLocal();
    auto input_db   = AMP::Database::parseInputFile( input_file );
    input_db->print( AMP::plog );

    // create the Mesh
    const auto meshAdapter = createMesh( input_db );

    // create a power source from neutronics
    auto PowerInWattsVec = constructNeutronicsPowerSource( input_db, meshAdapter );

    // create a nonlinear BVP operator for nonlinear thermal diffusion
    AMP_INSIST( input_db->keyExists( "testNonlinearThermalOperator" ), "key missing!" );
    std::shared_ptr<AMP::Operator::ElementPhysicsModel> thermalTransportModel;
    auto nonlinearThermalOperator = std::dynamic_pointer_cast<AMP::Operator::NonlinearBVPOperator>(
        AMP::Operator::OperatorBuilder::createOperator(
            meshAdapter, "testNonlinearThermalOperator", input_db, thermalTransportModel ) );

    auto thermalVariable = nonlinearThermalOperator->getOutputVariable();
    auto nodalDofMap     = PowerInWattsVec->getDOFManager();

    // create solution, rhs, and residual vectors
    auto solVec = AMP::LinearAlgebra::createVector( nodalDofMap, thermalVariable );
    auto rhsVec = AMP::LinearAlgebra::createVector( nodalDofMap, thermalVariable );
    auto resVec = AMP::LinearAlgebra::createVector( nodalDofMap, thermalVariable );

    // Initial guess
    solVec->setToScalar( 400. );
    auto initialGuessNorm = static_cast<double>( solVec->L2Norm() );
    AMP::pout << "initial guess norm = " << initialGuessNorm << "\n";

    initialGuessNorm = static_cast<double>( solVec->L2Norm() );
    AMP::pout << "initial guess norm  after apply = " << initialGuessNorm << "\n";

    rhsVec->copyVector( PowerInWattsVec );

    nonlinearThermalOperator->modifyInitialSolutionVector( solVec );
    nonlinearThermalOperator->modifyRHSvector( rhsVec );

    // Create the solver
    auto nonlinearSolver = AMP::Solver::Test::buildSolver(
        "NonlinearSolver", input_db, globalComm, solVec, nonlinearThermalOperator );

    nonlinearSolver->setZeroInitialGuess( false );

    nonlinearSolver->apply( rhsVec, solVec );

    nonlinearThermalOperator->residual( rhsVec, solVec, resVec );

    auto finalResidualNorm = static_cast<double>( resVec->L2Norm() );
    auto finalSolutionNorm = static_cast<double>( solVec->L2Norm() );

    AMP::pout << "Final Residual Norm: " << finalResidualNorm << std::endl;
    AMP::pout << "Final Solution Norm: " << finalSolutionNorm << std::endl;

    auto referenceSolutionNorm = input_db->getScalar<double>( "referenceSolutionNorm" );

    if ( fabs( finalResidualNorm ) > 1e-8 )
        ut->failure( "the Final Residual is larger than the tolerance" );
    if ( !AMP::Utilities::approx_equal( referenceSolutionNorm, finalSolutionNorm, 1e-5 ) )
        ut->failure( "the Final Solution Norm has changed." );
    if ( N_error0 == ut->NumFailLocal() )
        ut->passes( inputName );
    else
        ut->failure( inputName );
}

int main( int argc, char *argv[] )
{
    AMP::AMPManager::startup( argc, argv );
    AMP::UnitTest ut;
    PROFILE_ENABLE( 8 );

    std::vector<std::string> inputNames;

    if ( argc > 1 ) {
        inputNames.push_back( argv[1] );
    } else {

#ifdef AMP_USE_PETSC
        inputNames.emplace_back( "input_PetscSNESSolver-NonlinearThermal-cylinder_MATPRO" );
        inputNames.emplace_back( "input_PetscSNESSolver-NonlinearThermal-cylinder_MATPRO2" );

    #ifdef AMP_USE_HYPRE
        inputNames.emplace_back(
            "input_PetscSNESSolver-BoomerAMG-NonlinearThermal-cylinder_MATPRO" );
        inputNames.emplace_back(
            "input_PetscSNESSolver-BoomerAMG-NonlinearThermal-cylinder_MATPRO2" );
    #endif
    #ifdef AMP_USE_TRILINOS_ML
        inputNames.emplace_back(
            "input_PetscSNESSolver-TrilinosML-NonlinearThermal-cylinder_MATPRO" );
    #endif
#endif

#ifdef AMP_USE_HYPRE
        inputNames.emplace_back( "input_NKA-BoomerAMG-NonlinearThermal-cylinder_MATPRO" );
#endif

#ifdef AMP_USE_TRILINOS_ML
        inputNames.emplace_back( "input_NKA-TrilinosML-NonlinearThermal-cylinder_MATPRO" );
#endif

#ifdef AMP_USE_TRILINOS_NOX
    #ifdef AMP_USE_TRILINOS_ML
        inputNames.emplace_back( "input_TrilinosNOX-TrilinosML-NonlinearThermal-cylinder_MATPRO" );
    #endif
#endif
    }

    for ( auto &inputName : inputNames )
        myTest( &ut, inputName );

    ut.report();

    PROFILE_SAVE( "testNonlinearSolvers-NonlinearThermal-1" );

    int num_failed = ut.NumFailGlobal();
    AMP::AMPManager::shutdown();
    return num_failed;
}
