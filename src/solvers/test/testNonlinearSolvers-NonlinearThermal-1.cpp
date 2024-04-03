#include "AMP/IO/PIO.h"
#include "AMP/IO/Writer.h"
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

#include <memory>

void myTest( AMP::UnitTest *ut, const std::string &inputName )
{
    std::string input_file = inputName;

    AMP::AMP_MPI globalComm( AMP_COMM_WORLD );
    size_t N_error0 = ut->NumFailLocal();
    auto input_db   = AMP::Database::parseInputFile( input_file );
    input_db->print( AMP::plog );

    //   Create the Mesh.
    AMP_INSIST( input_db->keyExists( "Mesh" ), "Key ''Mesh'' is missing!" );
    auto mesh_db   = input_db->getDatabase( "Mesh" );
    auto mgrParams = std::make_shared<AMP::Mesh::MeshParameters>( mesh_db );
    mgrParams->setComm( AMP::AMP_MPI( AMP_COMM_WORLD ) );
    auto meshAdapter = AMP::Mesh::MeshFactory::create( mgrParams );

    // Create a DOF manager for a nodal vector
    int DOFsPerNode          = 1;
    int DOFsPerElement       = 8;
    int nodalGhostWidth      = 1;
    int gaussPointGhostWidth = 1;
    bool split               = true;
    auto nodalDofMap         = AMP::Discretization::simpleDOFManager::create(
        meshAdapter, AMP::Mesh::GeomType::Vertex, nodalGhostWidth, DOFsPerNode, split );
    auto gaussPointDofMap = AMP::Discretization::simpleDOFManager::create(
        meshAdapter, AMP::Mesh::GeomType::Cell, gaussPointGhostWidth, DOFsPerElement, split );

    // create a nonlinear BVP operator for nonlinear thermal diffusion
    AMP_INSIST( input_db->keyExists( "testNonlinearThermalOperator" ), "key missing!" );
    std::shared_ptr<AMP::Operator::ElementPhysicsModel> thermalTransportModel;
    auto nonlinearThermalOperator = std::dynamic_pointer_cast<AMP::Operator::NonlinearBVPOperator>(
        AMP::Operator::OperatorBuilder::createOperator(
            meshAdapter, "testNonlinearThermalOperator", input_db, thermalTransportModel ) );

    auto thermalVariable = nonlinearThermalOperator->getOutputVariable();

    // create solution, rhs, and residual vectors
    auto solVec = AMP::LinearAlgebra::createVector( nodalDofMap, thermalVariable );
    auto rhsVec = AMP::LinearAlgebra::createVector( nodalDofMap, thermalVariable );
    auto resVec = AMP::LinearAlgebra::createVector( nodalDofMap, thermalVariable );

    // create the following shared pointers for ease of use
    AMP::LinearAlgebra::Vector::shared_ptr nullVec;

    // Initial guess
    solVec->setToScalar( 400. );
    double initialGuessNorm = static_cast<double>( solVec->L2Norm() );
    std::cout << "initial guess norm = " << initialGuessNorm << "\n";

    initialGuessNorm = static_cast<double>( solVec->L2Norm() );
    std::cout << "initial guess norm  after apply = " << initialGuessNorm << "\n";

    //  CREATE THE NEUTRONICS SOURCE  //
    AMP_INSIST( input_db->keyExists( "NeutronicsOperator" ),
                "Key ''NeutronicsOperator'' is missing!" );
    auto neutronicsOp_db = input_db->getDatabase( "NeutronicsOperator" );
    auto neutronicsParams =
        std::make_shared<AMP::Operator::NeutronicsRhsParameters>( neutronicsOp_db );
    auto neutronicsOperator = std::make_shared<AMP::Operator::NeutronicsRhs>( neutronicsParams );
    auto SpecificPowerVar   = neutronicsOperator->getOutputVariable();
    auto SpecificPowerVec = AMP::LinearAlgebra::createVector( gaussPointDofMap, SpecificPowerVar );

    neutronicsOperator->apply( nullVec, SpecificPowerVec );

    //  Integrate Nuclear Rhs over Desnity * GeomType::Cell //
    AMP_INSIST( input_db->keyExists( "VolumeIntegralOperator" ), "key missing!" );
    std::shared_ptr<AMP::Operator::ElementPhysicsModel> stransportModel;
    auto sourceOperator = std::dynamic_pointer_cast<AMP::Operator::VolumeIntegralOperator>(
        AMP::Operator::OperatorBuilder::createOperator(
            meshAdapter, "VolumeIntegralOperator", input_db, stransportModel ) );

    // Create the power (heat source) vector.
    auto PowerInWattsVar = sourceOperator->getOutputVariable();
    auto PowerInWattsVec = AMP::LinearAlgebra::createVector( nodalDofMap, PowerInWattsVar );
    PowerInWattsVec->zero();

    // convert the vector of specific power to power for a given basis.
    sourceOperator->apply( SpecificPowerVec, PowerInWattsVec );

    rhsVec->copyVector( PowerInWattsVec );

    nonlinearThermalOperator->modifyInitialSolutionVector( solVec );
    nonlinearThermalOperator->modifyRHSvector( rhsVec );

    // Create the solver
    auto nonlinearSolver = AMP::Solver::Test::buildSolver(
        "NonlinearSolver", input_db, globalComm, solVec, nonlinearThermalOperator );

    nonlinearSolver->apply( rhsVec, solVec );

    nonlinearThermalOperator->residual( rhsVec, solVec, resVec );

    double finalResidualNorm = static_cast<double>( resVec->L2Norm() );
    double finalSolutionNorm = static_cast<double>( solVec->L2Norm() );

    AMP::pout << "Final Residual Norm: " << finalResidualNorm << std::endl;
    AMP::pout << "Final Solution Norm: " << finalSolutionNorm << std::endl;

    if ( fabs( finalResidualNorm ) > 1e-8 )
        ut->failure( "the Final Residual is larger than the tolerance" );
    if ( !AMP::Utilities::approx_equal( 45431.3, finalSolutionNorm, 1e-5 ) )
        ut->failure( "the Final Solution Norm has changed." );

#ifdef AMP_USE_SILO
    auto siloWriter = AMP::IO::Writer::buildWriter( "Silo" );
    siloWriter->registerMesh( meshAdapter );
    siloWriter->registerVector( solVec, meshAdapter, AMP::Mesh::GeomType::Vertex, "Solution" );
    siloWriter->registerVector( resVec, meshAdapter, AMP::Mesh::GeomType::Vertex, "Residual" );
    siloWriter->writeFile( input_file, 0 );
#endif

    if ( N_error0 == ut->NumFailLocal() )
        ut->passes( inputName );
    else
        ut->failure( inputName );
}

int main( int argc, char *argv[] )
{
    AMP::AMPManager::startup( argc, argv );
    AMP::UnitTest ut;
    AMP::Solver::registerSolverFactories();

    std::vector<std::string> inputNames;

    if ( argc > 1 ) {
        inputNames.push_back( argv[1] );
    } else {

#ifdef AMP_USE_PETSC
        inputNames.emplace_back( "input_PetscSNESSolver-NonlinearThermal-cylinder_MATPRO" );
    #ifdef AMP_USE_HYPRE
        inputNames.emplace_back(
            "input_PetscSNESSolver-BoomerAMG-NonlinearThermal-cylinder_MATPRO" );
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
    }

    for ( auto &inputName : inputNames )
        myTest( &ut, inputName );

    ut.report();

    int num_failed = ut.NumFailGlobal();
    AMP::AMPManager::shutdown();
    return num_failed;
}
