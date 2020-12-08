#include "AMP/ampmesh/Mesh.h"
#include "AMP/discretization/DOF_Manager.h"
#include "AMP/discretization/simpleDOF_Manager.h"
#include "AMP/materials/Material.h"
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
#include "AMP/solvers/ColumnSolver.h"
#include "AMP/solvers/trilinos/ml/TrilinosMLSolver.h"
#include "AMP/solvers/trilinos/nox/TrilinosNOXSolver.h"
#include "AMP/solvers/trilinos/nox/TrilinosNOXSolverParameters.h"
#include "AMP/utils/AMPManager.h"
#include "AMP/utils/AMP_MPI.h"
#include "AMP/utils/Database.h"
#include "AMP/utils/PIO.h"
#include "AMP/utils/UnitTest.h"
#include "AMP/utils/Utilities.h"
#include "AMP/utils/Writer.h"
#include "AMP/vectors/VectorBuilder.h"
#include <memory>

#include <iostream>
#include <string>


static void myTest( AMP::UnitTest *ut, const std::string &exeName )
{
    std::string input_file = "input_" + exeName;

    AMP::AMP_MPI globalComm( AMP_COMM_WORLD );

    size_t N_error0 = ut->NumFailLocal();

    auto input_db = AMP::Database::parseInputFile( input_file );
    input_db->print( AMP::plog );

    // Create the Mesh
    AMP_INSIST( input_db->keyExists( "Mesh" ), "Key ''Mesh'' is missing!" );
    auto mesh_db   = input_db->getDatabase( "Mesh" );
    auto mgrParams = std::make_shared<AMP::Mesh::MeshParameters>( mesh_db );
    mgrParams->setComm( AMP::AMP_MPI( AMP_COMM_WORLD ) );
    auto meshAdapter = AMP::Mesh::Mesh::buildMesh( mgrParams );

    // Create a DOF manager for a nodal vector
    int DOFsPerNode          = 1;
    int DOFsPerElement       = 8;
    int nodalGhostWidth      = 1;
    int gaussPointGhostWidth = 1;
    bool split               = true;
    auto nodalDofMap         = AMP::Discretization::simpleDOFManager::create(
        meshAdapter, AMP::Mesh::GeomType::Vertex, nodalGhostWidth, DOFsPerNode, split );
    auto gaussPointDofMap = AMP::Discretization::simpleDOFManager::create(
        meshAdapter, AMP::Mesh::GeomType::Volume, gaussPointGhostWidth, DOFsPerElement, split );

    // create a nonlinear BVP operator for nonlinear thermal diffusion
    AMP_INSIST( input_db->keyExists( "testNonlinearThermalOperator" ), "key missing!" );
    std::shared_ptr<AMP::Operator::ElementPhysicsModel> thermalTransportModel;
    auto nonlinearThermalOperator = std::dynamic_pointer_cast<AMP::Operator::NonlinearBVPOperator>(
        AMP::Operator::OperatorBuilder::createOperator(
            meshAdapter, "testNonlinearThermalOperator", input_db, thermalTransportModel ) );

    // initialize the input variable
    auto thermalVolumeOperator =
        std::dynamic_pointer_cast<AMP::Operator::DiffusionNonlinearFEOperator>(
            nonlinearThermalOperator->getVolumeOperator() );

    auto thermalVariable = thermalVolumeOperator->getOutputVariable();

    // create solution, rhs, and residual vectors
    auto solVec = AMP::LinearAlgebra::createVector( nodalDofMap, thermalVariable );
    auto rhsVec = AMP::LinearAlgebra::createVector( nodalDofMap, thermalVariable );
    auto resVec = AMP::LinearAlgebra::createVector( nodalDofMap, thermalVariable );

    // create the following shared pointers for ease of use
    AMP::LinearAlgebra::Vector::shared_ptr nullVec;

    // now construct the linear BVP operator for thermal
    AMP_INSIST( input_db->keyExists( "testLinearThermalOperator" ), "key missing!" );
    auto linearThermalOperator = std::dynamic_pointer_cast<AMP::Operator::LinearBVPOperator>(
        AMP::Operator::OperatorBuilder::createOperator(
            meshAdapter, "testLinearThermalOperator", input_db, thermalTransportModel ) );

    // Initial guess
    solVec->setToScalar( 400. );
    // solVec->setToScalar(557.2);
    std::cout << "initial guess norm = " << solVec->L2Norm() << "\n";
    nonlinearThermalOperator->modifyInitialSolutionVector( solVec );
    std::cout << "initial guess norm  after apply = " << solVec->L2Norm() << "\n";

    // CREATE THE NEUTRONICS SOURCE
    AMP_INSIST( input_db->keyExists( "NeutronicsOperator" ),
                "Key ''NeutronicsOperator'' is missing!" );
    auto neutronicsOp_db = input_db->getDatabase( "NeutronicsOperator" );
    auto neutronicsParams =
        std::make_shared<AMP::Operator::NeutronicsRhsParameters>( neutronicsOp_db );
    auto neutronicsOperator = std::make_shared<AMP::Operator::NeutronicsRhs>( neutronicsParams );

    auto SpecificPowerVar = neutronicsOperator->getOutputVariable();
    auto SpecificPowerVec = AMP::LinearAlgebra::createVector( gaussPointDofMap, SpecificPowerVar );

    // Set the SpecificPowerVec to the residual of the neutronicsOperator
    neutronicsOperator->apply( nullVec, SpecificPowerVec );
    SpecificPowerVec->scale( -1.0 );
    SpecificPowerVec->makeConsistent( AMP::LinearAlgebra::VectorData::ScatterType::CONSISTENT_SET );

    //  Integrate Nuclear Rhs over Desnity * GeomType::Volume
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
    sourceOperator->residual( nullVec, SpecificPowerVec, PowerInWattsVec );

    rhsVec->copyVector( PowerInWattsVec );

    nonlinearThermalOperator->modifyRHSvector( rhsVec );

    double initialRhsNorm = static_cast<double>( rhsVec->L2Norm() );
    std::cout << "rhs norm  after modifyRHSvector = " << initialRhsNorm << "\n";
    double expectedVal = 0.688628;
    if ( !AMP::Utilities::approx_equal( expectedVal, initialRhsNorm, 1e-5 ) )
        ut->failure( "the rhs norm after modifyRHSvector has changed." );

    // Get the solver databases
    auto nonlinearSolver_db = input_db->getDatabase( "NonlinearSolver" );
    auto linearSolver_db    = nonlinearSolver_db->getDatabase( "LinearSolver" );

    // Create the preconditioner
    auto thermalPreconditioner_db = linearSolver_db->getDatabase( "Preconditioner" );
    auto thermalPreconditionerParams =
        std::make_shared<AMP::Solver::SolverStrategyParameters>( thermalPreconditioner_db );
    thermalPreconditionerParams->d_pOperator = linearThermalOperator;
    auto linearThermalPreconditioner =
        std::make_shared<AMP::Solver::TrilinosMLSolver>( thermalPreconditionerParams );

    // Crete the solvers
    auto nonlinearSolverParams =
        std::make_shared<AMP::Solver::TrilinosNOXSolverParameters>( nonlinearSolver_db );
    nonlinearSolverParams->d_comm            = globalComm;
    nonlinearSolverParams->d_pOperator       = nonlinearThermalOperator;
    nonlinearSolverParams->d_pLinearOperator = nonlinearThermalOperator;
    nonlinearSolverParams->d_pInitialGuess   = solVec;
    nonlinearSolverParams->d_preconditioner  = linearThermalPreconditioner;
    auto nonlinearSolver =
        std::make_shared<AMP::Solver::TrilinosNOXSolver>( nonlinearSolverParams );

    nonlinearThermalOperator->residual( rhsVec, solVec, resVec );
    double initialResidualNorm = static_cast<double>( resVec->L2Norm() );

    AMP::pout << "Initial Residual Norm: " << initialResidualNorm << std::endl;
    expectedVal = 3625.84;
    if ( !AMP::Utilities::approx_equal( expectedVal, initialResidualNorm, 1e-5 ) ) {
        ut->failure( "the Initial Residual Norm has changed." );
    }

    nonlinearSolver->setZeroInitialGuess( false );

    nonlinearSolver->solve( rhsVec, solVec );

    solVec->makeConsistent( AMP::LinearAlgebra::VectorData::ScatterType::CONSISTENT_SET );
    resVec->makeConsistent( AMP::LinearAlgebra::VectorData::ScatterType::CONSISTENT_SET );

    nonlinearThermalOperator->residual( rhsVec, solVec, resVec );

    double finalResidualNorm = static_cast<double>( resVec->L2Norm() );
    double finalSolutionNorm = static_cast<double>( solVec->L2Norm() );
    double finalRhsNorm      = static_cast<double>( rhsVec->L2Norm() );

    std::cout << "Final Residual Norm: " << finalResidualNorm << std::endl;
    std::cout << "Final Solution Norm: " << finalSolutionNorm << std::endl;
    std::cout << "Final Rhs Norm: " << finalRhsNorm << std::endl;

    double ansSolutionNorm = 45431.3;
    if ( fabs( finalResidualNorm ) > 1e-9 )
        ut->failure( "the Final Residual is larger than the tolerance" );
    if ( !AMP::Utilities::approx_equal( ansSolutionNorm, finalSolutionNorm, 1e-5 ) )
        ut->failure( "the Final Solution Norm has changed." );
    if ( !AMP::Utilities::approx_equal( initialRhsNorm, finalRhsNorm, 1e-12 ) )
        ut->failure( "the Final Rhs Norm has changed." );

#ifdef USE_EXT_SILO
    auto siloWriter = AMP::Utilities::Writer::buildWriter( "Silo" );
    siloWriter->registerMesh( meshAdapter );
    siloWriter->registerVector( solVec, meshAdapter, AMP::Mesh::GeomType::Vertex, "Solution" );
    siloWriter->registerVector( resVec, meshAdapter, AMP::Mesh::GeomType::Vertex, "Residual" );
    siloWriter->writeFile( input_file, 0 );
#endif

    if ( N_error0 == ut->NumFailLocal() )
        ut->passes( exeName );
    else
        ut->failure( exeName );
}

int testTrilinosNOX_NonlinearThermal( int argc, char *argv[] )
{
    AMP::AMPManager::startup( argc, argv );
    AMP::UnitTest ut;

    std::vector<std::string> exeNames;
    exeNames.emplace_back( "testTrilinosNOX-NonlinearThermal-cylinder_MATPRO" );

    for ( const auto &exeName : exeNames )
        myTest( &ut, exeName );

    ut.report();

    int num_failed = ut.NumFailGlobal();
    AMP::AMPManager::shutdown();
    return num_failed;
}
