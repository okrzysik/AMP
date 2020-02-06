#include "AMP/ampmesh/Mesh.h"
#include "AMP/discretization/DOF_Manager.h"
#include "AMP/discretization/simpleDOF_Manager.h"
#include "AMP/materials/Material.h"
#include "AMP/operators/LinearBVPOperator.h"
#include "AMP/operators/NeutronicsRhs.h"
#include "AMP/operators/OperatorBuilder.h"
#include "AMP/operators/diffusion/DiffusionLinearElement.h"
#include "AMP/operators/diffusion/DiffusionLinearFEOperator.h"
#include "AMP/operators/diffusion/DiffusionTransportModel.h"
#include "AMP/operators/libmesh/MassLinearElement.h"
#include "AMP/operators/libmesh/MassLinearFEOperator.h"
#include "AMP/operators/libmesh/VolumeIntegralOperator.h"
#include "AMP/solvers/trilinos/ml/TrilinosMLSolver.h"
#include "AMP/time_integrators/BackwardEulerTimeIntegrator.h"
#include "AMP/time_integrators/BackwardEulerTimeOperator.h"
#include "AMP/time_integrators/ImplicitTimeIntegratorParameters.h"
#include "AMP/utils/AMPManager.h"
#include "AMP/utils/AMP_MPI.h"
#include "AMP/utils/Database.h"
#include "AMP/utils/PIO.h"
#include "AMP/utils/UnitTest.h"
#include "AMP/utils/Utilities.h"
#include "AMP/utils/Writer.h"
#include "AMP/vectors/Variable.h"
#include "AMP/vectors/Vector.h"
#include "AMP/vectors/VectorBuilder.h"

#include <memory>
#include <string>


static inline double fun( double x, double y, double z, double t )
{
    return ( 750.0 + 10000.0 * ( 0.5 + x ) * ( 0.5 - x ) * ( 0.5 + y ) * ( 0.5 - y ) * ( 0.5 + z ) *
                         ( 0.5 - z ) );
}


static void BackwardEulerTimeIntegrator( AMP::UnitTest *ut )
{
    std::string input_file = "input_testBackwardEulerTimeIntegrator";
    std::string log_file   = "output_testBackwardEulerTimeIntegrator";
    AMP::PIO::logOnlyNodeZero( log_file );

    // Read the input file
    auto input_db = AMP::Database::parseInputFile( input_file );

    // Get the Mesh database and create the mesh parameters
    auto database = input_db->getDatabase( "Mesh" );
    auto params   = std::make_shared<AMP::Mesh::MeshParameters>( database );
    params->setComm( AMP::AMP_MPI( AMP_COMM_WORLD ) );

    // Create the meshes from the input database
    auto manager     = AMP::Mesh::Mesh::buildMesh( params );
    auto meshAdapter = manager->Subset( "Cube" );

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

    // create a linear BVP operator
    std::shared_ptr<AMP::Operator::ElementPhysicsModel> elementModel;
    auto diffusionOperator = std::dynamic_pointer_cast<AMP::Operator::LinearBVPOperator>(
        AMP::Operator::OperatorBuilder::createOperator(
            meshAdapter, "LinearOperator", input_db, elementModel ) );

    // create a mass linear BVP operator
    std::shared_ptr<AMP::Operator::ElementPhysicsModel> massElementModel;
    auto massOperator = std::dynamic_pointer_cast<AMP::Operator::LinearBVPOperator>(
        AMP::Operator::OperatorBuilder::createOperator(
            meshAdapter, "MassLinearOperator", input_db, massElementModel ) );

    //  create neutronics source
    AMP_INSIST( input_db->keyExists( "NeutronicsOperator" ),
                "Key ''NeutronicsOperator'' is missing!" );
    auto neutronicsOp_db = input_db->getDatabase( "NeutronicsOperator" );
    auto neutronicsParams =
        std::make_shared<AMP::Operator::NeutronicsRhsParameters>( neutronicsOp_db );
    neutronicsParams->d_Mesh = meshAdapter;
    auto neutronicsOperator  = std::make_shared<AMP::Operator::NeutronicsRhs>( neutronicsParams );

    auto SpecificPowerVar = neutronicsOperator->getOutputVariable();
    auto SpecificPowerVec = AMP::LinearAlgebra::createVector( gaussPointDofMap, SpecificPowerVar );

    // create the following shared pointers for ease of use
    AMP::LinearAlgebra::Vector::shared_ptr nullVec;

    neutronicsOperator->apply( nullVec, SpecificPowerVec );

    //  Integrate Nuclear Rhs over Density * GeomType::Volume //
    AMP_INSIST( input_db->keyExists( "VolumeIntegralOperator" ), "key missing!" );

    std::shared_ptr<AMP::Operator::ElementPhysicsModel> sourceTransportModel;
    auto sourceOperator = std::dynamic_pointer_cast<AMP::Operator::VolumeIntegralOperator>(
        AMP::Operator::OperatorBuilder::createOperator(
            meshAdapter, "VolumeIntegralOperator", input_db, sourceTransportModel ) );

    // Create the power (heat source) vector.
    auto powerInWattsVar = sourceOperator->getOutputVariable();
    auto powerInWattsVec = AMP::LinearAlgebra::createVector( nodalDofMap, powerInWattsVar );
    powerInWattsVec->zero();

    // convert the vector of specific power to power for a given basis.
    sourceOperator->apply( SpecificPowerVec, powerInWattsVec );

    auto outputVar = diffusionOperator->getOutputVariable();

    auto initialCondition = AMP::LinearAlgebra::createVector( nodalDofMap, outputVar );
    auto rhsVec           = AMP::LinearAlgebra::createVector( nodalDofMap, outputVar );

    //----------------------------------------------------------------------------------------------------------------------------------------------//
    // set initial conditions, initialize created vectors
    int zeroGhostWidth = 0;
    auto node          = meshAdapter->getIterator( AMP::Mesh::GeomType::Vertex, zeroGhostWidth );
    auto end_node      = node.end();
    for ( ; node != end_node; ++node ) {
        std::vector<size_t> gid;
        nodalDofMap->getDOFs( node->globalID(), gid );

        double px = ( node->coord() )[0];
        double py = ( node->coord() )[1];
        double pz = ( node->coord() )[2];

        double val = fun( px, py, pz, 0 );
        for ( auto &elem : gid ) {
            initialCondition->setValueByGlobalID( elem, val );
        } // end for i
    }     // end for node

    diffusionOperator->modifyRHSvector( rhsVec );

    auto pcSolver_db    = input_db->getDatabase( "Solver" );
    auto pcSolverParams = std::make_shared<AMP::Solver::SolverStrategyParameters>( pcSolver_db );
    pcSolverParams->d_pOperator = diffusionOperator;
    auto pcSolver               = std::make_shared<AMP::Solver::TrilinosMLSolver>( pcSolverParams );

    auto timeIntegrator_db = input_db->getDatabase( "BDFTimeIntegrator" );
    auto time_Params = std::make_shared<AMP::TimeIntegrator::ImplicitTimeIntegratorParameters>(
        timeIntegrator_db );
    time_Params->d_pMassOperator = massOperator;
    time_Params->d_operator      = diffusionOperator;
    time_Params->d_solver        = pcSolver;

    time_Params->d_ic_vector = initialCondition;

    time_Params->d_pSourceTerm = rhsVec;
    time_Params->d_object_name = "ImplicitTimeIntegratorParameters";

    auto BDFTimeIntegrator =
        std::make_shared<AMP::TimeIntegrator::BackwardEulerTimeIntegrator>( time_Params );

    if ( BDFTimeIntegrator.get() == nullptr ) {
        ut->failure( "Testing BDFTimeIntegrator's constructor" );
    } else {
        ut->passes( "Tested BDFTimeIntegrator's constructor" );
    }

    double current_time = 0, max, min;
    int j               = 0;
    while ( BDFTimeIntegrator->getCurrentTime() < BDFTimeIntegrator->getFinalTime() ) {
        BDFTimeIntegrator->advanceSolution( BDFTimeIntegrator->getCurrentDt(), false );
        current_time = BDFTimeIntegrator->getCurrentTime();

        std::cout << j++ << "-th timestep" << std::endl;

        max = BDFTimeIntegrator->getCurrentSolution()->max();
        min = BDFTimeIntegrator->getCurrentSolution()->min();

        std::cout << "current_time = " << current_time << std::endl;
        std::cout << "max val of the current solution = " << max << std::endl;
        std::cout << "min val of the current solution = " << min << std::endl;
    }

    if ( ut->NumFailLocal() == 0 ) {
        ut->passes( "test Backward Euler Time Intgrator successful" );
    }
}


//---------------------------------------------------------------------------//

int testBackwardEulerTimeIntegrator( int argc, char *argv[] )
{
    AMP::AMPManager::startup( argc, argv );
    AMP::UnitTest ut;

    BackwardEulerTimeIntegrator( &ut );

    ut.report();

    int num_failed = ut.NumFailGlobal();
    AMP::AMPManager::shutdown();
    return num_failed;
}
