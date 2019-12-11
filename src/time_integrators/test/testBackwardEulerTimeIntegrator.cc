#include <memory>
#include <string>

#include "AMP/utils/AMPManager.h"
#include "AMP/utils/AMP_MPI.h"
#include "AMP/utils/Database.h"
#include "AMP/utils/PIO.h"
#include "AMP/utils/UnitTest.h"
#include "AMP/utils/Utilities.h"

#include "AMP/ampmesh/Mesh.h"
#include "AMP/discretization/DOF_Manager.h"
#include "AMP/discretization/simpleDOF_Manager.h"
#include "AMP/utils/Writer.h"
#include "AMP/vectors/VectorBuilder.h"

#include "AMP/materials/Material.h"

#include "AMP/vectors/Variable.h"
#include "AMP/vectors/Vector.h"

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

#define ITFAILS ut.failure( __LINE__ );
#define UNIT_TEST( a ) \
    if ( !( a ) )      \
        ut.failure( __LINE__ );

#define __PI__ 3.14159265
#define __INIT_FN__( x, y, z, t )                                                             \
    ( 750.0 + 10000.0 * ( 0.5 + x ) * ( 0.5 - x ) * ( 0.5 + y ) * ( 0.5 - y ) * ( 0.5 + z ) * \
                  ( 0.5 - z ) )

static void BackwardEulerTimeIntegrator( AMP::UnitTest *ut )
{
    std::string input_file = "input_testBackwardEulerTimeIntegrator";
    std::string log_file   = "output_testBackwardEulerTimeIntegrator";
    AMP::PIO::logOnlyNodeZero( log_file );

    // Read the input file
    auto input_db = AMP::Database::parseInputFile( input_file );

    // Get the Mesh database and create the mesh parameters
    std::shared_ptr<AMP::Database> database = input_db->getDatabase( "Mesh" );
    std::shared_ptr<AMP::Mesh::MeshParameters> params( new AMP::Mesh::MeshParameters( database ) );
    params->setComm( AMP::AMP_MPI( AMP_COMM_WORLD ) );

    // Create the meshes from the input database
    AMP::Mesh::Mesh::shared_ptr manager     = AMP::Mesh::Mesh::buildMesh( params );
    AMP::Mesh::Mesh::shared_ptr meshAdapter = manager->Subset( "Cube" );

    //--------------------------------------------------
    // Create a DOF manager for a nodal vector
    //--------------------------------------------------
    int DOFsPerNode          = 1;
    int DOFsPerElement       = 8;
    int nodalGhostWidth      = 1;
    int gaussPointGhostWidth = 1;
    bool split               = true;
    AMP::Discretization::DOFManager::shared_ptr nodalDofMap =
        AMP::Discretization::simpleDOFManager::create(
            meshAdapter, AMP::Mesh::GeomType::Vertex, nodalGhostWidth, DOFsPerNode, split );
    AMP::Discretization::DOFManager::shared_ptr gaussPointDofMap =
        AMP::Discretization::simpleDOFManager::create(
            meshAdapter, AMP::Mesh::GeomType::Volume, gaussPointGhostWidth, DOFsPerElement, split );

    // create a linear BVP operator
    std::shared_ptr<AMP::Operator::ElementPhysicsModel> elementModel;
    std::shared_ptr<AMP::Operator::LinearBVPOperator> diffusionOperator;
    diffusionOperator = std::dynamic_pointer_cast<AMP::Operator::LinearBVPOperator>(
        AMP::Operator::OperatorBuilder::createOperator(
            meshAdapter, "LinearOperator", input_db, elementModel ) );

    // create a mass linear BVP operator
    std::shared_ptr<AMP::Operator::ElementPhysicsModel> massElementModel;
    std::shared_ptr<AMP::Operator::LinearBVPOperator> massOperator;
    massOperator = std::dynamic_pointer_cast<AMP::Operator::LinearBVPOperator>(
        AMP::Operator::OperatorBuilder::createOperator(
            meshAdapter, "MassLinearOperator", input_db, massElementModel ) );

    //  create neutronics source
    AMP_INSIST( input_db->keyExists( "NeutronicsOperator" ),
                "Key ''NeutronicsOperator'' is missing!" );
    std::shared_ptr<AMP::Database> neutronicsOp_db = input_db->getDatabase( "NeutronicsOperator" );
    std::shared_ptr<AMP::Operator::NeutronicsRhsParameters> neutronicsParams(
        new AMP::Operator::NeutronicsRhsParameters( neutronicsOp_db ) );
    neutronicsParams->d_Mesh = meshAdapter;
    std::shared_ptr<AMP::Operator::NeutronicsRhs> neutronicsOperator(
        new AMP::Operator::NeutronicsRhs( neutronicsParams ) );

    AMP::LinearAlgebra::Variable::shared_ptr SpecificPowerVar =
        neutronicsOperator->getOutputVariable();
    AMP::LinearAlgebra::Vector::shared_ptr SpecificPowerVec =
        AMP::LinearAlgebra::createVector( gaussPointDofMap, SpecificPowerVar );

    // create the following shared pointers for ease of use
    AMP::LinearAlgebra::Vector::shared_ptr nullVec;

    neutronicsOperator->apply( nullVec, SpecificPowerVec );

    //  Integrate Nuclear Rhs over Density * GeomType::Volume //
    AMP_INSIST( input_db->keyExists( "VolumeIntegralOperator" ), "key missing!" );

    std::shared_ptr<AMP::Operator::ElementPhysicsModel> sourceTransportModel;
    std::shared_ptr<AMP::Operator::VolumeIntegralOperator> sourceOperator =
        std::dynamic_pointer_cast<AMP::Operator::VolumeIntegralOperator>(
            AMP::Operator::OperatorBuilder::createOperator(
                meshAdapter, "VolumeIntegralOperator", input_db, sourceTransportModel ) );

    // Create the power (heat source) vector.
    AMP::LinearAlgebra::Variable::shared_ptr powerInWattsVar = sourceOperator->getOutputVariable();
    AMP::LinearAlgebra::Vector::shared_ptr powerInWattsVec =
        AMP::LinearAlgebra::createVector( nodalDofMap, powerInWattsVar );
    powerInWattsVec->zero();

    // convert the vector of specific power to power for a given basis.
    sourceOperator->apply( SpecificPowerVec, powerInWattsVec );

    AMP::LinearAlgebra::Variable::shared_ptr outputVar = diffusionOperator->getOutputVariable();

    AMP::LinearAlgebra::Vector::shared_ptr initialCondition =
        AMP::LinearAlgebra::createVector( nodalDofMap, outputVar );
    AMP::LinearAlgebra::Vector::shared_ptr rhsVec =
        AMP::LinearAlgebra::createVector( nodalDofMap, outputVar );

    //----------------------------------------------------------------------------------------------------------------------------------------------//
    // set initial conditions, initialize created vectors
    int zeroGhostWidth = 0;
    AMP::Mesh::MeshIterator node =
        meshAdapter->getIterator( AMP::Mesh::GeomType::Vertex, zeroGhostWidth );
    AMP::Mesh::MeshIterator end_node = node.end();

    for ( ; node != end_node; ++node ) {
        std::vector<size_t> gid;
        nodalDofMap->getDOFs( node->globalID(), gid );

        double px = ( node->coord() )[0];
        double py = ( node->coord() )[1];
        double pz = ( node->coord() )[2];

        double val = __INIT_FN__( px, py, pz, 0 );
        for ( auto &elem : gid ) {
            initialCondition->setValueByGlobalID( elem, val );
        } // end for i
    }     // end for node

    diffusionOperator->modifyRHSvector( rhsVec );

    std::shared_ptr<AMP::Database> pcSolver_db = input_db->getDatabase( "Solver" );
    std::shared_ptr<AMP::Solver::SolverStrategyParameters> pcSolverParams(
        new AMP::Solver::SolverStrategyParameters( pcSolver_db ) );
    pcSolverParams->d_pOperator = diffusionOperator;
    std::shared_ptr<AMP::Solver::TrilinosMLSolver> pcSolver(
        new AMP::Solver::TrilinosMLSolver( pcSolverParams ) );

    std::shared_ptr<AMP::Database> timeIntegrator_db = input_db->getDatabase( "BDFTimeIntegrator" );
    std::shared_ptr<AMP::TimeIntegrator::ImplicitTimeIntegratorParameters> time_Params(
        new AMP::TimeIntegrator::ImplicitTimeIntegratorParameters( timeIntegrator_db ) );
    time_Params->d_pMassOperator = massOperator;
    time_Params->d_operator      = diffusionOperator;
    time_Params->d_solver        = pcSolver;

    time_Params->d_ic_vector = initialCondition;

    time_Params->d_pSourceTerm = rhsVec;
    time_Params->d_object_name = "ImplicitTimeIntegratorParameters";

    std::shared_ptr<AMP::TimeIntegrator::BackwardEulerTimeIntegrator> BDFTimeIntegrator(
        new AMP::TimeIntegrator::BackwardEulerTimeIntegrator( time_Params ) );

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
