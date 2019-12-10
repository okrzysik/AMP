#include "AMP/ampmesh/Mesh.h"
#include "AMP/discretization/DOF_Manager.h"
#include "AMP/discretization/simpleDOF_Manager.h"
#include "AMP/materials/Material.h"
#include "AMP/operators/LinearBVPOperator.h"
#include "AMP/operators/NeutronicsRhs.h"
#include "AMP/operators/NonlinearBVPOperator.h"
#include "AMP/operators/diffusion/DiffusionLinearFEOperator.h"
#include "AMP/operators/libmesh/MassLinearFEOperator.h"
#include "AMP/operators/libmesh/VolumeIntegralOperator.h"
#include "AMP/solvers/trilinos/ml/TrilinosMLSolver.h"
#include "AMP/time_integrators/sundials/IDATimeIntegrator.h"
#include "AMP/time_integrators/sundials/IDATimeOperator.h"
#include "AMP/utils/AMPManager.h"
#include "AMP/utils/AMP_MPI.h"
#include "AMP/utils/Database.h"
#include "AMP/utils/PIO.h"
#include "AMP/utils/UnitTest.h"
#include "AMP/utils/Utilities.h"
#include "AMP/utils/Writer.h"
#include "AMP/utils/shared_ptr.h"
#include "AMP/vectors/Variable.h"
#include "AMP/vectors/Vector.h"
#include "AMP/vectors/VectorBuilder.h"

#include <string>


#define ITFAILS ut.failure( __LINE__ );
#define UNIT_TEST( a ) \
    if ( !( a ) )      \
        ut.failure( __LINE__ );

//---------------------------------------------------------------------------//
// TESTS
//---------------------------------------------------------------------------//

#define __PI__ 3.14159265

//#define __INIT_FN__(x, y, z, t) ( exp(- __PI__ * __PI__ * t) * sin(__PI__ * x) * sin(__PI__ * y) *
// sin(__PI__ * z) )
//#define __INIT_FN__(x,y,z,t) ( exp(-0.015 *  __PI__ * __PI__ * t) * cos(0.1 * __PI__ * x) *
// cos(0.1 * __PI__ * y) *
// cos(0.05 * __PI__ * z) )
#define __INIT_FN__( x, y, z, t )                                                             \
    ( 750.0 + 10000.0 * ( 0.5 + x ) * ( 0.5 - x ) * ( 0.5 + y ) * ( 0.5 - y ) * ( 0.5 + z ) * \
                  ( 0.5 - z ) )

static void IDATimeIntegratorTest( AMP::UnitTest *ut )
{
    std::string input_file = "input_testIDA-NonlinearBVPOperator-1";
    std::string log_file   = "output_testIDA-NonlinearBVPOperator-1";

    AMP::PIO::logOnlyNodeZero( log_file );


    AMP::AMP_MPI globalComm( AMP_COMM_WORLD );
    auto input_db = AMP::Database::parseInputFile( input_file );
    input_db->print( AMP::plog );

    AMP_INSIST( input_db->keyExists( "Mesh" ), "Key ''Mesh'' is missing!" );
    AMP::shared_ptr<AMP::Database> mesh_db = input_db->getDatabase( "Mesh" );
    AMP::shared_ptr<AMP::Mesh::MeshParameters> mgrParams(
        new AMP::Mesh::MeshParameters( mesh_db ) );
    mgrParams->setComm( AMP::AMP_MPI( AMP_COMM_WORLD ) );
    AMP::shared_ptr<AMP::Mesh::Mesh> meshAdapter = AMP::Mesh::Mesh::buildMesh( mgrParams );

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

    //----------------------------------------------------------------------------------------------------------------------------------------------//
    // create a nonlinear BVP operator for nonlinear BVP operator
    AMP_INSIST( input_db->keyExists( "NonlinearOperator" ), "key missing!" );

    AMP::shared_ptr<AMP::Operator::ElementPhysicsModel> elementModel;
    AMP::shared_ptr<AMP::Operator::NonlinearBVPOperator> nonlinearOperator =
        AMP::dynamic_pointer_cast<AMP::Operator::NonlinearBVPOperator>(
            AMP::Operator::OperatorBuilder::createOperator(
                meshAdapter, "NonlinearOperator", input_db, elementModel ) );

    AMP::LinearAlgebra::Variable::shared_ptr outputVar = nonlinearOperator->getOutputVariable();
    // ---------------------------------------------------------------------------------------
    // create a linear BVP operator
    AMP::shared_ptr<AMP::Operator::LinearBVPOperator> linearOperator =
        AMP::dynamic_pointer_cast<AMP::Operator::LinearBVPOperator>(
            AMP::Operator::OperatorBuilder::createOperator(
                meshAdapter, "LinearOperator", input_db, elementModel ) );

    // ---------------------------------------------------------------------------------------
    // create a mass linear BVP operator
    AMP::shared_ptr<AMP::Operator::ElementPhysicsModel> massElementModel;
    AMP::shared_ptr<AMP::Operator::LinearBVPOperator> massOperator =
        AMP::dynamic_pointer_cast<AMP::Operator::LinearBVPOperator>(
            AMP::Operator::OperatorBuilder::createOperator(
                meshAdapter, "MassLinearOperator", input_db, massElementModel ) );

    // ---------------------------------------------------------------------------------------
    //  create neutronics source
    AMP_INSIST( input_db->keyExists( "NeutronicsOperator" ),
                "Key ''NeutronicsOperator'' is missing!" );
    AMP::shared_ptr<AMP::Database> neutronicsOp_db = input_db->getDatabase( "NeutronicsOperator" );
    AMP::shared_ptr<AMP::Operator::NeutronicsRhsParameters> neutronicsParams(
        new AMP::Operator::NeutronicsRhsParameters( neutronicsOp_db ) );
    AMP::shared_ptr<AMP::Operator::NeutronicsRhs> neutronicsOperator(
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

    AMP::shared_ptr<AMP::Operator::ElementPhysicsModel> sourceTransportModel;
    AMP::shared_ptr<AMP::Operator::VolumeIntegralOperator> sourceOperator =
        AMP::dynamic_pointer_cast<AMP::Operator::VolumeIntegralOperator>(
            AMP::Operator::OperatorBuilder::createOperator(
                meshAdapter, "VolumeIntegralOperator", input_db, sourceTransportModel ) );

    // Create the power (heat source) vector.
    AMP::LinearAlgebra::Variable::shared_ptr powerInWattsVar = sourceOperator->getOutputVariable();
    AMP::LinearAlgebra::Vector::shared_ptr powerInWattsVec =
        AMP::LinearAlgebra::createVector( nodalDofMap, powerInWattsVar );
    powerInWattsVec->zero();

    // convert the vector of specific power to power for a given basis.
    sourceOperator->apply( SpecificPowerVec, powerInWattsVec );

    // ---------------------------------------------------------------------------------------
    // create vectors for initial conditions (IC) and time derivative at IC

    AMP::LinearAlgebra::Vector::shared_ptr initialCondition =
        AMP::LinearAlgebra::createVector( nodalDofMap, outputVar );
    AMP::LinearAlgebra::Vector::shared_ptr initialConditionPrime =
        AMP::LinearAlgebra::createVector( nodalDofMap, outputVar );
    AMP::LinearAlgebra::Vector::shared_ptr f =
        AMP::LinearAlgebra::createVector( nodalDofMap, outputVar );

    //----------------------------------------------------------------------------------------------------------------------------------------------//
    // set initial conditions, initialize created vectors

    AMP::Mesh::MeshIterator node     = meshAdapter->getIterator( AMP::Mesh::GeomType::Vertex, 0 );
    AMP::Mesh::MeshIterator end_node = node.end();

    int counter = 0;
    for ( ; node != end_node; ++node ) {
        counter++;

        std::vector<size_t> bndGlobalIds;
        nodalDofMap->getDOFs( node->globalID(), bndGlobalIds );

        auto pt   = node->coord();
        double px = pt[0];
        double py = pt[1];
        double pz = pt[2];

        double val = __INIT_FN__( px, py, pz, 0 );

        for ( auto &bndGlobalId : bndGlobalIds ) {
            initialCondition->setValueByGlobalID( bndGlobalId, val );
            // ** please do not set the time derivative to be non-zero!!
            // ** as this causes trouble with the boundary - BP, 07/16/2010
            initialConditionPrime->setValueByGlobalID( bndGlobalId, 0.0 );

        } // end for i
    }     // end for node
    initialCondition->makeConsistent( AMP::LinearAlgebra::Vector::ScatterType::CONSISTENT_SET );
    initialConditionPrime->makeConsistent(
        AMP::LinearAlgebra::Vector::ScatterType::CONSISTENT_SET );

    std::cout << "With Counter " << counter << " Max initial temp " << initialCondition->max()
              << " Min initial temp " << initialCondition->min() << std::endl;

    // create a copy of the rhs which can be modified at each time step (maybe)
    f->copyVector( powerInWattsVec );
    // modify the rhs to take into account boundary conditions
    nonlinearOperator->modifyRHSvector( f );
    nonlinearOperator->modifyInitialSolutionVector( initialCondition );

    // ---------------------------------------------------------------------------------------
    // create a linear time operator
    AMP::shared_ptr<AMP::Database> timeOperator_db( new AMP::Database( "TimeOperatorDatabase" ) );
    timeOperator_db->putScalar<double>( "CurrentDt", 0.01 );
    timeOperator_db->putScalar<std::string>( "name", "TimeOperator" );
    timeOperator_db->putScalar<bool>( "bLinearMassOperator", true );
    timeOperator_db->putScalar<bool>( "bLinearRhsOperator", false );
    timeOperator_db->putScalar<double>( "ScalingFactor", 1.0 / 0.01 );
    timeOperator_db->putScalar<double>( "CurrentTime", .0 );

    AMP::shared_ptr<AMP::TimeIntegrator::TimeOperatorParameters> timeOperatorParameters(
        new AMP::TimeIntegrator::TimeOperatorParameters( timeOperator_db ) );
    timeOperatorParameters->d_pRhsOperator  = linearOperator;
    timeOperatorParameters->d_pMassOperator = massOperator;
    // timeOperatorParameters->d_pMassOperator = massLinearOperator;
    AMP::shared_ptr<AMP::Operator::Operator> linearTimeOperator(
        new AMP::TimeIntegrator::LinearTimeOperator( timeOperatorParameters ) );

    AMP::LinearAlgebra::Vector::shared_ptr residualVec =
        AMP::LinearAlgebra::createVector( nodalDofMap, outputVar );

    linearOperator->apply( initialCondition, residualVec );
    std::cout << "Residual Norm of linearTimeOp apply : " << residualVec->L2Norm() << std::endl;

    massOperator->apply( initialCondition, residualVec );
    std::cout << "Residual Norm of linearTimeOp apply : " << residualVec->L2Norm() << std::endl;

    linearTimeOperator->apply( initialCondition, residualVec );
    std::cout << "Residual Norm of linearTimeOp apply : " << residualVec->L2Norm() << std::endl;

    // ---------------------------------------------------------------------------------------
    // create a preconditioner

    // get the ida database
    AMP_INSIST( input_db->keyExists( "IDATimeIntegrator" ),
                "Key ''IDATimeIntegrator'' is missing!" );
    AMP::shared_ptr<AMP::Database> ida_db      = input_db->getDatabase( "IDATimeIntegrator" );
    AMP::shared_ptr<AMP::Database> pcSolver_db = ida_db->getDatabase( "Preconditioner" );
    AMP::shared_ptr<AMP::Solver::SolverStrategyParameters> pcSolverParams(
        new AMP::Solver::SolverStrategyParameters( pcSolver_db ) );
    pcSolverParams->d_pOperator = linearTimeOperator;

    if ( pcSolverParams.get() == nullptr ) {
        ut->failure( "Testing SolverStrategyParameters's constructor: FAIL" );
    } else {
        ut->passes( "Testing SolverStrategyParameters's constructor: PASS" );
    }

    AMP::shared_ptr<AMP::Solver::TrilinosMLSolver> pcSolver(
        new AMP::Solver::TrilinosMLSolver( pcSolverParams ) );

    if ( pcSolver.get() == nullptr ) {
        ut->failure( "Testing TrilinosMLSolver's constructor: FAIL" );
    } else {
        ut->passes( "Testing TrilinosMLSolver's constructor: PASS" );
    }

#ifdef USE_EXT_SILO
    AMP::Utilities::Writer::shared_ptr siloWriter = AMP::Utilities::Writer::buildWriter( "Silo" );
    siloWriter->registerMesh( meshAdapter );

    siloWriter->registerVector(
        initialCondition, meshAdapter, AMP::Mesh::GeomType::Vertex, "InitialSolution" );

    siloWriter->writeFile( input_file, 0 );
#endif

    // ---------------------------------------------------------------------------------------
    // create the IDA time integrator
    AMP::shared_ptr<AMP::TimeIntegrator::IDATimeIntegratorParameters> time_Params(
        new AMP::TimeIntegrator::IDATimeIntegratorParameters( ida_db ) );

    if ( ( time_Params.get() ) == nullptr ) {
        ut->failure( "Testing IDATimeIntegratorParameters' Constructor" );
    } else {
        ut->passes( "Testing IDATimeIntegratorParameters' Constructor" );
    }

    time_Params->d_pMassOperator = massOperator;
    // time_Params->d_pMassOperator = massLinearOperator;
    time_Params->d_operator        = nonlinearOperator;
    time_Params->d_pPreconditioner = pcSolver;

    time_Params->d_ic_vector       = initialCondition;
    time_Params->d_ic_vector_prime = initialConditionPrime;

    time_Params->d_pSourceTerm = f;
    time_Params->d_object_name = "IDATimeIntegratorParameters";

    std::cout << "Before IDATimeIntegrator" << std::endl;
    AMP::shared_ptr<AMP::TimeIntegrator::IDATimeIntegrator> pIDATimeIntegrator(
        new AMP::TimeIntegrator::IDATimeIntegrator( time_Params ) );

    if ( pIDATimeIntegrator.get() == nullptr ) {
        ut->failure( "Testing IDATimeIntegrator's constructor" );
    } else {
        ut->passes( "Tested IDATimeIntegrator's constructor" );
    }

    // ---------------------------------------------------------------------------------------
    // step in time
    int retval          = 0;
    double current_time = 0;
    double max          = 0;
    double min          = 0;
    int j               = 1;
    while ( pIDATimeIntegrator->getCurrentTime() < pIDATimeIntegrator->getFinalTime() ) {
        retval = pIDATimeIntegrator->advanceSolution( pIDATimeIntegrator->getCurrentDt(), false );
        // pIDATimeIntegrator->updateSolution();
        current_time = pIDATimeIntegrator->getCurrentTime();

        std::cout << j++ << "-th timestep" << std::endl;
        if ( retval == 0 ) {
            ut->passes( "Testing IDATimeIntegrator's advanceSolution. PASS!!" );
        } else {
            ut->failure( "Tested IDATimeIntegrator's advanceSolution. FAIL!!" );
        }

        max = pIDATimeIntegrator->getCurrentSolution()->max();
        min = pIDATimeIntegrator->getCurrentSolution()->min();

        std::cout << "current_time = " << current_time << std::endl;
        std::cout << "max val of the current solution = " << max << std::endl;
        std::cout << "min val of the current solution = " << min << std::endl;
    }


#ifdef USE_EXT_SILO

    AMP::LinearAlgebra::Vector::shared_ptr pSolution = pIDATimeIntegrator->getCurrentSolution();
    siloWriter->registerVector( pSolution, meshAdapter, AMP::Mesh::GeomType::Vertex, "Solution" );

    siloWriter->writeFile( input_file, 1 );
#endif

    if ( ut->NumFailLocal() == 0 ) {
        ut->passes( "testIDATimeIntegrator successful" );
    }
}


//---------------------------------------------------------------------------//

int testIDA_NonlinearBVPOperator( int argc, char *argv[] )
{
    AMP::AMPManager::startup( argc, argv );
    AMP::UnitTest ut;

    IDATimeIntegratorTest( &ut );

    ut.report();

    int num_failed = ut.NumFailGlobal();
    AMP::AMPManager::shutdown();
    return num_failed;
}

//---------------------------------------------------------------------------//
//                        end of SundialsVectorTest.cc
//---------------------------------------------------------------------------//
