#include "AMP/materials/Material.h"
#include "AMP/operators/NeutronicsRhs.h"
#include "AMP/utils/AMPManager.h"
#include "AMP/utils/AMP_MPI.h"
#include "AMP/utils/Database.h"
#include "AMP/utils/InputDatabase.h"
#include "AMP/utils/InputManager.h"
#include "AMP/utils/PIO.h"
#include "AMP/utils/Utilities.h"
#include "AMP/utils/shared_ptr.h"
#include "AMP/vectors/Variable.h"
#include <string>

#include "AMP/operators/NeutronicsRhs.h"
#include "AMP/operators/NonlinearBVPOperator.h"
#include "AMP/operators/libmesh/VolumeIntegralOperator.h"
#include "AMP/utils/Writer.h"
#include "AMP/vectors/Vector.h"
/* libMesh files */
#include "dense_matrix.h"
#include "dof_map.h"
#include "elem.h"
#include "equation_systems.h"
#include "fe.h"
#include "libmesh.h"
#include "linear_implicit_system.h"
#include "mesh.h"
#include "mesh_generation.h"
#include "petsc_matrix.h"
#include "petsc_vector.h"
#include "quadrature_gauss.h"
#include "sparse_matrix.h"

#include "AMP/operators/diffusion/DiffusionLinearFEOperator.h"
#include "AMP/operators/libmesh/MassLinearFEOperator.h"

#include "AMP/operators/DirichletMatrixCorrection.h"
#include "AMP/operators/DirichletVectorCorrection.h"
#include "AMP/operators/IsotropicElasticModel.h"
#include "AMP/operators/LinearBVPOperator.h"
#include "AMP/operators/MassMatrixCorrection.h"
#include "AMP/operators/MechanicsLinearElement.h"
#include "AMP/operators/MechanicsLinearFEOperator.h"

#include "AMP/solvers/trilinos/ml/TrilinosMLSolver.h"
#include "AMP/time_integrators/sundials/IDATimeIntegrator.h"
#include "AMP/time_integrators/sundials/IDATimeOperator.h"

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
#define __INIT_FN__( x, y, z )                                            \
    ( 1 - sqrt( x * x + y * y + z * z ) * sqrt( x * x + y * y + z * z ) * \
              sqrt( x * x + y * y + z * z ) )
#define __INIT_PRIME_FN__( x, y, z )                                  \
    ( sqrt( x * x + y * y + z * z ) * sqrt( x * x + y * y + z * z ) * \
      sqrt( x * x + y * y + z * z ) )

void IDATimeIntegratorTest( AMP::UnitTest *ut )
{
    std::string input_file = "input_testIDA-NonlinearOperator-long";
    std::string log_file   = "output_testIDA-NonlinearOperator-long";

    AMP::PIO::logOnlyNodeZero( log_file );

    AMP::shared_ptr<AMP::InputDatabase> input_db( new AMP::InputDatabase( "input_db" ) );
    AMP::AMP_MPI globalComm( AMP_COMM_WORLD );
    AMP::InputManager::getManager()->parseInputFile( input_file, input_db );
    input_db->printClassData( AMP::plog );

    AMP::Mesh::MeshManagerParameters::shared_ptr meshmgrParams(
        new AMP::Mesh::MeshManagerParameters( input_db ) );
    AMP::Mesh::MeshManager::shared_ptr manager( new AMP::Mesh::MeshManager( meshmgrParams ) );
    AMP::Mesh::MeshManager::Adapter::shared_ptr meshAdapter = manager->getMesh( "ida" );

    //----------------------------------------------------------------------------------------------------------------------------------------------//
    // create a nonlinear BVP operator for nonlinear BVP operator
    AMP_INSIST( input_db->keyExists( "NonlinearOperator" ), "key missing!" );

    AMP::shared_ptr<AMP::Operator::ElementPhysicsModel> elementModel;
    AMP::shared_ptr<AMP::Operator::DiffusionNonlinearFEOperator> nonlinearOperator;
    AMP::shared_ptr<AMP::Database> nonlinearDatabase = input_db->getDatabase( "NonlinearOperator" );
    AMP::shared_ptr<AMP::Operator::Operator> genericOperator =
        AMP::Operator::OperatorBuilder::createOperator(
            meshAdapter, nonlinearDatabase, elementModel );
    nonlinearOperator =
        AMP::dynamic_pointer_cast<AMP::Operator::DiffusionNonlinearFEOperator>( genericOperator );

    // ---------------------------------------------------------------------------------------
    // create a linear BVP operator
    AMP::shared_ptr<AMP::Operator::DiffusionLinearFEOperator> linearOperator;
    AMP::shared_ptr<AMP::InputDatabase> linearOp_db =
        AMP::dynamic_pointer_cast<AMP::InputDatabase>( input_db->getDatabase( "LinearOperator" ) );
    linearOperator = AMP::dynamic_pointer_cast<AMP::Operator::DiffusionLinearFEOperator>(
        AMP::Operator::OperatorBuilder::createOperator( meshAdapter, linearOp_db, elementModel ) );


    // ---------------------------------------------------------------------------------------
    // create a mass linear BVP operator
    AMP::shared_ptr<AMP::Operator::ElementPhysicsModel> massElementModel;
    AMP::shared_ptr<AMP::Operator::MassLinearFEOperator> massOperator;
    AMP::shared_ptr<AMP::InputDatabase> massOp_db = AMP::dynamic_pointer_cast<AMP::InputDatabase>(
        input_db->getDatabase( "MassLinearOperator" ) );
    massOperator = AMP::dynamic_pointer_cast<AMP::Operator::MassLinearFEOperator>(
        AMP::Operator::OperatorBuilder::createOperator(
            meshAdapter, massOp_db, massElementModel ) );

    // ---------------------------------------------------------------------------------------
    // create Source operators
    AMP_INSIST( input_db->keyExists( "VolumeIntegralOperator1" ), "key missing!" );
    AMP::shared_ptr<AMP::Operator::ElementPhysicsModel> manufacturedSourceModel1;
    AMP::shared_ptr<AMP::Database> sourceDatabase1 =
        input_db->getDatabase( "VolumeIntegralOperator1" );
    AMP::shared_ptr<AMP::Operator::VolumeIntegralOperator> sourceOperator1 =
        AMP::dynamic_pointer_cast<AMP::Operator::VolumeIntegralOperator>(
            AMP::Operator::OperatorBuilder::createOperator(
                meshAdapter, sourceDatabase1, manufacturedSourceModel1 ) );

    AMP_INSIST( input_db->keyExists( "VolumeIntegralOperator2" ), "key missing!" );
    AMP::shared_ptr<AMP::Operator::ElementPhysicsModel> manufacturedSourceModel2;
    AMP::shared_ptr<AMP::Database> sourceDatabase2 =
        input_db->getDatabase( "VolumeIntegralOperator2" );
    AMP::shared_ptr<AMP::Operator::VolumeIntegralOperator> sourceOperator2 =
        AMP::dynamic_pointer_cast<AMP::Operator::VolumeIntegralOperator>(
            AMP::Operator::OperatorBuilder::createOperator(
                meshAdapter, sourceDatabase2, manufacturedSourceModel2 ) );


    // create the following shared pointers for ease of use
    AMP::LinearAlgebra::Vector::shared_ptr nullVec;


    // ---------------------------------------------------------------------------------------
    // create vectors for initial conditions (IC) and time derivative at IC
    AMP::LinearAlgebra::Variable::shared_ptr inputVar =
        nonlinearOperator->getInputVariable( AMP::Operator::Diffusion::TEMPERATURE );
    AMP::LinearAlgebra::Variable::shared_ptr outputVar = nonlinearOperator->getOutputVariable();

    AMP::LinearAlgebra::Variable::shared_ptr inputVar_mass = massOperator->getInputVariable();
    if ( *inputVar == *inputVar_mass ) {
        cout << "IDARhsOperator and massOpeator have the same input variable" << endl;
    }

    AMP::LinearAlgebra::Variable::shared_ptr inputVarSource1 =
        sourceOperator1->getInputVariable( AMP::Operator::Diffusion::TEMPERATURE );
    cout << "IDARhsOp, sourceOp1 ? " << ( *inputVar == *inputVarSource1 ) << endl;

    AMP::LinearAlgebra::Vector::shared_ptr initialCondition = meshAdapter->createVector( inputVar );
    AMP::LinearAlgebra::Vector::shared_ptr initialConditionPrime =
        meshAdapter->createVector( inputVar );
    AMP::LinearAlgebra::Vector::shared_ptr scratchVec = meshAdapter->createVector( inputVar );
    AMP::LinearAlgebra::Vector::shared_ptr f          = meshAdapter->createVector( outputVar );

    cout << "*(initialCondition->getVariable()) == *inputVarSource1 = "
         << ( *( initialCondition->getVariable() ) == *inputVarSource1 ) << endl;
    cout << "initialCondition->getVariable() = " << initialCondition->getVariable() << endl;
    cout << "sourceOperator1->getInputVariable(AMP::Operator::Diffusion::TEMPERATURE) = "
         << sourceOperator1->getInputVariable( AMP::Operator::Diffusion::TEMPERATURE ) << endl;
    //----------------------------------------------------------------------------------------------------------------------------------------------//
    // set initial conditions, initialize created vectors

    AMP::Mesh::MeshManager::Adapter::NodeIterator node     = meshAdapter->beginNode();
    AMP::Mesh::MeshManager::Adapter::NodeIterator end_node = meshAdapter->endNode();
    AMP::Mesh::DOFMap::shared_ptr dof_map                  = meshAdapter->getDOFMap( inputVar );

    int counter = 0;
    for ( ; node != end_node; ++node ) {
        counter++;

        std::vector<unsigned int> bndGlobalIds;
        std::vector<unsigned int> d_dofIds;
        d_dofIds.resize( 0 );
        dof_map->getDOFs( *node, bndGlobalIds, d_dofIds );

        double px = ( *node ).x();
        double py = ( *node ).y();
        double pz = ( *node ).z();

        double r         = sqrt( px * px + py * py + pz * pz );
        double val       = __INIT_FN__( px, py, pz );
        double val_prime = __INIT_PRIME_FN__( px, py, pz );
        // cout << "val = " << val << endl;
        cout << "r = " << r << endl;
        for ( unsigned int i = 0; i < bndGlobalIds.size(); i++ ) {
            initialCondition->setValueByGlobalID( bndGlobalIds[i], val );
            // ** please do not set the time derivative to be non-zero!!
            // ** as this causes trouble with the boundary - BP, 07/16/2010
            initialConditionPrime->setValueByGlobalID( bndGlobalIds[i], val_prime );
            scratchVec->setValueByGlobalID( bndGlobalIds[i], r );
        } // end for i
    }     // end for node

    cout << "initialCondition->max() = " << initialCondition->max() << endl;
    cout << "initialCondition->min() = " << initialCondition->min() << endl;
    cout << "scratchVec->max() = " << scratchVec->max() << endl;
    cout << "scratchVec->min() = " << scratchVec->min() << endl;
    // ---------------------------------------------------------------------------------------
    // create a linear time operator
    AMP::shared_ptr<AMP::InputDatabase> timeOperator_db(
        new AMP::InputDatabase( "TimeOperatorDatabase" ) );
    timeOperator_db->putDouble( "CurrentDt", 0.01 );
    timeOperator_db->putString( "name", "TimeOperator" );
    timeOperator_db->putBool( "bLinearMassOperator", true );
    timeOperator_db->putBool( "bLinearRhsOperator", false );
    timeOperator_db->putDouble( "ScalingFactor", 1.0 / 0.01 );

    AMP::shared_ptr<AMP::TimeIntegrator::TimeOperatorParameters> timeOperatorParameters(
        new AMP::TimeIntegrator::TimeOperatorParameters( timeOperator_db ) );
    timeOperatorParameters->d_pRhsOperator  = linearOperator;
    timeOperatorParameters->d_pMassOperator = massOperator;
    timeOperatorParameters->d_pSourceOperators.resize( 2 );
    timeOperatorParameters->d_pSourceOperators[0] = sourceOperator1;
    timeOperatorParameters->d_pSourceOperators[1] = sourceOperator2;

    // timeOperatorParameters->d_pMassOperator = massLinearOperator;
    AMP::shared_ptr<AMP::Operator::Operator> linearTimeOperator(
        new AMP::TimeIntegrator::LinearTimeOperator( timeOperatorParameters ) );
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

    if ( pcSolverParams.get() == NULL ) {
        ut.failure( "Testing SolverStrategyParameters's constructor: FAIL" );
    } else {
        ut.passes( "Testing SolverStrategyParameters's constructor: PASS" );
    }

    AMP::shared_ptr<AMP::Solver::TrilinosMLSolver> pcSolver(
        new AMP::Solver::TrilinosMLSolver( pcSolverParams ) );

    if ( pcSolver.get() == NULL ) {
        ut.failure( "Testing TrilinosMLSolver's constructor: FAIL" );
    } else {
        ut.passes( "Testing TrilinosMLSolver's constructor: PASS" );
    }

    // ---------------------------------------------------------------------------------------
    // create the IDA time integrator
    AMP::shared_ptr<AMP::TimeIntegrator::IDATimeIntegratorParameters> time_Params(
        new AMP::TimeIntegrator::IDATimeIntegratorParameters( ida_db ) );

    if ( ( time_Params.get() ) == NULL ) {
        ut.failure( "Testing IDATimeIntegratorParameters' Constructor" );
    } else {
        ut.passes( "Testing IDATimeIntegratorParameters' Constructor" );
    }

    time_Params->d_pMassOperator = massOperator;
    // time_Params->d_pMassOperator = massLinearOperator;
    time_Params->d_operator        = nonlinearOperator;
    time_Params->d_pPreconditioner = pcSolver;

    // JL
    time_Params->d_pSourceOperators.resize( 2 );
    time_Params->d_pSourceOperators[0] = sourceOperator1;
    time_Params->d_pSourceOperators[1] = sourceOperator2;


    time_Params->d_ic_vector       = initialCondition;
    time_Params->d_ic_vector_prime = initialConditionPrime;

    time_Params->d_pSourceTerm = f;
    time_Params->d_object_name = "IDATimeIntegratorParameters";

    cout << "Before IDATimeIntegrator" << endl;
    AMP::shared_ptr<AMP::TimeIntegrator::IDATimeIntegrator> pIDATimeIntegrator(
        new AMP::TimeIntegrator::IDATimeIntegrator( time_Params ) );

    if ( pIDATimeIntegrator.get() == NULL ) {
        ut.failure( "Testing IDATimeIntegrator's constructor" );
    } else {
        ut.passes( "Tested IDATimeIntegrator's constructor" );
    }

    // ---------------------------------------------------------------------------------------
    // step in time
    int retval          = 0;
    double current_time = 0;
    double max          = 0;
    double min          = 0;
    int j               = 1;
    while ( pIDATimeIntegrator->getCurrentTime() < pIDATimeIntegrator->getFinalTime() ) {
        retval = pIDATimeIntegrator->advanceSolution( pIDATimeIntegrator->getCurrentDt(), 0 );
        // pIDATimeIntegrator->updateSolution();
        current_time = pIDATimeIntegrator->getCurrentTime();

        cout << j++ << "-th timestep" << endl;
        if ( retval == 0 ) {
            ut.passes( "Testing IDATimeIntegrator's advanceSolution. PASS!!" );
        } else {
            ut.failure( "Tested IDATimeIntegrator's advanceSolution. FAIL!!" );
        }

        max = pIDATimeIntegrator->getCurrentSolution()->max();
        min = pIDATimeIntegrator->getCurrentSolution()->min();

        cout << "current_time = " << current_time << endl;
        cout << "max val of the current solution = " << max << endl;
        cout << "min val of the current solution = " << min << endl;

        AMP::Mesh::MeshManager::Adapter::OwnedNodeIterator iterator = meshAdapter->beginOwnedNode();
        double x1, y1, z1, r1, sol, cal;

        for ( ; iterator != meshAdapter->endOwnedNode(); iterator++ ) {
            cal = pIDATimeIntegrator->getCurrentSolution()->getValueByGlobalID(
                dof_map->getGlobalID( iterator->globalID(), 0 ) );

            x1 = iterator->x();
            y1 = iterator->y();
            z1 = iterator->z();
            r1 = sqrt( x1 * x1 + y1 * y1 + z1 * z1 );

            sol = 1 - r1 * r1 * r1 * exp( -current_time );
            // cout << "fabs(cal-sol)/fabs(sol) = " << fabs(cal - sol)/fabs(sol) << endl;

            scratchVec->setValueByGlobalID( dof_map->getGlobalID( iterator->globalID(), 0 ),
                                            fabs( cal - sol ) / fabs( sol ) );

            if ( fabs( cal - sol ) / fabs( sol ) > 0.2 ) {
                cout << "x1 = " << x1 << " y1 = " << y1 << " z1 = " << z1 << endl;
            }

            /*
            if( fabs(cal - sol) > cal*1e-3 ) {
                passes = 0;
                ITFAILS;
            }
             */
        }
        cout << "max(rel_error) = " << scratchVec->max() << endl;
        cout << "min(rel_error) = " << scratchVec->min() << endl;
    }


#ifdef USE_EXT_SILO
    AMP::LinearAlgebra::Vector::shared_ptr pSolution = pIDATimeIntegrator->getCurrentSolution();
    meshAdapter->registerVectorAsData( pSolution, "ThemalSolutionVector" );
    manager->writeFile<AMP::SiloIO>( "IDA-NonlinearBVP", 1 );
#endif
    AMP::AMPManager::shutdown();

    if ( ut.numFails == 0 ) {
        ut.passes( "testIDATimeIntegrator successful" );
    }
}


//---------------------------------------------------------------------------//

int main( int argc, char *argv[] )
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
