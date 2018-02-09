#include "AMP/materials/Material.h"
#include "AMP/operators/NeutronicsRhs.h"
#include "AMP/utils/AMPManager.h"
#include "AMP/utils/AMP_MPI.h"
#include "AMP/utils/Database.h"
#include "AMP/utils/InputDatabase.h"
#include "AMP/utils/InputManager.h"
#include "AMP/utils/PIO.h"
#include "AMP/utils/UnitTest.h"
#include "AMP/utils/Utilities.h"
#include "AMP/utils/shared_ptr.h"
#include "AMP/vectors/Variable.h"
#include <string>

#include "AMP/operators/NeutronicsRhs.h"
#include "AMP/operators/diffusion/DiffusionLinearElement.h"
#include "AMP/operators/diffusion/DiffusionTransportModel.h"
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
#include "AMP/operators/libmesh/MassLinearElement.h"
#include "AMP/operators/libmesh/MassLinearFEOperator.h"

#include "AMP/operators/DirichletMatrixCorrection.h"
#include "AMP/operators/DirichletVectorCorrection.h"
#include "AMP/operators/IsotropicElasticModel.h"
#include "AMP/operators/LinearBVPOperator.h"
#include "AMP/operators/MechanicsLinearElement.h"
#include "AMP/operators/MechanicsLinearFEOperator.h"
#include "AMP/operators/NeumannVectorCorrection.h"

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
#define __INIT_FN__( x, y, z, t )                                           \
    ( 300 + exp( -0.015 * __PI__ * __PI__ * t ) * cos( 0.1 * __PI__ * x ) * \
                cos( 0.1 * __PI__ * y ) * cos( 0.05 * __PI__ * z ) )

void IDATimeIntegratorTest( AMP::UnitTest *ut )
{
    std::string input_file = "input-LinearIDA-1";
    std::string log_file   = "output-LinearIDA-1";

    AMP::AMP_MPI::initialize();
    AMP::AMPManager::startup();

    AMP::PIO::logOnlyNodeZero( log_file );

    AMP::shared_ptr<AMP::InputDatabase> input_db( new AMP::InputDatabase( "input_db" ) );
    AMP::InputManager::getManager()->parseInputFile( input_file, input_db );
    input_db->printClassData( AMP::plog );


    AMP::Mesh::MeshManagerParameters::shared_ptr mgrParams(
        new AMP::Mesh::MeshManagerParameters( input_db ) );
    AMP::Mesh::MeshManager manager( mgrParams );

    AMP::Mesh::MeshManager::Adapter::shared_ptr meshAdapter = manager.getMesh( "ida" );

    AMP_INSIST( input_db->keyExists( "DiffusionLinearElement" ),
                "Key ''DiffusionLinearElement'' is missing!" );
    AMP::shared_ptr<AMP::Database> elemOp_db = input_db->getDatabase( "DiffusionLinearElement" );
    AMP::shared_ptr<AMP::Operator::ElementOperationParameters> elemOpParams(
        new AMP::Operator::ElementOperationParameters( elemOp_db ) );
    AMP::shared_ptr<AMP::Operator::DiffusionLinearElement> thermLinElem(
        new AMP::Operator::DiffusionLinearElement( elemOpParams ) );

    AMP_INSIST( input_db->keyExists( "DiffusionTransportModel" ),
                "Key ''DiffusionTransportModel'' is missing!" );
    AMP::shared_ptr<AMP::Database> thermModel_db =
        input_db->getDatabase( "DiffusionTransportModel" );
    AMP::shared_ptr<AMP::Operator::DiffusionTransportModelParameters> thermModelParams(
        new AMP::Operator::DiffusionTransportModelParameters( thermModel_db ) );
    cout << "test 1" << endl;
    AMP::shared_ptr<AMP::Operator::DiffusionTransportModel> thermalDiffnModel(
        new AMP::Operator::DiffusionTransportModel( thermModelParams ) );


    // Creating the thermal operator. BVP.
    AMP::shared_ptr<AMP::Operator::ElementPhysicsModel> transportModel0;
    AMP::shared_ptr<AMP::InputDatabase> bvpDatabase =
        AMP::dynamic_pointer_cast<AMP::InputDatabase>( input_db->getDatabase( "LinearOperator" ) );
    AMP::shared_ptr<AMP::Operator::LinearBVPOperator> diffusionOperator =
        AMP::dynamic_pointer_cast<AMP::Operator::LinearBVPOperator>(
            AMP::Operator::OperatorBuilder::createOperator(
                meshAdapter, "LinearOperator", input_db, transportModel0 ) );

    AMP::LinearAlgebra::Vector::shared_ptr TemperatureInKelvinVec =
        meshAdapter->createVector( diffusionOperator->getInputVariable() );
    AMP::LinearAlgebra::Vector::shared_ptr RightHandSideVec =
        meshAdapter->createVector( diffusionOperator->getOutputVariable() );
    AMP::LinearAlgebra::Vector::shared_ptr ResidualVec =
        meshAdapter->createVector( diffusionOperator->getOutputVariable() );


    AMP_INSIST( input_db->keyExists( "DiffusionLinearFEOperator" ),
                "Key ''DiffusionLinearFEOperator'' is missing!" );
    AMP::shared_ptr<AMP::Database> thermAssembly_db =
        input_db->getDatabase( "DiffusionLinearFEOperator" );
    AMP::shared_ptr<AMP::Operator::DiffusionLinearFEOperatorParameters> thermOpParams(
        new AMP::Operator::DiffusionLinearFEOperatorParameters( thermAssembly_db ) );

    thermOpParams->d_elemOp         = thermLinElem;
    thermOpParams->d_transportModel = thermalDiffnModel;
    thermOpParams->d_MeshAdapter    = manager.getMesh( "ida" );
    // thermOpParams->d_MeshAdapter     = meshAdapter;
    AMP::shared_ptr<AMP::Operator::DiffusionLinearFEOperator> thermOp(
        new AMP::Operator::DiffusionLinearFEOperator( thermOpParams ) );
    AMP::shared_ptr<AMP::Operator::DiffusionLinearFEOperator> thermOp2(
        new AMP::Operator::DiffusionLinearFEOperator( thermOpParams ) );


    AMP_INSIST( input_db->keyExists( "DirichletVectorCorrection" ),
                "Key ''DirichletVectorCorrection'' is missing!" );
    AMP::shared_ptr<AMP::Database> temp1_db = input_db->getDatabase( "DirichletVectorCorrection" );
    AMP::shared_ptr<AMP::Operator::DirichletVectorCorrectionParameters> temp1_dirichletVecParams(
        new AMP::Operator::DirichletVectorCorrectionParameters( temp1_db ) );

    /// put source with mass matrix
    AMP_INSIST( input_db->keyExists( "MassLinearElement" ),
                "Key ''MassLinearElement'' is missing!" );
    AMP::shared_ptr<AMP::Database> srcmass_db = input_db->getDatabase( "MassLinearElement" );
    AMP::shared_ptr<AMP::Operator::ElementOperationParameters> srcmassParams(
        new AMP::Operator::ElementOperationParameters( srcmass_db ) );
    AMP::shared_ptr<AMP::MassLinearElement> massLinElem(
        new AMP::MassLinearElement( srcmassParams ) );

    AMP_INSIST( input_db->keyExists( "MassDensityModel" ), "Key ''MassDensityModel'' is missing!" );
    AMP::shared_ptr<AMP::Database> massModel_db = input_db->getDatabase( "MassDensityModel" );
    AMP::shared_ptr<AMP::Operator::MassDensityModelParameters> massModelParams(
        new AMP::Operator::MassDensityModelParameters( massModel_db ) );
    cout << "test 2 " << endl;
    AMP::shared_ptr<AMP::Operator::MassDensityModel> massDensityModel(
        new AMP::Operator::MassDensityModel( massModelParams ) );

    AMP_INSIST( input_db->keyExists( "MassLinearFEOperator" ),
                "Key ''MassLinearFEOperator'' is missing!" );
    AMP::shared_ptr<AMP::Database> srcAssembly_db = input_db->getDatabase( "MassLinearFEOperator" );
    AMP::shared_ptr<AMP::Operator::MassLinearFEOperatorParameters> srcOpParams(
        new AMP::Operator::MassLinearFEOperatorParameters( srcAssembly_db ) );


    // get the ida database
    AMP_INSIST( input_db->keyExists( "IDATimeIntegrator" ),
                "Key ''IDATimeIntegrator'' is missing!" );
    AMP::shared_ptr<AMP::Database> ida_db      = input_db->getDatabase( "IDATimeIntegrator" );
    AMP::shared_ptr<AMP::Database> pcSolver_db = ida_db->getDatabase( "Preconditioner" );
    AMP::shared_ptr<AMP::Solver::SolverStrategyParameters> pcSolverParams(
        new AMP::Solver::SolverStrategyParameters( pcSolver_db ) );


    // create the following shared pointers for ease of use
    AMP::LinearAlgebra::Vector::shared_ptr nullVec;
    ////////////////////////////////////
    //  CREATE THE NEUTRONICS SOURCE  //
    ////////////////////////////////////
    AMP_INSIST( input_db->keyExists( "NeutronicsOperator" ),
                "Key ''NeutronicsOperator'' is missing!" );
    AMP::shared_ptr<AMP::Database> neutronicsOp_db = input_db->getDatabase( "NeutronicsOperator" );
    AMP::shared_ptr<AMP::Operator::NeutronicsRhsParameters> neutronicsParams(
        new AMP::Operator::NeutronicsRhsParameters( neutronicsOp_db ) );
    neutronicsParams->d_MeshAdapter = meshAdapter;
    AMP::shared_ptr<AMP::Operator::NeutronicsRhs> neutronicsOperator(
        new AMP::Operator::NeutronicsRhs( neutronicsParams ) );

    AMP::LinearAlgebra::Variable::shared_ptr SpecificPowerVar =
        neutronicsOperator->getOutputVariable();
    AMP::LinearAlgebra::Vector::shared_ptr SpecificPowerVec =
        meshAdapter->createVector( SpecificPowerVar );

    neutronicsOperator->apply( nullVec, SpecificPowerVec );


    /////////////////////////////////////////////////////
    //  Integrate Nuclear Rhs over Desnity * GeomType::Volume //
    /////////////////////////////////////////////////////

    AMP_INSIST( input_db->keyExists( "VolumeIntegralOperator" ), "key missing!" );

    AMP::shared_ptr<AMP::Operator::ElementPhysicsModel> stransportModel;
    AMP::shared_ptr<AMP::Database> sourceDatabase =
        input_db->getDatabase( "VolumeIntegralOperator" );
    AMP::shared_ptr<AMP::Operator::VolumeIntegralOperator> sourceOperator =
        AMP::dynamic_pointer_cast<AMP::Operator::VolumeIntegralOperator>(
            AMP::Operator::OperatorBuilder::createOperator(
                meshAdapter, "VolumeIntegralOperator", input_db, stransportModel ) );

    // Create the power (heat source) vector.
    AMP::LinearAlgebra::Variable::shared_ptr PowerInWattsVar = sourceOperator->getOutputVariable();
    AMP::LinearAlgebra::Vector::shared_ptr PowerInWattsVec =
        meshAdapter->createVector( PowerInWattsVar );
    PowerInWattsVec->zero();

    // convert the vector of specific power to power for a given basis.
    sourceOperator->apply( SpecificPowerVec, PowerInWattsVec );

    double t1;
    t1 = SpecificPowerVec->L2Norm();
    cout << "n1 = " << t1 << std::endl;
    t1 = PowerInWattsVec->L2Norm();
    cout << "n1 = " << t1 << std::endl;


    RightHandSideVec->setToScalar( 0.0 );

    RightHandSideVec->copyVector( PowerInWattsVec );

    diffusionOperator->modifyRHSvector( RightHandSideVec );


    cout << "PowerInWattsVec->max() = " << PowerInWattsVec->max() << endl;
    cout << "PowerInWattsVec->min() = " << PowerInWattsVec->min() << endl;
    //
    // PowerInWattsVec->zero();
    // PowerInWattsVec->setToScalar(100.0);
    //----------------------------------------------------------------------------------------------------------------------------------------------//


    // ---------------------------------------------------------------------------------------
    // create a mass linear BVP operator
    AMP::shared_ptr<AMP::Operator::ElementPhysicsModel> massElementModel;
    AMP::shared_ptr<AMP::Operator::LinearBVPOperator> massOperator;
    AMP::shared_ptr<AMP::InputDatabase> massOp_db = AMP::dynamic_pointer_cast<AMP::InputDatabase>(
        input_db->getDatabase( "MassLinearOperator" ) );
    massOperator = AMP::dynamic_pointer_cast<AMP::Operator::LinearBVPOperator>(
        AMP::Operator::OperatorBuilder::createOperator(
            meshAdapter, "MassLinearOperator", input_db, massElementModel ) );


    // do we really need dynamic cast here?
    srcOpParams->d_elemOp =
        AMP::dynamic_pointer_cast<AMP::Operator::ElementOperation>( massLinElem );
    // srcOpParams->d_elemOp = massLinElem;
    srcOpParams->d_densityModel = massDensityModel;
    // srcOpParams->d_MeshAdapter     = meshAdapter; //manager.getMesh("ida");
    srcOpParams->d_MeshAdapter = manager.getMesh( "ida" );
    AMP::shared_ptr<AMP::Operator::MassLinearFEOperator> srcOp(
        new AMP::Operator::MassLinearFEOperator( srcOpParams ) );


    if ( pcSolverParams.get() == NULL ) {
        ut->failure( "Testing SolverStrategyParameters's constructor: FAIL" );
    } else {
        ut->passes( "Testing SolverStrategyParameters's constructor: PASS" );
    }

    // is this really necessary?

    // AMP::LinearAlgebra::Variable::shared_ptr ThermalInpVar = thermOp->getInputVariable();
    // AMP::LinearAlgebra::Variable::shared_ptr ThermalOutVar = thermOp->getOutputVariable();
    AMP::LinearAlgebra::Variable::shared_ptr ThermalInpVar = diffusionOperator->getInputVariable();
    AMP::LinearAlgebra::Variable::shared_ptr ThermalOutVar = diffusionOperator->getOutputVariable();


    AMP::LinearAlgebra::Vector::shared_ptr ThermalInpVec_ic =
        meshAdapter->createVector( ThermalInpVar );
    AMP::LinearAlgebra::Vector::shared_ptr ThermalInpVec_ic_prime =
        meshAdapter->createVector( ThermalInpVar );

    AMP::LinearAlgebra::Vector::shared_ptr f = meshAdapter->createVector( ThermalOutVar );


    // ThermalInpVec->zero();
    AMP::LinearAlgebra::Vector::shared_ptr ic_vector =
        AMP::LinearAlgebra::SundialsVector::view( ThermalInpVec_ic );
    bool isICVecsVarNull = ( ( ic_vector->getVariable() ).get() == 0 );
    cout << "isICVecsVarNull = " << isICVecsVarNull << endl;

    AMP::LinearAlgebra::Vector::shared_ptr ic_vector_prime =
        AMP::LinearAlgebra::SundialsVector::view( ThermalInpVec_ic_prime );


    AMP::LinearAlgebra::Vector::shared_ptr ic_vector_top =
        AMP::dynamic_pointer_cast<AMP::LinearAlgebra::Vector>( ic_vector );

    AMP::Mesh::MeshManager::Adapter::NodeIterator node = srcOpParams->d_MeshAdapter->beginNode();
    AMP::Mesh::MeshManager::Adapter::NodeIterator end_node = srcOpParams->d_MeshAdapter->endNode();

    // AMP::Mesh::DOFMap::shared_ptr dof_map = srcOpParams->d_MeshAdapter->getDOFMap(ThermalInpVar);
    AMP::Mesh::DOFMap::shared_ptr dof_map = meshAdapter->getDOFMap( ThermalInpVar );

    int counter = 0;
    for ( ; node != end_node; ++node ) {
        counter += 1;

        std::vector<unsigned int> bndGlobalIds;
        std::vector<unsigned int> d_dofIds;
        d_dofIds.resize( 0 );
        // d_dofIds[0]=1;
        dof_map->getDOFs( *node, bndGlobalIds, d_dofIds );


        double px = ( *node ).x();
        double py = ( *node ).y();
        double pz = ( *node ).z();

        double val = __INIT_FN__( px, py, pz, 0 );
        cout << "val = " << val << endl;

        cout << "counter = " << counter << "bndGlobalIds.size() = " << bndGlobalIds.size() << endl;
        for ( int i = 0; i < bndGlobalIds.size(); i++ ) {
            ThermalInpVec_ic->setValueByGlobalID( bndGlobalIds[i], val );
            // ThermalInpVec_ic_prime->setValueByGlobalID(bndGlobalIds[i], 0.0);
            // PowerInWattsVec->setValueByGlobalID(bndGlobalIds[i], val);
        } // end for i
    }     // end for node
    ThermalInpVec_ic_prime->zero();

    cout << "PowerInWattsVec->max() = " << PowerInWattsVec->max() << endl;
    cout << "PowerInWattsVec->min() = " << PowerInWattsVec->min() << endl;

    cout << "counter, the number of the nodes = " << counter << endl;

    cout << "ic_vector->max() = " << ic_vector->max() << endl;

    AMP::shared_ptr<AMP::Solver::TrilinosMLSolver> pcSolver(
        new AMP::Solver::TrilinosMLSolver( pcSolverParams ) );

    if ( pcSolver.get() == NULL ) {
        ut->failure( "Testing TrilinosMLSolver's constructor: FAIL" );
    } else {
        utvpasses( "Testing TrilinosMLSolver's constructor: PASS" );
    }


    AMP::shared_ptr<AMP::TimeIntegrator::IDATimeIntegratorParameters> time_Params(
        new AMP::TimeIntegrator::IDATimeIntegratorParameters( ida_db ) );

    if ( ( time_Params.get() ) == NULL ) {
        ut->failure( "Testing IDATimeIntegratorParameters' Constructor" );
    } else {
        ut->passes( "Testing IDATimeIntegratorParameters' Constructor" );
    }

    if ( ( ( time_Params->d_db ).get() ) != ( ida_db.get() ) ) {
        ut->failure( "Testing IDATimeIntegratorParameters::d_db" );

        ut->passes( "Testing IDATimeIntegratorParameters::d_db" );
    }

    if ( ( time_Params->d_db )->keyExists( "relative_tolerance" ) ) {
        ut->passes( "Testing IDATimeIntegratorParameters::d_db keyExists" );
    } else {
        ut->failure( "Testing IDATimeIntegratorParameters::d_db keyExists" );
    }
    if ( ( time_Params->d_db )->keyExists( "initial_time" ) ) {
        ut->passes( "Testing IDATimeIntegratorParameters::d_db keyExists" );
    } else {
        ut->failure( "Testing IDATimeIntegratorParameters::d_db keyExists" );
    }

    time_Params->d_ic_vector       = ic_vector;
    time_Params->d_ic_vector_prime = ic_vector_prime;
    time_Params->d_operator        = diffusionOperator;
    time_Params->d_pMassOperator   = massOperator;
    time_Params->d_pPreconditioner = pcSolver;
    time_Params->d_object_name     = "IDATimeIntegratorParameters";
    RightHandSideVec->zero();
    time_Params->d_pSourceTerm = RightHandSideVec;

    // time_Params->d_operator = thermOp;
    // time_Params->d_operator = diffusionOperator;

    cout << "Before IDATimeIntegrator" << endl;
    AMP::TimeIntegrator::IDATimeIntegrator *pIDAOp =
        new AMP::TimeIntegrator::IDATimeIntegrator( time_Params );

    if ( pIDAOp == NULL ) {
        ut->failure( "Testing IDATimeIntegrator's constructor" );
    } else {
        ut->passes( "Tested IDATimeIntegrator's constructor" );
    }


    AMP::shared_ptr<AMP::TimeIntegrator::IDATimeIntegrator> idaOp( pIDAOp );

    if ( idaOp == NULL ) {
        ut->failure( "Testing IDATimeIntegrator's constructor: part 2" );
    } else {
        ut->passes( "Tested IDATimeIntegrator's constructor: part 2" );
    }


    // int retval = AMP::TimeIntegrator::IDATimeIntegrator::advanceSolution(0.1, 0);
    // error: cannot call member function ‘virtual int
    // AMP::TimeIntegrator::IDATimeIntegrator::advanceSolution(double,
    // bool)’ without object

    int retval          = 0;
    double current_time = 0;
    double max          = 0;
    double abs_error    = 0.0;
    double min          = 0;
    double rel_error    = 0.0;
    double exact_sol    = 0.0;
    /*
    int j=1;
    while(idaOp->getCurrentTime() < idaOp->getFinalTime())
    {
        retval = idaOp->advanceSolution(idaOp->getCurrentDt(), 0);
        //pIDATimeIntegrator->updateSolution();
        current_time = idaOp->getCurrentTime();

        cout << j++ << "-th timestep" << endl;
        if(retval == 0) {
            ut.passes("Testing IDATimeIntegrator's advanceSolution. PASS!!");
        } else {
            ut.failure("Tested IDATimeIntegrator's advanceSolution. FAIL!!");
        }

        max = idaOp->getCurrentSolution()->max();
        min = idaOp->getCurrentSolution()->min();

        //      exact_sol = exp(-0.015 *__PI__ * __PI__ * current_time);
        //exact_sol = exp(
        //      abs_error = exact_sol-max;
        //      rel_error = abs_error/exact_sol;

        cout << "current_time = " << current_time << endl;
        cout << "max val of the current solution = " << max << endl;
        cout << "min val of the current solution = " << min << endl;
        //      cout << "exact solution = " << exact_sol << endl;
        //      cout << "absolute error = " << abs_error << endl;
        //      cout << "relative error = " << rel_error << endl;
    }
   */

    for ( int j = 0; j < 100; j++ ) {
        retval = idaOp->advanceSolution( 0.1, 0 );
        // idaOp->updateSolution();
        current_time = idaOp->getCurrentTime();

        cout << j << "th try" << endl;
        if ( retval == 0 ) {
            ut->passes( "Testing IDATimeIntegrator's advanceSolution. PASS!!" );
        } else {
            ut->failure( "Tested IDATimeIntegrator's advanceSolution. FAIL!!" );
        }

        // double max1 = idaOp->getCurrentSolution()->max();
        max = idaOp->getCurrentSolution()->max();
        min = idaOp->getCurrentSolution()->min();

        exact_sol = exp( -0.015 * __PI__ * __PI__ * current_time );
        // exact_sol = exp(
        abs_error = exact_sol - max;
        rel_error = abs_error / exact_sol;

        cout << "current_time = " << current_time << endl;
        cout << "max val of the current solution = " << max << endl;
        cout << "min val of the current solution = " << min << endl;
        cout << "exact solution = " << exact_sol << endl;
        cout << "absolute error = " << abs_error << endl;
        cout << "relative error = " << rel_error << endl;
    }


    AMP::AMPManager::shutdown();

    if ( ut->NumFailLocal() == 0 ) {
        ut->passes( "testIDATimeIntegrator successful" );
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
