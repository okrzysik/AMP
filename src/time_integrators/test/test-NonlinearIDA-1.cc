#include "materials/Material.h"
#include "operators/NeutronicsRhs.h"
#include "utils/AMPManager.h"
#include "utils/AMP_MPI.h"
#include "utils/Database.h"
#include "utils/InputDatabase.h"
#include "utils/InputManager.h"
#include "utils/PIO.h"
#include "utils/UnitTest.h"
#include "utils/Utilities.h"
#include "utils/shared_ptr.h"
#include "vectors/Variable.h"
#include <string>

#include "operators/diffusion/DiffusionLinearElement.h"
#include "operators/diffusion/DiffusionNonlinearElement.h"
#include "operators/diffusion/DiffusionTransportModel.h"
#include "utils/Writer.h"
#include "vectors/Vector.h"

#include "operators/NeutronicsRhs.h"
#include "operators/libmesh/VolumeIntegralOperator.h"

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

#include "operators/diffusion/DiffusionLinearFEOperator.h"
#include "operators/diffusion/DiffusionNonlinearFEOperator.h"
#include "operators/libmesh/MassLinearElement.h"
#include "operators/libmesh/MassLinearFEOperator.h"

#include "operators/DirichletMatrixCorrection.h"
#include "operators/DirichletVectorCorrection.h"
#include "operators/IsotropicElasticModel.h"
#include "operators/LinearBVPOperator.h"
#include "operators/MechanicsLinearElement.h"
#include "operators/MechanicsLinearFEOperator.h"

#include "solvers/trilinos/TrilinosMLSolver.h"
#include "time_integrators/LinearTimeOperator.h"
#include "time_integrators/sundials/IDATimeIntegrator.h"
#include "time_integrators/sundials/IDATimeOperator.h"

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
#define __INIT_FN__( x, y, z, t )                                                               \
    ( 400 +                                                                                     \
      exp( -0.015 * __PI__ * __PI__ * t ) * cos( 0.1 * __PI__ * x ) * cos( 0.1 * __PI__ * y ) * \
          cos( 0.05 * __PI__ * z ) )

void IDATimeIntegratorTest( AMP::UnitTest *ut )
{
    std::string input_file = "input-NonlinearIDA-1";
    std::string log_file   = "output-NonlinearIDA-1";

    AMP::AMP_MPI::initialize();
    AMP::AMPManager::startup();

    AMP::PIO::logOnlyNodeZero( log_file );

    AMP::shared_ptr<AMP::InputDatabase> input_db( new AMP::InputDatabase( "input_db" ) );
    AMP::InputManager::getManager()->parseInputFile( input_file, input_db );
    input_db->printClassData( AMP::plog );


    AMP::Mesh::MeshManagerParameters::shared_ptr mgrParams(
        new AMP::Mesh::MeshManagerParameters( input_db ) );
    AMP::Mesh::MeshManager manager( mgrParams );


    AMP_INSIST( input_db->keyExists( "DiffusionLinearElement" ),
                "Key ''DiffusionLinearElement'' is missing!" );
    AMP::shared_ptr<AMP::Database> elemOp_db = input_db->getDatabase( "DiffusionLinearElement" );
    AMP::shared_ptr<AMP::Operator::ElementOperationParameters> elemOpParams(
        new AMP::Operator::ElementOperationParameters( elemOp_db ) );
    AMP::shared_ptr<AMP::Operator::DiffusionLinearElement> thermLinElem(
        new AMP::Operator::DiffusionLinearElement( elemOpParams ) );

    AMP_INSIST( input_db->keyExists( "DiffusionNonlinearElement" ),
                "Key ''DiffusionNonlinearElement'' is missing!" );
    AMP::shared_ptr<AMP::Database> nonlinear_elemOp_db =
        input_db->getDatabase( "DiffusionNonlinearElement" );
    AMP::shared_ptr<AMP::Operator::ElementOperationParameters> nonlinear_elemOpParams(
        new AMP::Operator::ElementOperationParameters( nonlinear_elemOp_db ) );
    AMP::shared_ptr<AMP::Operator::DiffusionNonlinearElement> thermNonlinElem(
        new AMP::Operator::DiffusionNonlinearElement( nonlinear_elemOpParams ) );


    AMP_INSIST( input_db->keyExists( "DiffusionTransportModel_Nonlinear" ),
                "Key ''DiffusionTransportModel_Nonlinear'' is missing!" );
    AMP::shared_ptr<AMP::Database> thermModel_db_nonlinear =
        input_db->getDatabase( "DiffusionTransportModel_Nonlinear" );
    AMP::shared_ptr<AMP::Operator::DiffusionTransportModelParameters> thermModelParams_nonlinear(
        new AMP::Operator::DiffusionTransportModelParameters( thermModel_db_nonlinear ) );
    AMP::shared_ptr<AMP::Operator::DiffusionTransportModel> thermalDiffnModel_nonlinear(
        new AMP::Operator::DiffusionTransportModel( thermModelParams_nonlinear ) );

    AMP_INSIST( input_db->keyExists( "DiffusionTransportModel_Linear" ),
                "Key ''DiffusionTransportModel_Linear'' is missing!" );
    AMP::shared_ptr<AMP::Database> thermModel_db_linear =
        input_db->getDatabase( "DiffusionTransportModel_Linear" );
    AMP::shared_ptr<AMP::Operator::DiffusionTransportModelParameters> thermModelParams_linear(
        new AMP::Operator::DiffusionTransportModelParameters( thermModel_db_linear ) );
    AMP::shared_ptr<AMP::Operator::DiffusionTransportModel> thermalDiffnModel_linear(
        new AMP::Operator::DiffusionTransportModel( thermModelParams_linear ) );

    AMP_INSIST( input_db->keyExists( "DiffusionNonlinearFEOperator" ),
                "Key ''DiffusionNonlinearFEOperator'' is missing!" );
    AMP::shared_ptr<AMP::Database> thermNonlinAssembly_db =
        input_db->getDatabase( "DiffusionNonlinearFEOperator" );
    AMP::shared_ptr<AMP::Operator::DiffusionNonlinearFEOperatorParameters> thermNonlinOpParams(
        new AMP::Operator::DiffusionNonlinearFEOperatorParameters( thermNonlinAssembly_db ) );

    AMP_INSIST( input_db->keyExists( "DiffusionLinearFEOperator" ),
                "Key ''DiffusionLinearFEOperator'' is missing!" );
    AMP::shared_ptr<AMP::Database> thermLinAssembly_db =
        input_db->getDatabase( "DiffusionLinearFEOperator" );
    AMP::shared_ptr<AMP::Operator::DiffusionLinearFEOperatorParameters> thermLinOpParams(
        new AMP::Operator::DiffusionLinearFEOperatorParameters( thermLinAssembly_db ) );

    thermNonlinOpParams->d_elemOp = thermNonlinElem;
    // thermNonlinOpParams->d_elemOp_linear = thermLinElem;
    thermNonlinOpParams->d_transportModel = thermalDiffnModel_nonlinear;
    thermNonlinOpParams->d_MeshAdapter    = manager.getMesh( "ida" );

    thermLinOpParams->d_elemOp         = thermLinElem;
    thermLinOpParams->d_transportModel = thermalDiffnModel_linear;
    thermLinOpParams->d_MeshAdapter    = manager.getMesh( "ida" );

    AMP::shared_ptr<AMP::Operator::DiffusionNonlinearFEOperator> thermNonlinOp(
        new AMP::Operator::DiffusionNonlinearFEOperator( thermNonlinOpParams ) );
    AMP::shared_ptr<AMP::Operator::DiffusionLinearFEOperator> thermLinOp(
        new AMP::Operator::DiffusionLinearFEOperator( thermLinOpParams ) );


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
    AMP::shared_ptr<AMP::Operator::MassDensityModel> massDensityModel(
        new AMP::Operator::MassDensityModel( massModelParams ) );

    AMP_INSIST( input_db->keyExists( "MassLinearFEOperator" ),
                "Key ''MassLinearFEOperator'' is missing!" );
    AMP::shared_ptr<AMP::Database> srcAssembly_db = input_db->getDatabase( "MassLinearFEOperator" );
    AMP::shared_ptr<AMP::Operator::MassLinearFEOperatorParameters> srcOpParams(
        new AMP::Operator::MassLinearFEOperatorParameters( srcAssembly_db ) );

    // ida
    AMP_INSIST( input_db->keyExists( "IDATimeIntegrator" ),
                "Key ''IDATimeIntegrator'' is missing!" );
    AMP::shared_ptr<AMP::Database> ida_db = input_db->getDatabase( "IDATimeIntegrator" );

    // idaTimeOperator
    AMP::shared_ptr<AMP::Database> idaTimeOp_db    = input_db->getDatabase( "IDATimeOperator" );
    AMP::shared_ptr<AMP::Database> linearTimeOp_db = input_db->getDatabase( "LinearTimeOperator" );


    // Copied and pasted from testPetscSNESSolver-NonlinearThermal-1.cc
    ////////////////////////////////////
    //  CREATE THE NEUTRONICS SOURCE  //
    ////////////////////////////////////
    AMP::LinearAlgebra::Vector::shared_ptr nullVec;

    AMP_INSIST( input_db->keyExists( "NeutronicsOperator" ),
                "Key ''NeutronicsOperator'' is missing!" );
    AMP::shared_ptr<AMP::Database> neutronicsOp_db = input_db->getDatabase( "NeutronicsOperator" );
    AMP::shared_ptr<AMP::Operator::NeutronicsRhsParameters> neutronicsParams(
        new AMP::Operator::NeutronicsRhsParameters( neutronicsOp_db ) );
    // neutronicsParams->d_MeshAdapter = meshAdapter;
    AMP::Mesh::MeshManager::Adapter::shared_ptr meshAdapter = manager.getMesh( "ida" );
    neutronicsParams->d_MeshAdapter                         = manager.getMesh( "ida" );

    AMP::shared_ptr<AMP::Operator::NeutronicsRhs> neutronicsOperator(
        new AMP::Operator::NeutronicsRhs( neutronicsParams ) );

    AMP::LinearAlgebra::Variable::shared_ptr SpecificPowerVar =
        neutronicsOperator->getOutputVariable();
    AMP::LinearAlgebra::Vector::shared_ptr SpecificPowerVec =
        meshAdapter->createVector( SpecificPowerVar );

    neutronicsOperator->apply( nullVec, SpecificPowerVec );

    /////////////////////////////////////////////////////
    //  Integrate Nuclear Rhs over Desnity * Volume //
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

    cout << "PowerInWattsVec->max() = " << PowerInWattsVec->max() << endl;
    cout << "PowerInWattsVec->min() = " << PowerInWattsVec->min() << endl;


    // PowerInWattsVec->setToScalar(100.0);

    //----------------------------------------------------------------------------------------------------------------------------------------------//

    // do we really need dynamic cast here?
    srcOpParams->d_elemOp =
        AMP::dynamic_pointer_cast<AMP::Operator::ElementOperation>( massLinElem );
    // srcOpParams->d_elemOp = massLinElem;
    srcOpParams->d_densityModel = massDensityModel;
    srcOpParams->d_MeshAdapter  = manager.getMesh( "ida" );
    AMP::shared_ptr<AMP::Operator::MassLinearFEOperator> srcOp(
        new AMP::Operator::MassLinearFEOperator( srcOpParams ) );


    AMP::shared_ptr<AMP::Database> pcSolver_db = input_db->getDatabase( "Preconditioner" );
    AMP::shared_ptr<AMP::Solver::SolverStrategyParameters> pcSolverParams(
        new AMP::Solver::SolverStrategyParameters( pcSolver_db ) );

    if ( pcSolverParams.get() == NULL ) {
        ut->failure( "Testing SolverStrategyParameters's constructor: FAIL" );
    } else {
        ut->passes( "Testing SolverStrategyParameters's constructor: PASS" );
    }


    // AMP::Mesh::MeshManager::Adapter::shared_ptr  meshAdapter = manager.getMesh("ida");
    // AMP::LinearAlgebra::Variable::shared_ptr ThermalInpVar =
    // thermNonlinOp->createInputVariable("ThermalInputVariable",thermNonlinOpParams->d_PrincipalVariable);
    AMP::LinearAlgebra::Variable::shared_ptr ThermalInpVar =
        thermNonlinOp->getInputVariable( thermNonlinOp->getPrincipalVariableId() );
    // thermNonlinOp->setInputVariableName("ThermalInputVariable",thermNonlinOpParams->d_PrincipalVariable);
    AMP::LinearAlgebra::Variable::shared_ptr ThermalOutVar =
        thermNonlinOp->createOutputVariable( "ThermalOutputVariable", 0 );
    thermNonlinOp->setOutputVariableName( "ThermalOutputVariable" );

    AMP::LinearAlgebra::Variable::shared_ptr MassInpVar = srcOp->getInputVariable();
    AMP::LinearAlgebra::Variable::shared_ptr MassOutVar = srcOp->getOutputVariable();


    AMP::LinearAlgebra::Vector::shared_ptr ThermalInpVec_ic =
        meshAdapter->createVector( ThermalInpVar );
    AMP::LinearAlgebra::Vector::shared_ptr ThermalInpVec_ic_prime =
        meshAdapter->createVector( MassInpVar );
    AMP::LinearAlgebra::Vector::shared_ptr ThermalInpVec_ic_temp =
        meshAdapter->createVector( ThermalOutVar );
    AMP::LinearAlgebra::Vector::shared_ptr f = meshAdapter->createVector( ThermalInpVar );

    AMP::LinearAlgebra::Vector::shared_ptr ic_vector =
        AMP::LinearAlgebra::SundialsVector::view( ThermalInpVec_ic );

    AMP::LinearAlgebra::Vector::shared_ptr ic_vector_prime =
        AMP::LinearAlgebra::SundialsVector::view( ThermalInpVec_ic_prime );
    AMP::LinearAlgebra::Vector::shared_ptr temp_vector =
        AMP::LinearAlgebra::SundialsVector::view( ThermalInpVec_ic_temp );


    AMP::Mesh::MeshManager::Adapter::NodeIterator node = srcOpParams->d_MeshAdapter->beginNode();
    AMP::Mesh::MeshManager::Adapter::NodeIterator end_node = srcOpParams->d_MeshAdapter->endNode();

    AMP::Mesh::DOFMap::shared_ptr dof_map = meshAdapter->getDOFMap( ThermalInpVar );
    // int counter=0;
    for ( ; node != end_node; ++node ) {
        // counter+=1;

        std::vector<unsigned int> bndGlobalIds;
        std::vector<unsigned int> d_dofIds;
        d_dofIds.resize( 0 );
        // d_dofIds[0]=1;
        dof_map->getDOFs( *node, bndGlobalIds, d_dofIds );


        double px = ( *node ).x();
        double py = ( *node ).y();
        double pz = ( *node ).z();

        double val = __INIT_FN__( px, py, pz, 0 );
        // cout << "val = " << val << endl;

        // cout << "counter = " << counter << "bndGlobalIds.size() = " << bndGlobalIds.size() <<
        // endl;
        for ( int i = 0; i < bndGlobalIds.size(); i++ ) {
            ThermalInpVec_ic->setValueByGlobalID( bndGlobalIds[i], val );
            ThermalInpVec_ic_prime->setValueByGlobalID( bndGlobalIds[i], 0.0 );

            PowerInWattsVec->setValueByGlobalID( bndGlobalIds[i], val );
        } // end for i
    }     // end for node

    /*
    AMP::Mesh::MeshManager::Adapter::OwnedBoundaryNodeIterator curNode =
    meshAdapter->beginOwnedBoundary (8 );
        while ( curNode != meshAdapter->endOwnedBoundary (8 ) )
            {
                    PowerInWattsVec->setValueByGlobalID ( dof_map->getGlobalID ( curNode->globalID()
    , 0 ),0 );
                    curNode++;
            }
    */
    // pcSolverParams->d_pOperator = thermLinOp;

    AMP::shared_ptr<AMP::Solver::TrilinosMLSolver> pcSolver(
        new AMP::Solver::TrilinosMLSolver( pcSolverParams ) );

    if ( pcSolver.get() == NULL ) {
        ut->failure( "Testing TrilinosMLSolver's constructor: FAIL" );
    } else {
        ut->passes( "Testing TrilinosMLSolver's constructor: PASS" );
    }

    AMP::shared_ptr<AMP::Operator::OperatorParameters> TimeOp_Params(
        new AMP::TimeIntegrator::IDATimeOperatorParameters( idaTimeOp_db ) );
    AMP::shared_ptr<AMP::TimeIntegrator::IDATimeOperatorParameters> idaTimeOp_Params =
        AMP::dynamic_pointer_cast<AMP::TimeIntegrator::IDATimeOperatorParameters>( TimeOp_Params );
    idaTimeOp_Params->d_pRhsOperator  = thermNonlinOp;
    idaTimeOp_Params->d_pMassOperator = srcOp;

    AMP::shared_ptr<AMP::TimeIntegrator::IDATimeOperator> idaTimeOp(
        new AMP::TimeIntegrator::IDATimeOperator( TimeOp_Params ) );

    AMP::shared_ptr<AMP::TimeIntegrator::TimeOperatorParameters> linearTimeOp_Params(
        new AMP::TimeIntegrator::TimeOperatorParameters( linearTimeOp_db ) );

    AMP::shared_ptr<AMP::TimeIntegrator::TimeIntegratorParameters> time_Params(
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

    time_Params->d_ic_vector = ic_vector;

    time_Params->d_pMassOperator = srcOp;
    AMP::shared_ptr<AMP::TimeIntegrator::IDATimeIntegratorParameters> ida_Params =
        AMP::dynamic_pointer_cast<AMP::TimeIntegrator::IDATimeIntegratorParameters>( time_Params );
    ida_Params->d_ic_vector_prime               = ic_vector_prime;
    ida_Params->d_pPreconditioner               = pcSolver;
    ida_Params->d_pLinearOperator               = thermLinOp;
    ida_Params->d_pLinearTimeOperatorParameters = linearTimeOp_Params;
    time_Params->d_pSourceTerm                  = PowerInWattsVec;

    time_Params->d_object_name = "whatever";
    time_Params->d_operator    = thermNonlinOp;


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
        utvpasses( "Tested IDATimeIntegrator's constructor: part 2" );
    }


    int retval          = 0;
    double current_time = 0;
    double max          = 0;
    double min          = 0;

    for ( int j = 0; j < 100; j++ ) {
        retval       = idaOp->advanceSolution( 0.1, 0 );
        current_time = idaOp->getCurrentTime();

        cout << j << "th try" << endl;
        if ( retval == 0 ) {
            ut->passes( "Testing IDATimeIntegrator's advanceSolution. PASS!!" );
        } else {
            ut->failure( "Tested IDATimeIntegrator's advanceSolution. FAIL!!" );
        }

        max = idaOp->getCurrentSolution()->max();
        min = idaOp->getCurrentSolution()->min();

        cout << "current_time = " << current_time << endl;
        cout << "max val of the current solution = " << max << endl;
        cout << "min val of the current solution = " << min << endl;
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

    try {
        IDATimeIntegratorTest( &ut );
    } catch ( std::exception &err ) {
        std::cout << "ERROR: While testing " << argv[0] << err.what() << std::endl;
        ut.failure( "ERROR: While testing" );
    } catch ( ... ) {
        std::cout << "ERROR: While testing " << argv[0] << "An unknown exception was thrown."
                  << std::endl;
        ut.failure( "ERROR: While testing" );
    }

    ut.report();

    int num_failed = ut.NumFailGlobal();
    AMP::AMPManager::shutdown();
    return num_failed;
}

//---------------------------------------------------------------------------//
//                        end of SundialsVectorTest.cc
//---------------------------------------------------------------------------//
