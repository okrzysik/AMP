#include "AMP/ampmesh/Mesh.h"
#include "AMP/discretization/DOF_Manager.h"
#include "AMP/discretization/simpleDOF_Manager.h"
#include "AMP/materials/Material.h"
#include "AMP/operators/ElementPhysicsModelFactory.h"
#include "AMP/operators/ElementPhysicsModelParameters.h"
#include "AMP/operators/OperatorBuilder.h"
#include "AMP/operators/diffusion/DiffusionConstants.h"
#include "AMP/operators/diffusion/DiffusionLinearElement.h"
#include "AMP/operators/diffusion/DiffusionLinearFEOperator.h"
#include "AMP/operators/diffusion/DiffusionLinearFEOperatorParameters.h"
#include "AMP/operators/diffusion/DiffusionNonlinearElement.h"
#include "AMP/operators/diffusion/DiffusionNonlinearFEOperator.h"
#include "AMP/operators/diffusion/DiffusionNonlinearFEOperatorParameters.h"
#include "AMP/operators/diffusion/DiffusionTransportModel.h"
#include "AMP/operators/diffusion/FickSoretNonlinearFEOperator.h"
#include "AMP/operators/diffusion/FickSoretNonlinearFEOperatorParameters.h"
#include "AMP/utils/AMPManager.h"
#include "AMP/utils/AMP_MPI.h"
#include "AMP/utils/Database.h"
#include "AMP/utils/PIO.h"
#include "AMP/utils/UnitTest.h"
#include "AMP/utils/Utilities.h"
#include "AMP/vectors/VectorBuilder.h"
#include <memory>

#include "applyTests.h"

#include <iostream>
#include <string>


static void nonlinearTest( AMP::UnitTest *ut, const std::string &exeName )
{
    // Initialization
    std::string input_file = "input_" + exeName;
    std::string log_file   = "output_" + exeName;

    AMP::PIO::logOnlyNodeZero( log_file );

    std::cout << "testing with input file " << input_file << std::endl;
    std::cout.flush();

    // Test create

    auto input_db = AMP::Database::parseInputFile( input_file );
    input_db->print( AMP::plog );

    // Create the Mesh.
    AMP_INSIST( input_db->keyExists( "Mesh" ), "Key ''Mesh'' is missing!" );
    std::shared_ptr<AMP::Database> mesh_db = input_db->getDatabase( "Mesh" );
    std::shared_ptr<AMP::Mesh::MeshParameters> mgrParams(
        new AMP::Mesh::MeshParameters( mesh_db ) );
    mgrParams->setComm( AMP::AMP_MPI( AMP_COMM_WORLD ) );
    std::shared_ptr<AMP::Mesh::Mesh> meshAdapter = AMP::Mesh::Mesh::buildMesh( mgrParams );

    std::shared_ptr<AMP::Operator::FickSoretNonlinearFEOperator> fsOp;
    std::shared_ptr<AMP::Operator::ElementPhysicsModel> elementModel;
    std::shared_ptr<AMP::Database> fsOp_db =
        std::dynamic_pointer_cast<AMP::Database>( input_db->getDatabase( "NonlinearFickSoretOp" ) );
    std::shared_ptr<AMP::Operator::Operator> nonlinearOperator =
        AMP::Operator::OperatorBuilder::createOperator(
            meshAdapter, "NonlinearFickSoretOp", input_db, elementModel );
    fsOp =
        std::dynamic_pointer_cast<AMP::Operator::FickSoretNonlinearFEOperator>( nonlinearOperator );
    std::shared_ptr<AMP::Operator::DiffusionNonlinearFEOperator> fickOp  = fsOp->getFickOperator();
    std::shared_ptr<AMP::Operator::DiffusionNonlinearFEOperator> soretOp = fsOp->getSoretOperator();

    ut->passes( exeName + ": create" );
    std::cout.flush();

    // set up defaults for materials arguments and create transport model
    std::shared_ptr<AMP::Operator::DiffusionTransportModel> fickModel = fickOp->getTransportModel();
    std::shared_ptr<AMP::Operator::DiffusionTransportModel> soretModel =
        soretOp->getTransportModel();

    // create parameters for reset test
    std::shared_ptr<AMP::Operator::DiffusionNonlinearFEOperatorParameters> fickOpParams(
        new AMP::Operator::DiffusionNonlinearFEOperatorParameters( fsOp_db ) );
    std::shared_ptr<AMP::Operator::DiffusionNonlinearFEOperatorParameters> soretOpParams(
        new AMP::Operator::DiffusionNonlinearFEOperatorParameters( fsOp_db ) );
    fickOpParams->d_transportModel  = fickModel;
    soretOpParams->d_transportModel = soretModel;
    std::shared_ptr<AMP::Database> fsOpBase_db(
        std::dynamic_pointer_cast<AMP::Database>( fsOp_db ) );
    std::shared_ptr<AMP::Operator::FickSoretNonlinearFEOperatorParameters> fsOpParams(
        new AMP::Operator::FickSoretNonlinearFEOperatorParameters( fsOpBase_db ) );
    fsOpParams->d_FickParameters  = fickOpParams;
    fsOpParams->d_SoretParameters = soretOpParams;

    // create vectors for parameters
    AMP::LinearAlgebra::Variable::shared_ptr tVar( new AMP::LinearAlgebra::Variable( "temp" ) );
    AMP::LinearAlgebra::Variable::shared_ptr cVar( new AMP::LinearAlgebra::Variable( "conc" ) );
    AMP::LinearAlgebra::Variable::shared_ptr bVar( new AMP::LinearAlgebra::Variable( "burnup" ) );

    // Create a DOF manager for a nodal vector
    int DOFsPerNode     = 1;
    int nodalGhostWidth = 1;
    bool split          = true;
    AMP::Discretization::DOFManager::shared_ptr nodalDofMap =
        AMP::Discretization::simpleDOFManager::create(
            meshAdapter, AMP::Mesh::GeomType::Vertex, nodalGhostWidth, DOFsPerNode, split );

    // create solution, rhs, and residual vectors
    AMP::LinearAlgebra::Vector::shared_ptr tVec =
        AMP::LinearAlgebra::createVector( nodalDofMap, tVar );
    AMP::LinearAlgebra::Vector::shared_ptr cVec =
        AMP::LinearAlgebra::createVector( nodalDofMap, cVar );
    AMP::LinearAlgebra::Vector::shared_ptr bVec =
        AMP::LinearAlgebra::createVector( nodalDofMap, bVar );

    // set vectors in parameters
    fickOpParams->d_FrozenTemperature    = tVec;
    fickOpParams->d_FrozenConcentration  = cVec;
    fickOpParams->d_FrozenBurnup         = bVec;
    soretOpParams->d_FrozenTemperature   = tVec;
    soretOpParams->d_FrozenConcentration = cVec;
    soretOpParams->d_FrozenBurnup        = bVec;

    // Test reset
    {
        fsOp->reset( fsOpParams );
        ut->passes( exeName + ": reset" );
        std::cout.flush();
    }

    // set shift and scale for applyTests
    double shift[2];
    double scale[2];
    shift[0] = 0.;
    shift[1] = 0.;
    scale[0] = 1.;
    scale[1] = 1.;
    std::vector<double> range( 2 );
    std::vector<double> defaults;
    AMP::Materials::Material::shared_ptr matFick  = fickModel->getMaterial();  // compile error
    AMP::Materials::Material::shared_ptr matSoret = soretModel->getMaterial(); // compile error
    // the Soret has a principal variable of temperature
    if ( soretOp->getPrincipalVariableId() == AMP::Operator::Diffusion::TEMPERATURE ) {
        std::string property = "ThermalDiffusionCoefficient";
        if ( ( matSoret->property( property ) )->is_argument( "temperature" ) ) {
            range =
                ( matSoret->property( property ) )->get_arg_range( "temperature" ); // Compile error
            scale[0] = range[1] - range[0];
            shift[0] = range[0] + 0.001 * scale[0];
            scale[0] *= 0.999;
            defaults = ( matSoret->property( property ) )->get_defaults(); // compile error
        }
    }
    // the fick has a principal variable of oxygen
    if ( fickOp->getPrincipalVariableId() == AMP::Operator::Diffusion::CONCENTRATION ) {
        std::string property = "FickCoefficient";
        if ( ( matFick->property( property ) )->is_argument( "concentration" ) ) {
            range = ( matFick->property( property ) )
                        ->get_arg_range( "concentration" ); // Compile error
            scale[1] = range[1] - range[0];
            shift[1] = range[0] + 0.001 * scale[1];
            scale[1] *= 0.999;
            defaults = ( matFick->property( property ) )->get_defaults(); // compile error
        }
    }
    if ( defaults.size() > 0 )
        tVec->setToScalar( defaults[0], tVec ); // compile error
    if ( defaults.size() > 1 )
        cVec->setToScalar( defaults[1], cVec ); // compile error
    if ( defaults.size() > 2 )
        bVec->setToScalar( defaults[2], bVec ); // compile error
    // set up input multivariable and output variable
    std::shared_ptr<AMP::LinearAlgebra::MultiVariable> fsInpVar(
        new AMP::LinearAlgebra::MultiVariable( "fsInput" ) );
    fsInpVar->add( tVar );
    fsInpVar->add( cVar );
    fsInpVar->add( bVar );
    std::shared_ptr<AMP::LinearAlgebra::Variable> fsOutVar( fickOp->getOutputVariable() );

    std::string msgPrefix = exeName + ": apply ";
    AMP::LinearAlgebra::Vector::shared_ptr solVec =
        AMP::LinearAlgebra::createVector( nodalDofMap, fsInpVar );
    AMP::LinearAlgebra::Vector::shared_ptr rhsVec =
        AMP::LinearAlgebra::createVector( nodalDofMap, fsOutVar );
    AMP::LinearAlgebra::Vector::shared_ptr resVec =
        AMP::LinearAlgebra::createVector( nodalDofMap, fsOutVar );

    // set default values of input variables
    auto inTempVec = solVec->subsetVectorForVariable( tVar );
    auto inConcVec = solVec->subsetVectorForVariable( cVar );
    auto inBurnVec = solVec->subsetVectorForVariable( bVar );
    if ( defaults.size() > 0 )
        inTempVec->setToScalar( defaults[0], inTempVec ); // compile error
    if ( defaults.size() > 1 )
        inConcVec->setToScalar( defaults[1], inConcVec ); // compile error
    if ( defaults.size() > 2 )
        inBurnVec->setToScalar( defaults[2], inBurnVec ); // compile error

    AMP_INSIST( fsOp->isValidInput( solVec ), "input variable not set up correctly" );

    // test apply
    applyTests( ut, msgPrefix, nonlinearOperator, rhsVec, solVec, resVec, shift, scale, 2 );
    std::cout.flush();
}

int testNonlinearFickSoret( int argc, char *argv[] )
{
    AMP::AMPManagerProperties startup_properties;
    startup_properties.use_MPI_Abort = false;
    AMP::AMPManager::startup( argc, argv, startup_properties );

    AMP::UnitTest ut;
    const int NUMFILES          = 2;
    std::string files[NUMFILES] = { "FickSoret-TUI-1", "FickSoret-UO2MSRZC09-1" };

    for ( auto &file : files )
        nonlinearTest( &ut, file );

    ut.report();

    int num_failed = ut.NumFailGlobal();
    AMP::AMPManager::shutdown();
    return num_failed;
}
