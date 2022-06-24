#include "AMP/IO/PIO.h"
#include "AMP/discretization/DOF_Manager.h"
#include "AMP/discretization/simpleDOF_Manager.h"
#include "AMP/mesh/Mesh.h"
#include "AMP/mesh/MeshFactory.h"
#include "AMP/mesh/MeshParameters.h"
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
#include "AMP/utils/UnitTest.h"
#include "AMP/utils/Utilities.h"
#include "AMP/vectors/VectorBuilder.h"

#include "applyTests.h"

#include <iostream>
#include <memory>
#include <string>


static void nonlinearTest( AMP::UnitTest *ut, const std::string &exeName )
{
    // Initialization
    std::string input_file = "input_" + exeName;
    std::string log_file   = "output_" + exeName;

    AMP::logOnlyNodeZero( log_file );

    std::cout << "testing with input file " << input_file << std::endl;
    std::cout.flush();

    // Test create
    auto input_db = AMP::Database::parseInputFile( input_file );
    input_db->print( AMP::plog );

    // Create the Mesh.
    AMP_INSIST( input_db->keyExists( "Mesh" ), "Key ''Mesh'' is missing!" );
    auto mesh_db   = input_db->getDatabase( "Mesh" );
    auto mgrParams = std::make_shared<AMP::Mesh::MeshParameters>( mesh_db );
    mgrParams->setComm( AMP::AMP_MPI( AMP_COMM_WORLD ) );
    auto meshAdapter = AMP::Mesh::MeshFactory::create( mgrParams );

    std::shared_ptr<AMP::Operator::ElementPhysicsModel> elementModel;
    auto fsOp_db           = input_db->getDatabase( "NonlinearFickSoretOp" );
    auto nonlinearOperator = AMP::Operator::OperatorBuilder::createOperator(
        meshAdapter, "NonlinearFickSoretOp", input_db, elementModel );
    auto fsOp =
        std::dynamic_pointer_cast<AMP::Operator::FickSoretNonlinearFEOperator>( nonlinearOperator );
    auto fickOp  = fsOp->getFickOperator();
    auto soretOp = fsOp->getSoretOperator();

    ut->passes( exeName + ": create" );
    std::cout.flush();

    // set up defaults for materials arguments and create transport model
    auto fickModel  = fickOp->getTransportModel();
    auto soretModel = soretOp->getTransportModel();

    // create parameters for reset test
    auto fickOpParams =
        std::make_shared<AMP::Operator::DiffusionNonlinearFEOperatorParameters>( fsOp_db );
    auto soretOpParams =
        std::make_shared<AMP::Operator::DiffusionNonlinearFEOperatorParameters>( fsOp_db );
    fickOpParams->d_transportModel  = fickModel;
    soretOpParams->d_transportModel = soretModel;
    auto fsOpBase_db                = fsOp_db;
    auto fsOpParams =
        std::make_shared<AMP::Operator::FickSoretNonlinearFEOperatorParameters>( fsOpBase_db );
    fsOpParams->d_FickParameters  = fickOpParams;
    fsOpParams->d_SoretParameters = soretOpParams;

    // create vectors for parameters
    auto tVar = std::make_shared<AMP::LinearAlgebra::Variable>( "temp" );
    auto cVar = std::make_shared<AMP::LinearAlgebra::Variable>( "conc" );
    auto bVar = std::make_shared<AMP::LinearAlgebra::Variable>( "burnup" );

    // Create a DOF manager for a nodal vector
    int DOFsPerNode     = 1;
    int nodalGhostWidth = 1;
    bool split          = true;
    auto nodalDofMap    = AMP::Discretization::simpleDOFManager::create(
        meshAdapter, AMP::Mesh::GeomType::Vertex, nodalGhostWidth, DOFsPerNode, split );

    // create solution, rhs, and residual vectors
    auto tVec = AMP::LinearAlgebra::createVector( nodalDofMap, tVar );
    auto cVec = AMP::LinearAlgebra::createVector( nodalDofMap, cVar );
    auto bVec = AMP::LinearAlgebra::createVector( nodalDofMap, bVar );

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
    std::vector<double> defaults;
    auto matFick  = fickModel->getMaterial();  // compile error
    auto matSoret = soretModel->getMaterial(); // compile error
    // the Soret has a principal variable of temperature
    if ( soretOp->getPrincipalVariableId() == AMP::Operator::Diffusion::TEMPERATURE ) {
        std::string property = "ThermalDiffusionCoefficient";
        if ( ( matSoret->property( property ) )->is_argument( "temperature" ) ) {
            auto range =
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
            auto range = ( matFick->property( property ) )
                             ->get_arg_range( "concentration" ); // Compile error
            scale[1] = range[1] - range[0];
            shift[1] = range[0] + 0.001 * scale[1];
            scale[1] *= 0.999;
            defaults = ( matFick->property( property ) )->get_defaults(); // compile error
        }
    }
    if ( defaults.size() > 0 )
        tVec->setToScalar( defaults[0] ); // compile error
    if ( defaults.size() > 1 )
        cVec->setToScalar( defaults[1] ); // compile error
    if ( defaults.size() > 2 )
        bVec->setToScalar( defaults[2] ); // compile error
    // set up input multivariable and output variable
    auto fsInpVar = std::make_shared<AMP::LinearAlgebra::MultiVariable>( "fsInput" );
    fsInpVar->add( tVar );
    fsInpVar->add( cVar );
    fsInpVar->add( bVar );
    auto fsOutVar = fickOp->getOutputVariable();

    std::string msgPrefix = exeName + ": apply ";
    auto solVec           = AMP::LinearAlgebra::createVector( nodalDofMap, fsInpVar );
    auto rhsVec           = AMP::LinearAlgebra::createVector( nodalDofMap, fsOutVar );
    auto resVec           = AMP::LinearAlgebra::createVector( nodalDofMap, fsOutVar );

    // set default values of input variables
    auto inTempVec = solVec->subsetVectorForVariable( tVar );
    auto inConcVec = solVec->subsetVectorForVariable( cVar );
    auto inBurnVec = solVec->subsetVectorForVariable( bVar );
    if ( defaults.size() > 0 )
        inTempVec->setToScalar( defaults[0] ); // compile error
    if ( defaults.size() > 1 )
        inConcVec->setToScalar( defaults[1] ); // compile error
    if ( defaults.size() > 2 )
        inBurnVec->setToScalar( defaults[2] ); // compile error

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
