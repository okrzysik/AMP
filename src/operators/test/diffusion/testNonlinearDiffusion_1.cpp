#include "AMP/ampmesh/Mesh.h"
#include "AMP/ampmesh/MeshParameters.h"
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
#include "AMP/utils/AMPManager.h"
#include "AMP/utils/AMP_MPI.h"
#include "AMP/utils/Database.h"
#include "AMP/utils/PIO.h"
#include "AMP/utils/UnitTest.h"
#include "AMP/utils/Utilities.h"
#include "AMP/vectors/Variable.h"
#include "AMP/vectors/Vector.h"
#include "AMP/vectors/VectorBuilder.h"

#include "../applyTests.h"

#include <iostream>
#include <memory>
#include <string>


static void nonlinearTest( AMP::UnitTest *ut, const std::string &exeName )
{
    // Initialization
    std::string input_file = "input_" + exeName;
    std::string log_file   = "output_" + exeName;

    AMP::PIO::logOnlyNodeZero( log_file );
    AMP::AMP_MPI globalComm( AMP_COMM_WORLD );

    std::cout << "testing with input file " << input_file << std::endl;
    std::cout.flush();

    // Test create
    auto input_db = AMP::Database::parseInputFile( input_file );
    input_db->print( AMP::plog );

    // Get the Mesh database and create the mesh parameters
    auto database = input_db->getDatabase( "Mesh" );
    auto params   = std::make_shared<AMP::Mesh::MeshParameters>( database );
    params->setComm( globalComm );

    // Create the meshes from the input database
    auto meshAdapter = AMP::Mesh::Mesh::buildMesh( params );

    // nonlinear operator
    std::shared_ptr<AMP::Operator::ElementPhysicsModel> elementModel;
    auto diffFEOp_db =
        std::dynamic_pointer_cast<AMP::Database>( input_db->getDatabase( "NonlinearDiffusionOp" ) );
    auto nonlinearOperator = AMP::Operator::OperatorBuilder::createOperator(
        meshAdapter, "NonlinearDiffusionOp", input_db, elementModel );
    auto diffOp =
        std::dynamic_pointer_cast<AMP::Operator::DiffusionNonlinearFEOperator>( nonlinearOperator );

    // linear operator
    std::shared_ptr<AMP::Operator::ElementPhysicsModel> linElementModel;
    auto linearOperator = AMP::Operator::OperatorBuilder::createOperator(
        meshAdapter, "LinearDiffusionOp", input_db, linElementModel );
    auto linOp =
        std::dynamic_pointer_cast<AMP::Operator::DiffusionLinearFEOperator>( linearOperator );

    ut->passes( exeName + ": create" );
    std::cout.flush();

    // set up defaults for materials arguments and create transport model
    auto transportModel_db = input_db->getDatabase( "DiffusionTransportModel" );
    auto elementPhysicsModel =
        AMP::Operator::ElementPhysicsModelFactory::createElementPhysicsModel( transportModel_db );
    auto transportModel =
        std::dynamic_pointer_cast<AMP::Operator::DiffusionTransportModel>( elementPhysicsModel );

    double defTemp = transportModel_db->getWithDefault<double>( "Default_Temperature", 400.0 );
    double defConc = transportModel_db->getWithDefault<double>( "Default_Concentration", .33 );
    double defBurn = transportModel_db->getWithDefault<double>( "Default_Burnup", .5 );

    std::string property = transportModel_db->getString( "Property" );

    // create parameters for reset test
    auto diffOpParams =
        std::make_shared<AMP::Operator::DiffusionNonlinearFEOperatorParameters>( diffFEOp_db );

    // nullify vectors in parameters
    diffOpParams->d_FrozenTemperature.reset();
    diffOpParams->d_FrozenConcentration.reset();
    diffOpParams->d_FrozenBurnup.reset();

    // create vectors for parameters
    auto active_db = diffFEOp_db->getDatabase( "ActiveInputVariables" );
    auto tVar      = std::make_shared<AMP::LinearAlgebra::Variable>(
        active_db->getWithDefault<std::string>( "temperature", "not_specified" ) );
    auto cVar = std::make_shared<AMP::LinearAlgebra::Variable>(
        active_db->getWithDefault<std::string>( "concentration", "not_specified" ) );
    auto bVar = std::make_shared<AMP::LinearAlgebra::Variable>(
        active_db->getWithDefault<std::string>( "burnup", "not_specified" ) );

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
    tVec->setToScalar( defTemp );
    cVec->setToScalar( defConc );
    bVec->setToScalar( defBurn );

    // set principal variable vector and shift for applyTests
    double shift = 0., scale = 1.;
    auto mat = transportModel->getMaterial();
    if ( diffOp->getPrincipalVariableId() == AMP::Operator::Diffusion::TEMPERATURE ) {
        diffOpParams->d_FrozenTemperature = tVec;
        if ( ( mat->property( property ) )->is_argument( "temperature" ) ) {
            auto range =
                ( mat->property( property ) )->get_arg_range( "temperature" ); // Compile error
            scale = range[1] - range[0];
            shift = range[0] + 0.001 * scale;
            scale *= 0.999;
        }
    }
    if ( diffOp->getPrincipalVariableId() == AMP::Operator::Diffusion::CONCENTRATION ) {
        diffOpParams->d_FrozenConcentration = cVec;
        if ( ( mat->property( property ) )->is_argument( "concentration" ) ) {
            auto range =
                ( mat->property( property ) )->get_arg_range( "concentration" ); // Compile error
            scale = range[1] - range[0];
            shift = range[0] + 0.001 * scale;
            scale *= 0.999;
        }
    }
    if ( diffOp->getPrincipalVariableId() == AMP::Operator::Diffusion::BURNUP ) {
        AMP_INSIST( false, "do not know what to do" );
    }

    // set frozen vectors in parameters
    if ( diffFEOp_db->getWithDefault( "Freezetemperature", false ) )
        diffOpParams->d_FrozenTemperature = tVec;
    if ( diffFEOp_db->getWithDefault( "Freezeconcentration", false ) )
        diffOpParams->d_FrozenConcentration = cVec;
    if ( diffFEOp_db->getWithDefault( "Freezeburnup", false ) )
        diffOpParams->d_FrozenBurnup = bVec;

    // set transport model
    diffOpParams->d_transportModel = transportModel;

    // Test reset
    {
        diffOp->reset( diffOpParams );
        ut->passes( exeName + ": reset" );
        std::cout.flush();
    }

    // set up variables for apply tests
    // auto diffSolVar = diffOp->getInputVariable(diffOp->getPrincipalVariableId());
    auto diffSolVar             = diffOp->getOutputVariable();
    auto diffRhsVar             = diffOp->getOutputVariable();
    auto diffResVar             = diffOp->getOutputVariable();
    auto workVar                = std::make_shared<AMP::LinearAlgebra::Variable>( "work" );
    auto nonPrincIds            = diffOp->getNonPrincipalVariableIds();
    unsigned int numNonPrincIds = nonPrincIds.size();
    std::vector<AMP::LinearAlgebra::Variable::shared_ptr> nonPrincVars( numNonPrincIds );
    auto inputVar = diffOp->getInputVariable();
    for ( size_t i = 0; i < numNonPrincIds; i++ ) {
        // nonPrincVars[i] = diffOp->getInputVariable(nonPrincIds[i]);
        nonPrincVars[i] = std::dynamic_pointer_cast<AMP::LinearAlgebra::MultiVariable>( inputVar )
                              ->getVariable( i );
    }

    // Test apply
    {
        std::string msgPrefix = exeName + ": apply";
        auto diffSolVec       = AMP::LinearAlgebra::createVector( nodalDofMap, diffSolVar );
        auto diffRhsVec       = AMP::LinearAlgebra::createVector( nodalDofMap, diffRhsVar );
        auto diffResVec       = AMP::LinearAlgebra::createVector( nodalDofMap, diffResVar );
        std::vector<AMP::LinearAlgebra::Vector::shared_ptr> nonPrincVecs( numNonPrincIds );
        for ( unsigned int i = 0; i < numNonPrincIds; i++ ) {
            nonPrincVecs[i] = AMP::LinearAlgebra::createVector( nodalDofMap, nonPrincVars[i] );
            if ( nonPrincIds[i] == AMP::Operator::Diffusion::TEMPERATURE )
                nonPrincVecs[i]->setToScalar( defTemp );
            if ( nonPrincIds[i] == AMP::Operator::Diffusion::CONCENTRATION )
                nonPrincVecs[i]->setToScalar( defConc );
            if ( nonPrincIds[i] == AMP::Operator::Diffusion::BURNUP )
                nonPrincVecs[i]->setToScalar( defBurn );
        }
        diffRhsVec->setToScalar( 0.0 );
        applyTests(
            ut, msgPrefix, nonlinearOperator, diffRhsVec, diffSolVec, diffResVec, shift, scale );
        std::cout.flush();

        // Test getJacobianParameters and linear operator creation
        {
            diffSolVec->setRandomValues();
            adjust( diffSolVec, shift, scale );
            auto jacParams = diffOp->getParameters( "Jacobian", diffSolVec );
            linOp->reset(
                std::dynamic_pointer_cast<AMP::Operator::DiffusionLinearFEOperatorParameters>(
                    jacParams ) );
            ut->passes( exeName + ": getJacobianParameters" );
            std::cout.flush();
        }
    }

    // now run apply tests with multi-vectors
    auto auxInpVar  = diffSolVar->cloneVariable( "NonlinearDiffusionOperator-auxInpVar" );
    auto auxOutVar  = diffResVar->cloneVariable( "NonlinearDiffusionOperator-auxOutVar" );
    auto auxWorkVar = diffSolVar->cloneVariable( "NonlinearDiffusionOperator-auxWorkVar" );

    auto myMultiInpVar =
        std::make_shared<AMP::LinearAlgebra::MultiVariable>( "MultiInputVariable" );
    myMultiInpVar->add( diffSolVar );
    myMultiInpVar->add( auxInpVar );

    auto myMultiOutVar =
        std::make_shared<AMP::LinearAlgebra::MultiVariable>( "MultiOutputVariable" );
    myMultiOutVar->add( diffResVar );
    myMultiOutVar->add( auxOutVar );

    auto myMultiWorkVar =
        std::make_shared<AMP::LinearAlgebra::MultiVariable>( "MultiWorkVariable" );
    myMultiWorkVar->add( workVar );
    myMultiWorkVar->add( auxWorkVar );

    {
        std::string msgPrefix = exeName + ": apply MultiVector ";
        auto solVec           = AMP::LinearAlgebra::createVector( nodalDofMap, myMultiInpVar );
        auto rhsVec           = AMP::LinearAlgebra::createVector( nodalDofMap, myMultiOutVar );
        auto resVec           = AMP::LinearAlgebra::createVector( nodalDofMap, myMultiOutVar );

        // test apply with single variable vectors
        applyTests( ut, msgPrefix, nonlinearOperator, rhsVec, solVec, resVec, shift, scale );
        std::cout.flush();
    }

    // Test isValidInput function
    {
        auto testVec = AMP::LinearAlgebra::createVector( nodalDofMap, diffSolVar );

        testVec->setToScalar( -1000. );
        if ( not diffOp->isValidInput( testVec ) )
            ut->passes( exeName + ": validInput-1" );
        else {
            if ( ( diffOp->getPrincipalVariableId() == AMP::Operator::Diffusion::TEMPERATURE ) &&
                 ( ( mat->property( property ) )->is_argument( "temperature" ) ) ) {
                ut->failure( exeName + ": validInput-1" );
            } else if ( ( diffOp->getPrincipalVariableId() ==
                          AMP::Operator::Diffusion::CONCENTRATION ) &&
                        ( ( mat->property( property ) )->is_argument( "concentration" ) ) ) {
                ut->failure( exeName + ": validInput-1" );
            }
        }
        testVec->setToScalar( 1.e99 );
        if ( not diffOp->isValidInput( testVec ) )
            ut->passes( exeName + ": validInput-2" );
        else {
            if ( ( diffOp->getPrincipalVariableId() == AMP::Operator::Diffusion::TEMPERATURE ) &&
                 ( ( mat->property( property ) )->is_argument( "temperature" ) ) ) {
                ut->failure( exeName + ": validInput-2" );
            } else if ( ( diffOp->getPrincipalVariableId() ==
                          AMP::Operator::Diffusion::CONCENTRATION ) &&
                        ( ( mat->property( property ) )->is_argument( "concentration" ) ) ) {
                ut->failure( exeName + ": validInput-2" );
            }
        }
        testVec->setToScalar( 1.e99 );
        std::cout.flush();
    }
}

int testNonlinearDiffusion_1( int argc, char *argv[] )
{
    AMP::AMPManagerProperties startup_properties;
    startup_properties.use_MPI_Abort = false;
    AMP::AMPManager::startup( argc, argv, startup_properties );

    AMP::UnitTest ut;

    const int NUMFILES          = 14;
    std::string files[NUMFILES] = {
        "Diffusion-CylindricalFick-1",
        "Diffusion-TUI-Thermal-1",
        "Diffusion-TUI-Fick-1",
        "Diffusion-TUI-Soret-1",
        "Diffusion-UO2MSRZC09-Thermal-1",
        "Diffusion-UO2MSRZC09-Fick-1",
        "Diffusion-UO2MSRZC09-Soret-1",
        "Diffusion-TUI-Thermal-ActiveTemperatureAndConcentration-1",
        "Diffusion-TUI-Fick-ActiveTemperatureAndConcentration-1",
        "Diffusion-TUI-Soret-ActiveTemperatureAndConcentration-1",
        "Diffusion-UO2MSRZC09-Thermal-ActiveTemperatureAndConcentration-1",
        "Diffusion-UO2MSRZC09-Fick-ActiveTemperatureAndConcentration-1",
        "Diffusion-UO2MSRZC09-Soret-ActiveTemperatureAndConcentration-1",
        "Diffusion-TUI-TensorFick-1"
    };

    for ( auto &file : files )
        nonlinearTest( &ut, file );

    ut.report();

    int num_failed = ut.NumFailGlobal();
    AMP::AMPManager::shutdown();
    return num_failed;
}
