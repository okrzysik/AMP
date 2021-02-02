#include "AMP/ampmesh/Mesh.h"
#include "AMP/ampmesh/MeshParameters.h"
#include "AMP/discretization/DOF_Manager.h"
#include "AMP/discretization/simpleDOF_Manager.h"
#include "AMP/operators/BVPOperatorParameters.h"
#include "AMP/operators/LinearBVPOperator.h"
#include "AMP/operators/NonlinearBVPOperator.h"
#include "AMP/operators/OperatorBuilder.h"
#include "AMP/operators/boundary/DirichletMatrixCorrection.h"
#include "AMP/operators/boundary/DirichletVectorCorrection.h"
#include "AMP/operators/diffusion/DiffusionLinearElement.h"
#include "AMP/operators/diffusion/DiffusionLinearFEOperator.h"
#include "AMP/operators/diffusion/DiffusionNonlinearElement.h"
#include "AMP/operators/diffusion/DiffusionNonlinearFEOperator.h"
#include "AMP/utils/AMPManager.h"
#include "AMP/utils/AMP_MPI.h"
#include "AMP/utils/Database.h"
#include "AMP/utils/PIO.h"
#include "AMP/utils/UnitTest.h"
#include "AMP/utils/Utilities.h"
#include "AMP/vectors/VectorBuilder.h"

#include "../applyTests.h"

#include <iostream>
#include <memory>
#include <string>


static void bvpTest1( AMP::UnitTest *ut, const std::string &exeName )
{
    // Tests diffusion Dirchlet BVP operator for temperature

    // Initialization
    std::string input_file = "input_" + exeName;
    std::string log_file   = "output_" + exeName;
    AMP::PIO::logOnlyNodeZero( log_file );

    // Input database
    auto input_db = AMP::Database::parseInputFile( input_file );
    input_db->print( AMP::plog );

    //   Create the Mesh.
    AMP_INSIST( input_db->keyExists( "Mesh" ), "Key ''Mesh'' is missing!" );
    auto mesh_db   = input_db->getDatabase( "Mesh" );
    auto mgrParams = std::make_shared<AMP::Mesh::MeshParameters>( mesh_db );
    mgrParams->setComm( AMP::AMP_MPI( AMP_COMM_WORLD ) );
    auto meshAdapter = AMP::Mesh::Mesh::buildMesh( mgrParams );

    // Create nonlinear Diffusion BVP operator and access volume nonlinear Diffusion operator
    auto nbvp_db = std::dynamic_pointer_cast<AMP::Database>(
        input_db->getDatabase( "ThermalNonlinearBVPOperator" ) );
    auto nlinBVPOperator = AMP::Operator::OperatorBuilder::createOperator(
        meshAdapter, "ThermalNonlinearBVPOperator", input_db );
    auto nlinBVPOp =
        std::dynamic_pointer_cast<AMP::Operator::NonlinearBVPOperator>( nlinBVPOperator );
    auto nlinOp = std::dynamic_pointer_cast<AMP::Operator::DiffusionNonlinearFEOperator>(
        nlinBVPOp->getVolumeOperator() );
    auto elementPhysicsModel = nlinOp->getTransportModel();

    // use the linear BVP operator to create a thermal linear operator with bc's
    auto linBVPOperator = AMP::Operator::OperatorBuilder::createOperator(
        meshAdapter, "ThermalLinearBVPOperator", input_db, elementPhysicsModel );
    auto linBVPOp = std::dynamic_pointer_cast<AMP::Operator::LinearBVPOperator>( linBVPOperator );
    ut->passes( exeName + ": creation" );
    std::cout.flush();

    // Set up input and output vectors
    auto bvpSolVar = nlinOp->getOutputVariable();
    auto bvpRhsVar = nlinOp->getOutputVariable();
    auto bvpResVar = nlinOp->getOutputVariable();

    // Create a DOF manager for a nodal vector
    int DOFsPerNode     = 1;
    int nodalGhostWidth = 1;
    bool split          = true;
    auto nodalDofMap    = AMP::Discretization::simpleDOFManager::create(
        meshAdapter, AMP::Mesh::GeomType::Vertex, nodalGhostWidth, DOFsPerNode, split );

    // create solution, rhs, and residual vectors
    auto bvpSolVec = AMP::LinearAlgebra::createVector( nodalDofMap, bvpSolVar );
    auto bvpRhsVec = AMP::LinearAlgebra::createVector( nodalDofMap, bvpRhsVar );
    auto bvpResVec = AMP::LinearAlgebra::createVector( nodalDofMap, bvpResVar );

    bvpRhsVec->setToScalar( 0.0 );

    auto volOp_db = input_db->getDatabase( nbvp_db->getString( "VolumeOperator" ) );
    auto model_db = input_db->getDatabase( volOp_db->getString( "LocalModel" ) );
    auto property = model_db->getString( "Property" );

    // set shift, scale for applyTests
    double shift = 0., scale = 1.;
    auto transportModel =
        std::dynamic_pointer_cast<AMP::Operator::DiffusionTransportModel>( elementPhysicsModel );
    auto mat = transportModel->getMaterial();
    if ( nlinOp->getPrincipalVariableId() == AMP::Operator::Diffusion::TEMPERATURE ) {
        if ( ( mat->property( property ) )->is_argument( "temperature" ) ) {
            auto range =
                ( mat->property( property ) )->get_arg_range( "temperature" ); // Compile error
            scale = range[1] - range[0];
            shift = range[0] + 0.001 * scale;
            scale *= 0.999;
        }
    }
    if ( nlinOp->getPrincipalVariableId() == AMP::Operator::Diffusion::CONCENTRATION ) {
        if ( ( mat->property( property ) )->is_argument( "concentration" ) ) {
            auto range =
                ( mat->property( property ) )->get_arg_range( "concentration" ); // Compile error
            scale = range[1] - range[0];
            shift = range[0] + 0.001 * scale;
            scale *= 0.999;
        }
    }
    if ( nlinOp->getPrincipalVariableId() == AMP::Operator::Diffusion::BURNUP ) {
        AMP_INSIST( false, "do not know what to do" );
    }

    // Test apply
    std::string msgPrefix = exeName + ": apply";
    {
        applyTests( ut, msgPrefix, nlinBVPOperator, bvpRhsVec, bvpSolVec, bvpResVec, shift, scale );
    }
    std::cout.flush();

    // Test linear reset from getJacobianParameters
    for ( int i = 0; i < 3; i++ ) {
        bvpSolVec->setRandomValues();
        adjust( bvpSolVec, shift, scale );
        bvpRhsVec->setRandomValues();
        bvpResVec->setRandomValues();
        std::shared_ptr<AMP::Operator::OperatorParameters> jacparams =
            nlinBVPOp->getParameters( "Jacobian", bvpSolVec );
        linBVPOp->reset( jacparams );
    } // end for i

    ut->passes( exeName + ": getJacobianParameters" );
    std::cout.flush();
}

int testDiffusionBVP_1( int argc, char *argv[] )
{
    AMP::AMPManagerProperties startup_properties;
    startup_properties.use_MPI_Abort = false;
    AMP::AMPManager::startup( argc, argv, startup_properties );

    AMP::UnitTest ut;

    const int NUMFILES          = 6;
    std::string files[NUMFILES] = {
        "Diffusion-TUI-Thermal-2",        "Diffusion-TUI-Fick-2",        "Diffusion-TUI-Soret-2",
        "Diffusion-UO2MSRZC09-Thermal-2", "Diffusion-UO2MSRZC09-Fick-2",
        "Diffusion-UO2MSRZC09-Soret-2" //,
        //"Diffusion-TUI-Fick-3", "Diffusion-UO2MSRZC09-Fick-3"
    };

    for ( auto &file : files )
        bvpTest1( &ut, file );

    ut.report();

    int num_failed = ut.NumFailGlobal();
    AMP::AMPManager::shutdown();
    return num_failed;
}
