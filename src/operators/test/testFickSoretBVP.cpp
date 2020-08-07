#include "AMP/ampmesh/Mesh.h"
#include "AMP/discretization/DOF_Manager.h"
#include "AMP/discretization/simpleDOF_Manager.h"
#include "AMP/operators/BVPOperatorParameters.h"
#include "AMP/operators/LinearBVPOperator.h"
#include "AMP/operators/NonlinearBVPOperator.h"
#include "AMP/operators/OperatorBuilder.h"
#include "AMP/operators/boundary/DirichletMatrixCorrection.h"
#include "AMP/operators/boundary/DirichletVectorCorrection.h"
#include "AMP/operators/diffusion/DiffusionLinearFEOperator.h"
#include "AMP/operators/diffusion/FickSoretNonlinearFEOperator.h"
#include "AMP/utils/AMPManager.h"
#include "AMP/utils/AMP_MPI.h"
#include "AMP/utils/Database.h"
#include "AMP/utils/PIO.h"
#include "AMP/utils/UnitTest.h"
#include "AMP/utils/Utilities.h"
#include "AMP/vectors/MultiVariable.h"
#include "AMP/vectors/Variable.h"
#include "AMP/vectors/VectorBuilder.h"
#include <memory>

#include "applyTests.h"

#include <iostream>
#include <string>


static void bvpTest1( AMP::UnitTest *ut, const std::string &exeName )
{
    // Tests FickSoret Dirchlet BVP operator for temperature

    // Initialization
    std::string input_file = "input_" + exeName;
    std::string log_file   = "output_" + exeName;
    AMP::PIO::logOnlyNodeZero( log_file );

    // Input database
    auto input_db = AMP::Database::parseInputFile( input_file );
    input_db->print( AMP::plog );

    // Create the Mesh.
    AMP_INSIST( input_db->keyExists( "Mesh" ), "Key ''Mesh'' is missing!" );
    auto mesh_db   = input_db->getDatabase( "Mesh" );
    auto mgrParams = std::make_shared<AMP::Mesh::MeshParameters>( mesh_db );
    mgrParams->setComm( AMP::AMP_MPI( AMP_COMM_WORLD ) );
    auto meshAdapter = AMP::Mesh::Mesh::buildMesh( mgrParams );

    // Create nonlinear FickSoret BVP operator and access volume nonlinear FickSoret operator
    std::shared_ptr<AMP::Operator::ElementPhysicsModel> elementPhysicsModel;
    auto nlinBVPOperator = AMP::Operator::OperatorBuilder::createOperator(
        meshAdapter, "testFickSoretBVPOperator", input_db, elementPhysicsModel );
    auto nlinBVPOp =
        std::dynamic_pointer_cast<AMP::Operator::NonlinearBVPOperator>( nlinBVPOperator );
    auto nlinOp = std::dynamic_pointer_cast<AMP::Operator::FickSoretNonlinearFEOperator>(
        nlinBVPOp->getVolumeOperator() );
    auto fickOp = std::dynamic_pointer_cast<AMP::Operator::DiffusionNonlinearFEOperator>(
        nlinOp->getFickOperator() );
    auto soretOp = std::dynamic_pointer_cast<AMP::Operator::DiffusionNonlinearFEOperator>(
        nlinOp->getSoretOperator() );

    // use the linear BVP operator to create a Fick linear operator with bc's
    auto linBVPOperator = AMP::Operator::OperatorBuilder::createOperator(
        meshAdapter, "testLinearFickBVPOperator", input_db, elementPhysicsModel );
    auto linBVPOp = std::dynamic_pointer_cast<AMP::Operator::LinearBVPOperator>( linBVPOperator );

    ut->passes( exeName + ": creation" );
    std::cout.flush();

    // set shift, scale for applyTests
    double shift[2];
    double scale[2];
    shift[0] = 0.;
    shift[1] = 0.;
    scale[0] = 1.;
    scale[1] = 1.;
    std::vector<double> trange( 2 ), crange( 2 );
    auto fickTransportModel = fickOp->getTransportModel();
    std::vector<double> defaults( 2, 0 );
    auto fmat = fickTransportModel->getMaterial();
    // the Soret has a principal variable of temperature
    if ( soretOp->getPrincipalVariableId() == AMP::Operator::Diffusion::TEMPERATURE ) {
        std::string property = "ThermalDiffusionCoefficient";
        if ( ( fmat->property( property ) )->is_argument( "temperature" ) ) {
            trange   = ( fmat->property( property ) )->get_arg_range( "temperature" );
            scale[1] = trange[1] - trange[0];
            shift[1] = trange[0] + 0.001 * scale[1];
            scale[1] *= 0.999;
            defaults = ( fmat->property( property ) )->get_defaults();
        }
    }
    // the Fick has a principal variable of temperature
    if ( fickOp->getPrincipalVariableId() == AMP::Operator::Diffusion::CONCENTRATION ) {
        std::string property = "FickCoefficient";
        if ( ( fmat->property( property ) )->is_argument( "concentration" ) ) {
            crange   = ( fmat->property( property ) )->get_arg_range( "concentration" );
            scale[0] = crange[1] - crange[0];
            shift[0] = crange[0] + 0.001 * scale[0];
            scale[0] *= 0.999;
            defaults = ( fmat->property( property ) )->get_defaults();
        }
    }

    // Set up input and output vectors
    auto cVar =
        std::dynamic_pointer_cast<AMP::LinearAlgebra::MultiVariable>( fickOp->getInputVariable() )
            ->getVariable( AMP::Operator::Diffusion::CONCENTRATION );
    auto tVar =
        std::dynamic_pointer_cast<AMP::LinearAlgebra::MultiVariable>( soretOp->getInputVariable() )
            ->getVariable( AMP::Operator::Diffusion::TEMPERATURE );
    auto fsInpVar = std::make_shared<AMP::LinearAlgebra::MultiVariable>( "fsInput" );
    fsInpVar->add( cVar );
    fsInpVar->add( tVar );

    // Create a DOF manager for a nodal vector
    int DOFsPerNode     = 1;
    int nodalGhostWidth = 1;
    bool split          = true;
    auto nodalDofMap    = AMP::Discretization::simpleDOFManager::create(
        meshAdapter, AMP::Mesh::GeomType::Vertex, nodalGhostWidth, DOFsPerNode, split );

    // create solution, rhs, and residual vectors
    auto multivector =
        AMP::LinearAlgebra::MultiVector::create( "mulitvector", meshAdapter->getComm() );
    multivector->addVector( AMP::LinearAlgebra::createVector( nodalDofMap, cVar ) );
    multivector->addVector( AMP::LinearAlgebra::createVector( nodalDofMap, tVar ) );
    auto solVec = std::dynamic_pointer_cast<AMP::LinearAlgebra::Vector>( multivector );
    auto rhsVec = solVec->cloneVector();
    auto resVec = solVec->cloneVector();

    // set default values of input variables
    auto inConcVec = solVec->subsetVectorForVariable( cVar );
    auto inTempVec = solVec->subsetVectorForVariable( tVar );
    inConcVec->setToScalar( defaults[1], inConcVec ); // compile error
    inTempVec->setToScalar( defaults[0], inTempVec ); // compile error
    rhsVec->setToScalar( 0., rhsVec );

    AMP_INSIST( nlinOp->isValidInput( solVec ), "input variable not set up correctly" );

    // Test apply
    std::string msgPrefix = exeName + ": apply";
    {
        applyTests( ut, msgPrefix, nlinBVPOperator, rhsVec, solVec, resVec, shift, scale, 2 );
    }
    std::cout.flush();

    // Test linear reset from getJacobianParameters
    for ( int i = 0; i < 3; i++ ) {
        inConcVec->setRandomValues(inConcVec);
        inTempVec->setRandomValues(inTempVec);
        adjust( solVec, shift, scale, 2 );
        rhsVec->setRandomValues(rhsVec);
        resVec->setRandomValues(resVec);
        std::shared_ptr<AMP::Operator::OperatorParameters> jacparams =
            nlinBVPOp->getParameters( "Jacobian", solVec );
        linBVPOp->reset( jacparams );
    } // end for i

    ut->passes( exeName + ": getJacobianParameters" );
    std::cout.flush();
}


int testFickSoretBVP( int argc, char *argv[] )
{
    AMP::AMPManagerProperties startup_properties;
    startup_properties.use_MPI_Abort = false;
    AMP::AMPManager::startup( argc, argv, startup_properties );

    AMP::UnitTest ut;

    std::vector<std::string> files = { "FickSoret-BVP-TUI-1", "FickSoret-BVP-UO2MSRZC09-1" };
    for ( auto &file : files )
        bvpTest1( &ut, file );

    ut.report();

    int num_failed = ut.NumFailGlobal();
    AMP::AMPManager::shutdown();
    return num_failed;
}
