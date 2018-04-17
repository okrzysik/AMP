#include "AMP/utils/AMPManager.h"
#include "AMP/utils/UnitTest.h"
#include "AMP/utils/Utilities.h"
#include <string>

#include "AMP/materials/Material.h"
#include "AMP/utils/AMPManager.h"
#include "AMP/utils/Database.h"
#include "AMP/utils/InputDatabase.h"
#include "AMP/utils/InputManager.h"
#include "AMP/utils/PIO.h"
#include "AMP/utils/Utilities.h"
#include "AMP/utils/shared_ptr.h"


#include "AMP/vectors/SimpleVector.h"
#include "AMP/vectors/Variable.h"
#include "AMP/vectors/Vector.h"

#include "AMP/ampmesh/MeshVariable.h"
#include "AMP/utils/Writer.h"

#include "AMP/operators/ColumnBoundaryOperator.h"
#include "AMP/operators/ElementOperationFactory.h"
#include "AMP/operators/ElementPhysicsModelFactory.h"
#include "AMP/operators/diffusion/DiffusionLinearElement.h"
#include "AMP/operators/diffusion/DiffusionLinearFEOperator.h"
#include "AMP/operators/diffusion/DiffusionNonlinearFEOperator.h"
#include "AMP/operators/diffusion/DiffusionTransportModel.h"
#include "AMP/operators/libmesh/MassLinearElement.h"
#include "AMP/operators/libmesh/MassLinearFEOperator.h"

#include "AMP/operators/CoupledOperator.h"
#include "AMP/operators/CoupledOperatorParameters.h"
#include "AMP/operators/DirichletMatrixCorrection.h"
#include "AMP/operators/DirichletVectorCorrection.h"
#include "AMP/operators/LinearBVPOperator.h"
#include "AMP/operators/NeumannVectorCorrection.h"
#include "AMP/operators/NeumannVectorCorrectionParameters.h"
#include "AMP/operators/NeutronicsRhs.h"
#include "AMP/operators/NonlinearBVPOperator.h"
#include "AMP/operators/OperatorBuilder.h"
#include "AMP/operators/RobinMatrixCorrection.h"
#include "AMP/operators/RobinVectorCorrection.h"
#include "AMP/operators/libmesh/VolumeIntegralOperator.h"
#include "AMP/operators/map/MapOperatorParameters.h"
#include "AMP/operators/map/MapSurface.h"


void resetTests( AMP::UnitTest *ut,
                 std::string msgPrefix,
                 AMP::shared_ptr<AMP::Mesh::MeshManager::Adapter>,
                 AMP::shared_ptr<AMP::Operator::Operator>,
                 AMP::shared_ptr<AMP::InputDatabase> )
//             AMP::shared_ptr<AMP::Operator::OperatorParameters> &bcParameters)
{

    ut.passes( msgPrefix );
    std::cout.flush();
}

void adjust( const AMP::LinearAlgebra::Vector::shared_ptr vec,
             AMP::LinearAlgebra::Vector::shared_ptr work )
{
    work->setToScalar( 301. );
    AMP::LinearAlgebra::Vector &x = *vec;
    AMP::LinearAlgebra::Vector &y = *work;
    vec->add( x, y );
}

void applyTest( AMP::UnitTest *ut,
                std::string msgPrefix,
                AMP::shared_ptr<AMP::Operator::Operator> &testOperator,
                AMP::LinearAlgebra::Vector::shared_ptr rhsVec,
                AMP::LinearAlgebra::Vector::shared_ptr solVec,
                AMP::LinearAlgebra::Vector::shared_ptr resVec,
                AMP::LinearAlgebra::Vector::shared_ptr workVec )
{
    // first test for apply - random values in all three input vectors
    bool passed = true;
    try {
        for ( int j = 0; j < 3; j++ ) {
            solVec->setRandomValues();
            rhsVec->setRandomValues();
            resVec->setRandomValues();
            adjust( solVec, workVec );
            testOperator->apply( rhsVec, solVec, resVec, 1.0, -1.0 );
        } // end for j
    } catch ( std::exception ) {
        passed = false;
    }
    if ( passed ) {
        ut.passes( msgPrefix + " : apply with random f, u, r, a=1, b=-1.0" );
    } else {
        ut.numFails++;
        ut.failure( msgPrefix + " : apply with random f, u, r, a=1, b=-1.0" );
    }

    // second test for apply - f NULL, u, r, random values
    passed = true;
    try {
        for ( int j = 0; j < 3; j++ ) {
            AMP::LinearAlgebra::Vector::shared_ptr fVec;
            solVec->setRandomValues();
            resVec->setRandomValues();
            adjust( solVec, workVec );
            testOperator->apply( fVec, solVec, resVec, 1.0, -1.0 );
        } // end for j
    } catch ( std::exception ) {
        passed = false;
    }
    if ( passed ) {
        ut.passes( msgPrefix + " : apply with f NULL, random u, r, a=1, b=-1.0" );
    } else {
        ut.numFails++;
        ut.failure( msgPrefix + " : apply with f NULL, random u, r, a=1, b=-1.0" );
    }

    // R.S.: u is allowed to be NULL for some operators. For example, operators
    // with an in-place apply. However, this test is not meant to be used with those operators.
    // third test for apply - u NULL, f, r, random values
    passed = false;
    try {
        for ( int j = 0; j < 3; j++ ) {
            AMP::LinearAlgebra::Vector::shared_ptr uVec;
            rhsVec->setRandomValues();
            resVec->setRandomValues();
            testOperator->apply( rhsVec, uVec, resVec, 1.0, -1.0 );
        } // end for j
    } catch ( std::exception ) {
        passed = true;
    }
    if ( passed ) {
        ut.passes( msgPrefix +
                   " : apply with u NULL, random values in the vectors f,r, a=1, b=-1.0" );
    } else {
        ut.numFails++;
        ut.failure( msgPrefix +
                    " : apply with u NULL, random values in the vectors f,r, a=1, b=-1.0" );
    }

    // fourth test for apply - r NULL, f, u, random values
    passed = false;
    try {
        for ( int j = 0; j < 3; j++ ) {
            AMP::LinearAlgebra::Vector::shared_ptr rVec;
            solVec->setRandomValues();
            rhsVec->setRandomValues();
            adjust( solVec, workVec );
            testOperator->apply( rhsVec, solVec, rVec, 1.0, -1.0 );
        } // end for j
    } catch ( std::exception ) {
        passed = true;
    }
    if ( passed ) {
        ut.passes( msgPrefix +
                   " : apply with r NULL, random values in the vectors f,u, a=1, b=-1.0" );
    } else {
        ut.numFails++;
        ut.failure( msgPrefix +
                    " : apply with r NULL, random values in the vectors f,u, a=1, b=-1.0" );
    }

    // fifth test for apply - f NULL, u NULL, r, random values
    passed = false;
    try {
        for ( int j = 0; j < 3; j++ ) {
            AMP::LinearAlgebra::Vector::shared_ptr fVec;
            AMP::LinearAlgebra::Vector::shared_ptr uVec;
            resVec->setRandomValues();
            testOperator->apply( fVec, uVec, resVec, 1.0, -1.0 );
        } // end for j
    } catch ( std::exception ) {
        passed = true;
    }
    if ( passed ) {
        ut.passes( msgPrefix +
                   " : apply with f NULL, u NULL random values in the vector r, a=1, b=-1.0" );
    } else {
        ut.numFails++;
        ut.failure( msgPrefix +
                    " : apply with f NULL, u NULL random values in the vector r, a=1, b=-1.0" );
    }

    // sixth test for apply - u NULL, r NULL, f, random values
    passed = false;
    try {
        for ( int j = 0; j < 3; j++ ) {
            AMP::LinearAlgebra::Vector::shared_ptr uVec;
            AMP::LinearAlgebra::Vector::shared_ptr rVec;
            rhsVec->setRandomValues();
            testOperator->apply( rhsVec, uVec, rVec, 1.0, -1.0 );
        } // end for j
    } catch ( std::exception ) {
        passed = true;
    }
    if ( passed ) {
        ut.passes( msgPrefix +
                   " : apply with u NULL, r NULL, random values in the vector f, a=1, b=-1.0" );
    } else {
        ut.numFails++;
        ut.failure( msgPrefix +
                    " : apply with u NULL, r NULL, random values in the vector f, a=1, b=-1.0" );
    }

    // seventh test for apply - r NULL, f NULL, u random values
    passed = false;
    try {
        for ( int j = 0; j < 3; j++ ) {
            AMP::LinearAlgebra::Vector::shared_ptr rVec;
            AMP::LinearAlgebra::Vector::shared_ptr fVec;
            solVec->setRandomValues();
            adjust( solVec, workVec );
            testOperator->apply( fVec, solVec, rVec, 1.0, -1.0 );
        } // end for j
    } catch ( std::exception ) {
        passed = true;
    }
    if ( passed ) {
        ut.passes( msgPrefix +
                   " : apply with f, r NULL, random values in the vector u, a=1, b=-1.0" );
    } else {
        ut.numFails++;
        ut.failure( msgPrefix +
                    " : apply with f, r NULL, random values in the vector u, a=1, b=-1.0" );
    }

    // eighth test for apply - r NULL, f NULL, u NULL
    passed = false;
    try {
        for ( int j = 0; j < 3; j++ ) {
            AMP::LinearAlgebra::Vector::shared_ptr rVec;
            AMP::LinearAlgebra::Vector::shared_ptr fVec;
            AMP::LinearAlgebra::Vector::shared_ptr uVec;
            testOperator->apply( fVec, uVec, rVec, 1.0, -1.0 );
        } // end for j
    } catch ( std::exception ) {
        passed = true;
    }
    if ( passed ) {
        ut.passes( msgPrefix + " : apply with f, u, r NULL, a=1, b=-1.0" );
    } else {
        ut.numFails++;
        ut.failure( msgPrefix + " : apply with f, u, r NULL, a=1, b=-1.0" );
    }

    // ninth test for apply - random values in all three input vectors, a=0, b=1
    rhsVec->setRandomValues();
    resVec->setRandomValues();
    solVec->setRandomValues();
    adjust( solVec, workVec );
    testOperator->apply( rhsVec, solVec, resVec, 0.0, 1.0 );
    rhsVec->subtract( rhsVec, resVec );
    if ( AMP::Utilities::approx_equal( rhsVec->L2Norm(), 0.0 ) ) {
        ut.passes( msgPrefix + " : apply with random values in the vectors f,u,r, a=0.0, b=1.0" );
    } else {
        ut.numFails++;
        ut.failure( msgPrefix + " : apply with random values in the vectors f,u,r, a=0.0, b=1.0" );
    }

    // tenth test for apply - random values in all three input vectors, a=0, b=-1, to test scaling
    rhsVec->setRandomValues();
    resVec->setRandomValues();
    solVec->setRandomValues();
    adjust( solVec, workVec );
    testOperator->apply( rhsVec, solVec, resVec, 0.0, -1.0 );
    rhsVec->add( rhsVec, resVec );
    if ( AMP::Utilities::approx_equal( rhsVec->L2Norm(), 0.0 ) ) {
        ut.passes(
            msgPrefix +
            " : apply with random values in the vectors f,u,r, a=0.0, b=-1.0 (test scaling of f)" );
    } else {
        ut.numFails++;
        ut.failure(
            msgPrefix +
            " : apply with random values in the vectors f,u,r, a=0.0, b=-1.0 (test scaling of f)" );
    }
}


void thermalContactApplyTest( AMP::UnitTest *ut, const std::string &exeName )
{
    std::string input_file = "input_" + exeName;
    std::string log_file   = "output_" + exeName;

    //  AMP::AMPManager::startup();
    //

    auto input_db = AMP::make_shared<AMP::InputDatabase>( "input_db" );
    AMP::InputManager::getManager()->parseInputFile( input_file, input_db );
    input_db->printClassData( AMP::plog );

    AMP::PIO::logAllNodes( log_file );

    //  AMP_INSIST(input_db->keyExists("Mesh"), "Key ''Mesh'' is missing!");
    //  std::string mesh_file = input_db->getString("Mesh");

    auto mgrParams    = AMP::make_shared<AMP::Mesh::MeshManagerParameters>( input_db );
    auto manager      = AMP::make_shared<AMP::Mesh::MeshManager>( mgrParams );
    auto meshAdapter1 = manager->getMesh( "pellet" );
    auto meshAdapter2 = manager->getMesh( "clad" );

    AMP::LinearAlgebra::Vector::shared_ptr nullVec;
    auto TemperatureVar  = AMP::make_shared<AMP::Mesh::NodalScalarVariable>( "Temperature" );
    auto inputVariable1  = AMP::make_shared<AMP::Mesh::NodalScalarVariable>( "Temperature", meshAdapter1 );
    auto inputVariable2  = AMP::make_shared<AMP::Mesh::NodalScalarVariable>( "Temperature", meshAdapter2 );
    auto outputVariable1 = AMP::make_shared<AMP::Mesh::NodalScalarVariable>( "Temperature", meshAdapter1 );
    auto outputVariable2 = AMP::make_shared<AMP::Mesh::NodalScalarVariable>( "Temperature", meshAdapter2 );

    double intguess = input_db->getDoubleWithDefault( "InitialGuess", 400 );

    auto TemperatureInKelvin = manager->createVector( TemperatureVar );
    auto RightHandSideVec    = manager->createVector( TemperatureVar );
    auto ResidualVec         = manager->createVector( TemperatureVar );
    auto WorkVec             = manager->createVector( TemperatureVar );

    TemperatureInKelvin->setToScalar( intguess );
    manager->registerVectorAsData( TemperatureInKelvin, "Temperature" );

    //-----------------------------------------------
    //   CREATE THE NONLINEAR THERMAL OPERATOR 1 ----
    //-----------------------------------------------

    AMP_INSIST( input_db->keyExists( "NonlinearThermalOperator1" ), "key missing!" );

    auto thermalTransportModel1;
    auto nonlinearThermalDatabase1 = input_db->getDatabase( "NonlinearThermalOperator1" );
    auto nonlinearThermalOperator1 = AMP::dynamic_pointer_cast<AMP::Operator::NonlinearBVPOperator>(
        AMP::Operator::OperatorBuilder::createOperator(
            meshAdapter1, nonlinearThermalDatabase1, thermalTransportModel1 ) );

    // initialize the input variable
    auto thermalVolumeOperator1 = AMP::dynamic_pointer_cast<AMP::Operator::DiffusionNonlinearFEOperator>(
        nonlinearThermalOperator1->getVolumeOperator() );

    // initialize the output variable
    // auto outputVariable1 = thermalVolumeOperator1->getOutputVariable();

    auto TemperatureInKelvinVec1 = TemperatureInKelvin->subsetVectorForVariable( inputVariable1 );
    auto RightHandSideVec1 = RightHandSideVec->subsetVectorForVariable( outputVariable1 );
    auto ResidualVec1 = ResidualVec->subsetVectorForVariable( outputVariable1 );

    //-------------------------------------
    //   CREATE THE LINEAR THERMAL OPERATOR 1 ----
    //-------------------------------------

    AMP::shared_ptr<AMP::Operator::ElementPhysicsModel> transportModel1;
    auto bvpDatabase1 = AMP::dynamic_pointer_cast<AMP::InputDatabase>(
        input_db->getDatabase( "LinearThermalOperator1" ) );
    auto linearThermalOperator1 = AMP::dynamic_pointer_cast<AMP::Operator::LinearBVPOperator>(
        AMP::Operator::OperatorBuilder::createOperator(
            meshAdapter1, bvpDatabase1, transportModel1 ) );

    //-------------------------------------
    //  CREATE THE NEUTRONICS SOURCE  //
    //-------------------------------------
    AMP_INSIST( input_db->keyExists( "NeutronicsOperator" ),
                "Key ''NeutronicsOperator'' is missing!" );
    auto neutronicsOp_db = input_db->getDatabase( "NeutronicsOperator" );
    auto neutronicsParams = AMP::make_shared<AMP::Operator::NeutronicsRhsParameters>( neutronicsOp_db );
    neutronicsParams->d_MeshAdapter = meshAdapter1;
    auto neutronicsOperator = AMP::make_shared<AMP::Operator::NeutronicsRhs>( neutronicsParams );

    auto SpecificPowerVar = neutronicsOperator->getOutputVariable();
    auto SpecificPowerVec = meshAdapter1->createVector( SpecificPowerVar );

    neutronicsOperator->apply( nullVec, nullVec, SpecificPowerVec, 1., 0. );

    //----------------------------------------------------------
    //  Integrate Nuclear Rhs over Desnity * GeomType::Volume //
    //----------------------------------------------------------

    AMP_INSIST( input_db->keyExists( "VolumeIntegralOperator" ), "key missing!" );

    AMP::shared_ptr<AMP::Operator::ElementPhysicsModel> stransportModel;
    auto sourceDatabase = input_db->getDatabase( "VolumeIntegralOperator" );
    auto sourceOperator = AMP::dynamic_pointer_cast<AMP::Operator::VolumeIntegralOperator>(
        AMP::Operator::OperatorBuilder::createOperator(
            meshAdapter1, sourceDatabase, stransportModel ) );

    // Create the power (heat source) vector.
    auto PowerInWattsVar = sourceOperator->getOutputVariable();
    auto PowerInWattsVec = meshAdapter1->createVector( PowerInWattsVar );
    PowerInWattsVec->zero();

    // convert the vector of specific power to power for a given basis.
    sourceOperator->apply( nullVec, SpecificPowerVec, PowerInWattsVec, 1., 0. );

    // copy the power to pellet RHS vector
    RightHandSideVec1->copyVector( PowerInWattsVec );

    //--------------------------------------------
    //   CREATE THE NONLINEAR THERMAL OPERATOR 2 ----
    //--------------------------------------------

    AMP_INSIST( input_db->keyExists( "NonlinearThermalOperator2" ), "key missing!" );

    AMP::shared_ptr<AMP::Operator::ElementPhysicsModel> thermalTransportModel2;
    auto nonlinearThermalDatabase2 = input_db->getDatabase( "NonlinearThermalOperator2" );
    auto nonlinearThermalOperator2 = AMP::dynamic_pointer_cast<AMP::Operator::NonlinearBVPOperator>(
        AMP::Operator::OperatorBuilder::createOperator(
            meshAdapter2, nonlinearThermalDatabase2, thermalTransportModel2 ) );

    //----------------------------------------------------------------------------------------------------------------------------------------------//
    auto thermalVolumeOperator2 = AMP::dynamic_pointer_cast<AMP::Operator::DiffusionNonlinearFEOperator>(
        nonlinearThermalOperator2->getVolumeOperator() );

    // initialize the output variable
    // AMP::LinearAlgebra::Variable::shared_ptr outputVariable2 =
    // thermalVolumeOperator2->getOutputVariable();

    auto TemperatureInKelvinVec2 = TemperatureInKelvin->subsetVectorForVariable( inputVariable2 );
    auto RightHandSideVec2 = RightHandSideVec->subsetVectorForVariable( outputVariable2 );
    auto ResidualVec2 = ResidualVec->subsetVectorForVariable( outputVariable2 );

    //--------------------------------------------
    //   CREATE THE LINEAR THERMAL OPERATOR 2 ----
    //--------------------------------------------

    AMP::shared_ptr<AMP::Operator::ElementPhysicsModel> transportModel2;
    auto bvpDatabase2 = AMP::dynamic_pointer_cast<AMP::InputDatabase>(
        input_db->getDatabase( "LinearThermalOperator2" ) );
    auto linearThermalOperator2 = AMP::dynamic_pointer_cast<AMP::Operator::LinearBVPOperator>(
        AMP::Operator::OperatorBuilder::createOperator(
            meshAdapter2, bvpDatabase2, transportModel2 ) );

    //-------------------------------------
    auto mapcladtopellet_db = AMP::dynamic_pointer_cast<AMP::InputDatabase>(
        input_db->getDatabase( "MapCladtoPellet" ) );
    auto mapcladtopellet = AMP::dynamic_pointer_cast<AMP::Operator::MapSurface>(
        AMP::Operator::OperatorBuilder::createOperator(
            meshAdapter2, meshAdapter1, mapcladtopellet_db ) );

    auto mappellettoclad_db = AMP::dynamic_pointer_cast<AMP::InputDatabase>(
        input_db->getDatabase( "MapPellettoClad" ) );
    auto mappellettoclad = AMP::dynamic_pointer_cast<AMP::Operator::MapSurface>(
        AMP::Operator::OperatorBuilder::createOperator(
            meshAdapter1, meshAdapter2, mappellettoclad_db ) );

    //------------------------------------------

    AMP::Operator::Operator::shared_ptr boundaryOp1;
    boundaryOp1 = nonlinearThermalOperator1->getBoundaryOperator();

    AMP::Operator::Operator::shared_ptr robinBoundaryOp1;
    robinBoundaryOp1 =
        ( AMP::dynamic_pointer_cast<AMP::Operator::BoundaryOperator>( boundaryOp1 ) );

    AMP::shared_ptr<AMP::InputDatabase> boundaryDatabase1 =
        AMP::dynamic_pointer_cast<AMP::InputDatabase>(
            nonlinearThermalDatabase1->getDatabase( "BoundaryOperator" ) );
    AMP::shared_ptr<AMP::InputDatabase> robinboundaryDatabase1 =
        AMP::dynamic_pointer_cast<AMP::InputDatabase>( boundaryDatabase1 );

    robinboundaryDatabase1->putBool( "constant_flux", false );
    AMP::shared_ptr<AMP::Operator::NeumannVectorCorrectionParameters> correctionParameters1 = AMP::make_shared<AMP::Operator::NeumannVectorCorrectionParameters>( robinboundaryDatabase1 );


    //------------------------------------------

    AMP::Operator::Operator::shared_ptr boundaryOp2;
    boundaryOp2 = nonlinearThermalOperator2->getBoundaryOperator();

    AMP::Operator::Operator::shared_ptr robinBoundaryOp2;
    robinBoundaryOp2 =
        ( AMP::dynamic_pointer_cast<AMP::Operator::ColumnBoundaryOperator>( boundaryOp2 ) )
            ->getBoundaryOperator( 0 );

    auto boundaryDatabase2 = AMP::dynamic_pointer_cast<AMP::InputDatabase>(
        nonlinearThermalDatabase2->getDatabase( "BoundaryOperator" ) );
    auto robinboundaryDatabase2 = AMP::dynamic_pointer_cast<AMP::InputDatabase>(
        boundaryDatabase2->getDatabase( "RobinVectorCorrection" ) );

    robinboundaryDatabase2->putBool( "constant_flux", false );
    auto correctionParameters2 = AMP::make_shared<AMP::Operator::NeumannVectorCorrectionParameters>( robinboundaryDatabase2 );


    //-------------------------------------
    // Applying boundary conditions to the nonlinear BVP Operator

    nonlinearThermalOperator1->modifyRHSvector( RightHandSideVec1 );
    nonlinearThermalOperator1->modifyInitialSolutionVector( TemperatureInKelvinVec1 );

    nonlinearThermalOperator2->modifyRHSvector( RightHandSideVec2 );
    nonlinearThermalOperator2->modifyInitialSolutionVector( TemperatureInKelvinVec2 );

    //-------------------------------------
    // Coupling Map to the Nonlinear Operators
    auto tmp_db = AMP::make_shared<AMP::InputDatabase>( "Dummy" );
    auto coupledNonlinearPelletParams = AMP::make_shared<AMP::Operator::CoupledOperatorParameter>s( tmp_db );

    coupledNonlinearPelletParams->d_MapOperator = mapcladtopellet;
    coupledNonlinearPelletParams->d_BVPOperator = nonlinearThermalOperator1;
    auto coupledNonlinearPellet = AMP::make_shared<AMP::Operator::CoupledOperator>( coupledNonlinearPelletParams );
    //-------------------------------------
    auto coupledNonlinearCladParams = AMP::make_shared<AMP::Operator::CoupledOperatorParameters>( tmp_db );
    coupledNonlinearCladParams->d_MapOperator = mappellettoclad;
    coupledNonlinearCladParams->d_BVPOperator = nonlinearThermalOperator2;
    auto coupledNonlinearClad = AMP::make_shared<AMP::Operator::CoupledOperator>( coupledNonlinearCladParams );

    //-------------------------------------
    // Column of Coupled Operators
    auto nonlinearParams = AMP::make_shared<AMP::Operator::OperatorParameters( tmp_db );
    auto nonlinearCoupledOperator = AMP::make_shared<AMP::Operator::ColumnOperator( nonlinearParams ) );
    nonlinearCoupledOperator->append( coupledNonlinearPellet );
    nonlinearCoupledOperator->append( coupledNonlinearClad );

    //-------------------------------------
    // Coupling Map to the Linear Operators
    auto coupledLinearPelletParams = AMP::make_shared<AMP::Operator::CoupledOperatorParameters>( tmp_db );
    coupledLinearPelletParams->d_MapOperator = mapcladtopellet;
    coupledLinearPelletParams->d_BVPOperator = linearThermalOperator1;
    auto coupledLinearPellet = AMP::make_shared<AMP::Operator::CoupledOperator>( coupledLinearPelletParams );
    //-------------------------------------
    autocoupledLinearCladParams = AMP::make_shared<AMP::Operator::CoupledOperatorParameters>( tmp_db );
    coupledLinearCladParams->d_MapOperator = mappellettoclad;
    coupledLinearCladParams->d_BVPOperator = linearThermalOperator2;
    auto coupledLinearClad = AMP::make_shared<AMP::Operator::CoupledOperator>( coupledLinearCladParams );

    //-------------------------------------
    // Column of Coupled Operators
    auto linearParams = AMP::make_shared<AMP::Operator::OperatorParameters>( tmp_db );
    auto linearCoupledOperator = AMP::make_shared<AMP::Operator::ColumnOperator>( linearParams );
    linearCoupledOperator->append( coupledLinearPellet );
    linearCoupledOperator->append( coupledLinearClad );

    //-------------------------------------


    nonlinearCoupledOperator->apply(
        RightHandSideVec, TemperatureInKelvin, ResidualVec, 1.0, -1.0 );

    correctionParameters1->d_variableFlux = mapcladtopellet->getBoundaryVector();
    correctionParameters2->d_variableFlux = mappellettoclad->getBoundaryVector();

    robinBoundaryOp1->reset( correctionParameters1 );
    robinBoundaryOp2->reset( correctionParameters2 );

    ut.passes( exeName + " : create" );

    // test apply
    std::string msgPrefix = exeName + " : apply";
    auto testOperator = AMP::dynamic_pointer_cast<AMP::Operator::Operator>( nonlinearCoupledOperator );
    applyTest(
        ut, msgPrefix, testOperator, RightHandSideVec, TemperatureInKelvin, ResidualVec, WorkVec );
    ut.passes( msgPrefix );

    // test getJacobian
    AMP::shared_ptr<AMP::Operator::OperatorParameters> resetParams =
        nonlinearCoupledOperator->getParameters( "Jacobian", TemperatureInKelvin );
    ut.passes( exeName + " : getJacobianParameters" );

    // test reset
    linearCoupledOperator->reset( resetParams );
    ut.passes( exeName + " : Linear::reset" );

    input_db.reset();

    ut.passes( exeName );

    //  AMP::AMPManager::shutdown();
}

int main( int argc, char *argv[] )
{
    AMP::AMPManager::startup( argc, argv );
    AMP::UnitTest ut;

    thermalContactApplyTest( ut, "testNonlinearThermalCoupledContact-Apply" );

    ut.report();

    int num_failed = ut.NumFailGlobal();
    AMP::AMPManager::shutdown();
    return num_failed;
}
