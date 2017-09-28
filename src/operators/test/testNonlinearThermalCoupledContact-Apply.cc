#include "utils/AMPManager.h"
#include "utils/UnitTest.h"
#include "utils/Utilities.h"
#include <string>

#include "materials/Material.h"
#include "utils/AMPManager.h"
#include "utils/Database.h"
#include "utils/InputDatabase.h"
#include "utils/InputManager.h"
#include "utils/PIO.h"
#include "utils/Utilities.h"
#include "utils/shared_ptr.h"


#include "vectors/SimpleVector.h"
#include "vectors/Variable.h"
#include "vectors/Vector.h"

#include "ampmesh/MeshVariable.h"
#include "utils/Writer.h"

#include "operators/ColumnBoundaryOperator.h"
#include "operators/ElementOperationFactory.h"
#include "operators/ElementPhysicsModelFactory.h"
#include "operators/diffusion/DiffusionLinearElement.h"
#include "operators/diffusion/DiffusionLinearFEOperator.h"
#include "operators/diffusion/DiffusionNonlinearFEOperator.h"
#include "operators/diffusion/DiffusionTransportModel.h"
#include "operators/libmesh/MassLinearElement.h"
#include "operators/libmesh/MassLinearFEOperator.h"

#include "operators/CoupledOperator.h"
#include "operators/CoupledOperatorParameters.h"
#include "operators/DirichletMatrixCorrection.h"
#include "operators/DirichletVectorCorrection.h"
#include "operators/LinearBVPOperator.h"
#include "operators/NeumannVectorCorrection.h"
#include "operators/NeumannVectorCorrectionParameters.h"
#include "operators/NeutronicsRhs.h"
#include "operators/NonlinearBVPOperator.h"
#include "operators/OperatorBuilder.h"
#include "operators/RobinMatrixCorrection.h"
#include "operators/RobinVectorCorrection.h"
#include "operators/libmesh/VolumeIntegralOperator.h"
#include "operators/map/MapOperatorParameters.h"
#include "operators/map/MapSurface.h"

extern "C" {
#include "petsc.h"
}


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


void thermalContactApplyTest( AMP::UnitTest *ut, const std::string& exeName )
{
    std::string input_file = "input_" + exeName;
    std::string log_file   = "output_" + exeName;

    //  AMP::AMPManager::startup();
    //

    AMP::shared_ptr<AMP::InputDatabase> input_db( new AMP::InputDatabase( "input_db" ) );
    AMP::InputManager::getManager()->parseInputFile( input_file, input_db );
    input_db->printClassData( AMP::plog );

    AMP::PIO::logAllNodes( log_file );

    //  AMP_INSIST(input_db->keyExists("Mesh"), "Key ''Mesh'' is missing!");
    //  std::string mesh_file = input_db->getString("Mesh");

    AMP::Mesh::MeshManagerParameters::shared_ptr mgrParams(
        new AMP::Mesh::MeshManagerParameters( input_db ) );
    AMP::Mesh::MeshManager::shared_ptr manager( new AMP::Mesh::MeshManager( mgrParams ) );
    AMP::Mesh::MeshManager::Adapter::shared_ptr meshAdapter1 = manager->getMesh( "pellet" );
    AMP::Mesh::MeshManager::Adapter::shared_ptr meshAdapter2 = manager->getMesh( "clad" );

    AMP::LinearAlgebra::Vector::shared_ptr nullVec;

    AMP::LinearAlgebra::Variable::shared_ptr TemperatureVar(
        new AMP::Mesh::NodalScalarVariable( "Temperature" ) );
    AMP::LinearAlgebra::Variable::shared_ptr inputVariable1(
        new AMP::Mesh::NodalScalarVariable( "Temperature", meshAdapter1 ) );
    AMP::LinearAlgebra::Variable::shared_ptr inputVariable2(
        new AMP::Mesh::NodalScalarVariable( "Temperature", meshAdapter2 ) );

    AMP::LinearAlgebra::Variable::shared_ptr outputVariable1(
        new AMP::Mesh::NodalScalarVariable( "Temperature", meshAdapter1 ) );
    AMP::LinearAlgebra::Variable::shared_ptr outputVariable2(
        new AMP::Mesh::NodalScalarVariable( "Temperature", meshAdapter2 ) );

    double intguess = input_db->getDoubleWithDefault( "InitialGuess", 400 );

    AMP::LinearAlgebra::Vector::shared_ptr TemperatureInKelvin =
        manager->createVector( TemperatureVar );
    AMP::LinearAlgebra::Vector::shared_ptr RightHandSideVec =
        manager->createVector( TemperatureVar );
    AMP::LinearAlgebra::Vector::shared_ptr ResidualVec = manager->createVector( TemperatureVar );
    AMP::LinearAlgebra::Vector::shared_ptr WorkVec     = manager->createVector( TemperatureVar );

    TemperatureInKelvin->setToScalar( intguess );
    manager->registerVectorAsData( TemperatureInKelvin, "Temperature" );

    //-----------------------------------------------
    //   CREATE THE NONLINEAR THERMAL OPERATOR 1 ----
    //-----------------------------------------------

    AMP_INSIST( input_db->keyExists( "NonlinearThermalOperator1" ), "key missing!" );

    AMP::shared_ptr<AMP::Operator::ElementPhysicsModel> thermalTransportModel1;
    AMP::shared_ptr<AMP::Database> nonlinearThermalDatabase1 =
        input_db->getDatabase( "NonlinearThermalOperator1" );
    AMP::shared_ptr<AMP::Operator::NonlinearBVPOperator> nonlinearThermalOperator1 =
        AMP::dynamic_pointer_cast<AMP::Operator::NonlinearBVPOperator>(
            AMP::Operator::OperatorBuilder::createOperator(
                meshAdapter1, nonlinearThermalDatabase1, thermalTransportModel1 ) );

    // initialize the input variable
    AMP::shared_ptr<AMP::Operator::DiffusionNonlinearFEOperator> thermalVolumeOperator1 =
        AMP::dynamic_pointer_cast<AMP::Operator::DiffusionNonlinearFEOperator>(
            nonlinearThermalOperator1->getVolumeOperator() );

    // initialize the output variable
    // AMP::LinearAlgebra::Variable::shared_ptr outputVariable1 =
    // thermalVolumeOperator1->getOutputVariable();

    AMP::LinearAlgebra::Vector::shared_ptr TemperatureInKelvinVec1 =
        TemperatureInKelvin->subsetVectorForVariable( inputVariable1 );
    AMP::LinearAlgebra::Vector::shared_ptr RightHandSideVec1 =
        RightHandSideVec->subsetVectorForVariable( outputVariable1 );
    AMP::LinearAlgebra::Vector::shared_ptr ResidualVec1 =
        ResidualVec->subsetVectorForVariable( outputVariable1 );

    //-------------------------------------
    //   CREATE THE LINEAR THERMAL OPERATOR 1 ----
    //-------------------------------------

    AMP::shared_ptr<AMP::Operator::ElementPhysicsModel> transportModel1;
    AMP::shared_ptr<AMP::InputDatabase> bvpDatabase1 =
        AMP::dynamic_pointer_cast<AMP::InputDatabase>(
            input_db->getDatabase( "LinearThermalOperator1" ) );
    AMP::shared_ptr<AMP::Operator::LinearBVPOperator> linearThermalOperator1 =
        AMP::dynamic_pointer_cast<AMP::Operator::LinearBVPOperator>(
            AMP::Operator::OperatorBuilder::createOperator(
                meshAdapter1, bvpDatabase1, transportModel1 ) );

    //-------------------------------------
    //  CREATE THE NEUTRONICS SOURCE  //
    //-------------------------------------
    AMP_INSIST( input_db->keyExists( "NeutronicsOperator" ),
                "Key ''NeutronicsOperator'' is missing!" );
    AMP::shared_ptr<AMP::Database> neutronicsOp_db = input_db->getDatabase( "NeutronicsOperator" );
    AMP::shared_ptr<AMP::Operator::NeutronicsRhsParameters> neutronicsParams(
        new AMP::Operator::NeutronicsRhsParameters( neutronicsOp_db ) );
    neutronicsParams->d_MeshAdapter = meshAdapter1;
    AMP::shared_ptr<AMP::Operator::NeutronicsRhs> neutronicsOperator(
        new AMP::Operator::NeutronicsRhs( neutronicsParams ) );

    AMP::LinearAlgebra::Variable::shared_ptr SpecificPowerVar =
        neutronicsOperator->getOutputVariable();
    AMP::LinearAlgebra::Vector::shared_ptr SpecificPowerVec =
        meshAdapter1->createVector( SpecificPowerVar );

    neutronicsOperator->apply( nullVec, nullVec, SpecificPowerVec, 1., 0. );

    //----------------------------------------------------------
    //  Integrate Nuclear Rhs over Desnity * GeomType::Volume //
    //----------------------------------------------------------

    AMP_INSIST( input_db->keyExists( "VolumeIntegralOperator" ), "key missing!" );

    AMP::shared_ptr<AMP::Operator::ElementPhysicsModel> stransportModel;
    AMP::shared_ptr<AMP::Database> sourceDatabase =
        input_db->getDatabase( "VolumeIntegralOperator" );
    AMP::shared_ptr<AMP::Operator::VolumeIntegralOperator> sourceOperator =
        AMP::dynamic_pointer_cast<AMP::Operator::VolumeIntegralOperator>(
            AMP::Operator::OperatorBuilder::createOperator(
                meshAdapter1, sourceDatabase, stransportModel ) );

    // Create the power (heat source) vector.
    AMP::LinearAlgebra::Variable::shared_ptr PowerInWattsVar = sourceOperator->getOutputVariable();
    AMP::LinearAlgebra::Vector::shared_ptr PowerInWattsVec =
        meshAdapter1->createVector( PowerInWattsVar );
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
    AMP::shared_ptr<AMP::Database> nonlinearThermalDatabase2 =
        input_db->getDatabase( "NonlinearThermalOperator2" );
    AMP::shared_ptr<AMP::Operator::NonlinearBVPOperator> nonlinearThermalOperator2 =
        AMP::dynamic_pointer_cast<AMP::Operator::NonlinearBVPOperator>(
            AMP::Operator::OperatorBuilder::createOperator(
                meshAdapter2, nonlinearThermalDatabase2, thermalTransportModel2 ) );

    //----------------------------------------------------------------------------------------------------------------------------------------------//
    AMP::shared_ptr<AMP::Operator::DiffusionNonlinearFEOperator> thermalVolumeOperator2 =
        AMP::dynamic_pointer_cast<AMP::Operator::DiffusionNonlinearFEOperator>(
            nonlinearThermalOperator2->getVolumeOperator() );

    // initialize the output variable
    // AMP::LinearAlgebra::Variable::shared_ptr outputVariable2 =
    // thermalVolumeOperator2->getOutputVariable();

    AMP::LinearAlgebra::Vector::shared_ptr TemperatureInKelvinVec2 =
        TemperatureInKelvin->subsetVectorForVariable( inputVariable2 );
    AMP::LinearAlgebra::Vector::shared_ptr RightHandSideVec2 =
        RightHandSideVec->subsetVectorForVariable( outputVariable2 );
    AMP::LinearAlgebra::Vector::shared_ptr ResidualVec2 =
        ResidualVec->subsetVectorForVariable( outputVariable2 );

    //--------------------------------------------
    //   CREATE THE LINEAR THERMAL OPERATOR 2 ----
    //--------------------------------------------

    AMP::shared_ptr<AMP::Operator::ElementPhysicsModel> transportModel2;
    AMP::shared_ptr<AMP::InputDatabase> bvpDatabase2 =
        AMP::dynamic_pointer_cast<AMP::InputDatabase>(
            input_db->getDatabase( "LinearThermalOperator2" ) );
    AMP::shared_ptr<AMP::Operator::LinearBVPOperator> linearThermalOperator2 =
        AMP::dynamic_pointer_cast<AMP::Operator::LinearBVPOperator>(
            AMP::Operator::OperatorBuilder::createOperator(
                meshAdapter2, bvpDatabase2, transportModel2 ) );

    //-------------------------------------
    AMP::shared_ptr<AMP::InputDatabase> mapcladtopellet_db =
        AMP::dynamic_pointer_cast<AMP::InputDatabase>( input_db->getDatabase( "MapCladtoPellet" ) );
    AMP::shared_ptr<AMP::Operator::MapSurface> mapcladtopellet =
        AMP::dynamic_pointer_cast<AMP::Operator::MapSurface>(
            AMP::Operator::OperatorBuilder::createOperator(
                meshAdapter2, meshAdapter1, mapcladtopellet_db ) );

    AMP::shared_ptr<AMP::InputDatabase> mappellettoclad_db =
        AMP::dynamic_pointer_cast<AMP::InputDatabase>( input_db->getDatabase( "MapPellettoClad" ) );
    AMP::shared_ptr<AMP::Operator::MapSurface> mappellettoclad =
        AMP::dynamic_pointer_cast<AMP::Operator::MapSurface>(
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
    AMP::shared_ptr<AMP::Operator::NeumannVectorCorrectionParameters> correctionParameters1(
        new AMP::Operator::NeumannVectorCorrectionParameters( robinboundaryDatabase1 ) );


    //------------------------------------------

    AMP::Operator::Operator::shared_ptr boundaryOp2;
    boundaryOp2 = nonlinearThermalOperator2->getBoundaryOperator();

    AMP::Operator::Operator::shared_ptr robinBoundaryOp2;
    robinBoundaryOp2 =
        ( AMP::dynamic_pointer_cast<AMP::Operator::ColumnBoundaryOperator>( boundaryOp2 ) )
            ->getBoundaryOperator( 0 );

    AMP::shared_ptr<AMP::InputDatabase> boundaryDatabase2 =
        AMP::dynamic_pointer_cast<AMP::InputDatabase>(
            nonlinearThermalDatabase2->getDatabase( "BoundaryOperator" ) );
    AMP::shared_ptr<AMP::InputDatabase> robinboundaryDatabase2 =
        AMP::dynamic_pointer_cast<AMP::InputDatabase>(
            boundaryDatabase2->getDatabase( "RobinVectorCorrection" ) );

    robinboundaryDatabase2->putBool( "constant_flux", false );
    AMP::shared_ptr<AMP::Operator::NeumannVectorCorrectionParameters> correctionParameters2(
        new AMP::Operator::NeumannVectorCorrectionParameters( robinboundaryDatabase2 ) );


    //-------------------------------------
    // Applying boundary conditions to the nonlinear BVP Operator

    nonlinearThermalOperator1->modifyRHSvector( RightHandSideVec1 );
    nonlinearThermalOperator1->modifyInitialSolutionVector( TemperatureInKelvinVec1 );

    nonlinearThermalOperator2->modifyRHSvector( RightHandSideVec2 );
    nonlinearThermalOperator2->modifyInitialSolutionVector( TemperatureInKelvinVec2 );

    //-------------------------------------
    // Coupling Map to the Nonlinear Operators
    AMP::shared_ptr<AMP::InputDatabase> tmp_db( new AMP::InputDatabase( "Dummy" ) );
    AMP::shared_ptr<AMP::Operator::CoupledOperatorParameters> coupledNonlinearPelletParams(
        new AMP::Operator::CoupledOperatorParameters( tmp_db ) );

    coupledNonlinearPelletParams->d_MapOperator = mapcladtopellet;
    coupledNonlinearPelletParams->d_BVPOperator = nonlinearThermalOperator1;
    AMP::shared_ptr<AMP::Operator::CoupledOperator> coupledNonlinearPellet(
        new AMP::Operator::CoupledOperator( coupledNonlinearPelletParams ) );
    //-------------------------------------
    AMP::shared_ptr<AMP::Operator::CoupledOperatorParameters> coupledNonlinearCladParams(
        new AMP::Operator::CoupledOperatorParameters( tmp_db ) );
    coupledNonlinearCladParams->d_MapOperator = mappellettoclad;
    coupledNonlinearCladParams->d_BVPOperator = nonlinearThermalOperator2;
    AMP::shared_ptr<AMP::Operator::CoupledOperator> coupledNonlinearClad(
        new AMP::Operator::CoupledOperator( coupledNonlinearCladParams ) );

    //-------------------------------------
    // Column of Coupled Operators
    AMP::shared_ptr<AMP::Operator::OperatorParameters> nonlinearParams(
        new AMP::Operator::OperatorParameters( tmp_db ) );
    AMP::shared_ptr<AMP::Operator::ColumnOperator> nonlinearCoupledOperator(
        new AMP::Operator::ColumnOperator( nonlinearParams ) );
    nonlinearCoupledOperator->append( coupledNonlinearPellet );
    nonlinearCoupledOperator->append( coupledNonlinearClad );

    //-------------------------------------
    // Coupling Map to the Linear Operators
    AMP::shared_ptr<AMP::Operator::CoupledOperatorParameters> coupledLinearPelletParams(
        new AMP::Operator::CoupledOperatorParameters( tmp_db ) );
    coupledLinearPelletParams->d_MapOperator = mapcladtopellet;
    coupledLinearPelletParams->d_BVPOperator = linearThermalOperator1;
    AMP::shared_ptr<AMP::Operator::CoupledOperator> coupledLinearPellet(
        new AMP::Operator::CoupledOperator( coupledLinearPelletParams ) );
    //-------------------------------------
    AMP::shared_ptr<AMP::Operator::CoupledOperatorParameters> coupledLinearCladParams(
        new AMP::Operator::CoupledOperatorParameters( tmp_db ) );
    coupledLinearCladParams->d_MapOperator = mappellettoclad;
    coupledLinearCladParams->d_BVPOperator = linearThermalOperator2;
    AMP::shared_ptr<AMP::Operator::CoupledOperator> coupledLinearClad(
        new AMP::Operator::CoupledOperator( coupledLinearCladParams ) );

    //-------------------------------------
    // Column of Coupled Operators
    AMP::shared_ptr<AMP::Operator::OperatorParameters> linearParams(
        new AMP::Operator::OperatorParameters( tmp_db ) );
    AMP::shared_ptr<AMP::Operator::ColumnOperator> linearCoupledOperator(
        new AMP::Operator::ColumnOperator( linearParams ) );
    linearCoupledOperator->append( coupledLinearPellet );
    linearCoupledOperator->append( coupledLinearClad );

    //-------------------------------------

    bool testPassed = false;

    nonlinearCoupledOperator->apply(
        RightHandSideVec, TemperatureInKelvin, ResidualVec, 1.0, -1.0 );

    correctionParameters1->d_variableFlux = mapcladtopellet->getBoundaryVector();
    correctionParameters2->d_variableFlux = mappellettoclad->getBoundaryVector();

    robinBoundaryOp1->reset( correctionParameters1 );
    robinBoundaryOp2->reset( correctionParameters2 );

    testPassed = true;
    ut.passes( exeName + " : create" );

    // test apply
    std::string msgPrefix = exeName + " : apply";
    AMP::shared_ptr<AMP::Operator::Operator> testOperator =
        AMP::dynamic_pointer_cast<AMP::Operator::Operator>( nonlinearCoupledOperator );
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

    //-------------------------------------

    if ( testPassed ) {
        ut.passes( "Coupled Composite Nonlinear Operator Apply tests ." );
    } else {
        ITFAILS;
    }

    //} else {
    //  ut.expected_failure("parallel map3D-1D and map1D-3D fail in parallel, see bug #1219.");
    //}
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
