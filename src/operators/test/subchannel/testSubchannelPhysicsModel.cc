#include "utils/AMPManager.h"
#include "utils/UnitTest.h"
#include "utils/Utilities.h"
#include <iostream>
#include <string>

#include "utils/shared_ptr.h"

#include "utils/Database.h"
#include "utils/InputDatabase.h"
#include "utils/InputManager.h"

#include "utils/AMPManager.h"
#include "utils/AMP_MPI.h"
#include "utils/PIO.h"

#include "operators/ElementPhysicsModelFactory.h"
#include "operators/ElementPhysicsModelParameters.h"
#include "operators/subchannel/SubchannelPhysicsModel.h"

#include "materials/Material.h"

// function testing if two values are equal to within 1% and reporting if otherwise
bool areEqual( double result, double known, std::string quantity )
{
    bool passed;
    if ( !AMP::Utilities::approx_equal( result, known, 0.01 ) ) {
        AMP::pout << quantity << " result did not match known value.";
        AMP::pout << " Result = " << result << ", and Known = " << known << ".\n";
        passed = false;
    } else {
        passed = true;
    }
    return passed;
}

// function with tests to be performed on each input file
void nonlinearTest( AMP::UnitTest *ut, std::string exeName )
{

    //
    // test creation
    // =======================================================
    // initialize
    std::string input_file = "input_" + exeName;

    AMP::AMP_MPI globalComm( AMP_COMM_WORLD );

    std::cout << "\nInput file is " << input_file << std::endl;

    // parse input database file
    AMP::shared_ptr<AMP::InputDatabase> input_db( new AMP::InputDatabase( "input_db" ) );
    AMP::InputManager::getManager()->parseInputFile( input_file, input_db );
    input_db->printClassData( AMP::plog );

    AMP::shared_ptr<AMP::Database> subchannel_db =
        input_db->getDatabase( "SubchannelPhysicsModel" );
    AMP::shared_ptr<AMP::Operator::ElementPhysicsModelParameters> params(
        new AMP::Operator::ElementPhysicsModelParameters( subchannel_db ) );
    AMP::shared_ptr<AMP::Operator::SubchannelPhysicsModel> subchannelPhysicsModel(
        new AMP::Operator::SubchannelPhysicsModel( params ) );

    if ( subchannelPhysicsModel.get() ) {
        ut->passes( exeName + ": creation" );
    } else {
        ut->failure( exeName + ": creation" );
    }

    //
    // test function evaluations
    // =======================================================
    // create input argument maps
    const unsigned int n = 3;
    AMP::shared_ptr<std::vector<double>> enthalpyInput( new std::vector<double>( n ) );
    ( *enthalpyInput )[0] = 500.0e3;
    ( *enthalpyInput )[1] = 1.0e6;
    ( *enthalpyInput )[2] = 100.0e3;
    AMP::shared_ptr<std::vector<double>> pressureInput( new std::vector<double>( n ) );
    ( *pressureInput )[0] = 1.0e6;
    ( *pressureInput )[1] = 15.0e6;
    ( *pressureInput )[2] = 30.0e3;
    AMP::shared_ptr<std::vector<double>> temperatureInput( new std::vector<double>( n ) );
    ( *temperatureInput )[0] = 400.0;
    ( *temperatureInput )[1] = 600.0;
    ( *temperatureInput )[2] = 300.0;
    AMP::shared_ptr<std::vector<double>> temperatureInput2( new std::vector<double>( n ) );
    ( *temperatureInput2 )[0] = 392.140;
    ( *temperatureInput2 )[1] = 504.658;
    ( *temperatureInput2 )[2] = 297.004;
    AMP::shared_ptr<std::vector<double>> densityInput( new std::vector<double>( n ) );
    ( *densityInput )[0] = 937.871;
    ( *densityInput )[1] = 659.388;
    ( *densityInput )[2] = 996.526;
    std::map<std::string, AMP::shared_ptr<std::vector<double>>> temperatureArgs;
    temperatureArgs.insert( std::make_pair( "enthalpy", enthalpyInput ) );
    temperatureArgs.insert( std::make_pair( "pressure", pressureInput ) );
    std::map<std::string, AMP::shared_ptr<std::vector<double>>> liquidEnthalpyArgs;
    liquidEnthalpyArgs.insert( std::make_pair( "pressure", pressureInput ) );
    std::map<std::string, AMP::shared_ptr<std::vector<double>>> volumeArgs;
    volumeArgs.insert( std::make_pair( "enthalpy", enthalpyInput ) );
    volumeArgs.insert( std::make_pair( "pressure", pressureInput ) );
    std::map<std::string, AMP::shared_ptr<std::vector<double>>> conductivityArgs;
    conductivityArgs.insert( std::make_pair( "temperature", temperatureInput ) );
    conductivityArgs.insert( std::make_pair( "density", densityInput ) );
    std::map<std::string, AMP::shared_ptr<std::vector<double>>> viscosityArgs;
    viscosityArgs.insert( std::make_pair( "temperature", temperatureInput ) );
    viscosityArgs.insert( std::make_pair( "density", densityInput ) );
    std::map<std::string, AMP::shared_ptr<std::vector<double>>> enthalpyArgs;
    enthalpyArgs.insert( std::make_pair( "temperature", temperatureInput2 ) );
    enthalpyArgs.insert( std::make_pair( "pressure", pressureInput ) );

    // evaluate property functions
    std::vector<double> temperatureResult( n );
    subchannelPhysicsModel->getProperty( "Temperature", temperatureResult, temperatureArgs );
    std::vector<double> liquidEnthalpyResult( n );
    subchannelPhysicsModel->getProperty(
        "SaturatedLiquidEnthalpy", liquidEnthalpyResult, liquidEnthalpyArgs );
    std::vector<double> volumeResult( n );
    subchannelPhysicsModel->getProperty( "SpecificVolume", volumeResult, volumeArgs );
    std::vector<double> conductivityResult( n );
    subchannelPhysicsModel->getProperty(
        "ThermalConductivity", conductivityResult, conductivityArgs );
    std::vector<double> viscosityResult( n );
    subchannelPhysicsModel->getProperty( "DynamicViscosity", viscosityResult, viscosityArgs );
    std::vector<double> enthalpyResult( n );
    subchannelPhysicsModel->getProperty( "Enthalpy", enthalpyResult, enthalpyArgs );

    // compare property function results to known values
    bool passed                = true;
    double temperatureKnown[n] = { 392.140, 504.658, 297.004 };
    for ( unsigned int i = 0; i < n; ++i ) {
        passed = passed and areEqual( temperatureResult[i], temperatureKnown[i], "temperature" );
    }
    double liquidEnthalpyKnown[n] = { 762.683e3, 1610.15e3, 289.229e3 };
    for ( unsigned int i = 0; i < n; ++i ) {
        passed = passed and
                 areEqual( liquidEnthalpyResult[i], liquidEnthalpyKnown[i], "liquidEnthalpy" );
    }
    double volumeKnown[n] = { 0.00105925, 0.00119519, 0.00100259 };
    for ( unsigned int i = 0; i < n; ++i ) {
        passed = passed and areEqual( volumeResult[i], volumeKnown[i], "volume" );
    }
    double conductivityKnown[n] = { 0.684097, 0.503998, 0.610291 };
    for ( unsigned int i = 0; i < n; ++i ) {
        passed = passed and areEqual( conductivityResult[i], conductivityKnown[i], "conductivity" );
    }
    double viscosityKnown[n] = { 0.000218794, 7.72970e-5, 0.000853838 };
    for ( unsigned int i = 0; i < n; ++i ) {
        passed = passed and areEqual( viscosityResult[i], viscosityKnown[i], "viscosity" );
    }
    double enthalpyKnown[n] = { 500.0e3, 1.0e6, 100.0e3 };
    for ( unsigned int i = 0; i < n; ++i ) {
        passed = passed and areEqual( enthalpyResult[i], enthalpyKnown[i], "enthalpy" );
    }

    // determine if test passed or failed
    if ( passed )
        ut->passes( exeName + ": function evaluations" );
    else
        ut->failure( exeName + ": function evaluations" );

    //
    // test default argument function evaluations
    // =======================================================
    // create input argument maps
    enthalpyInput->resize( 1 );
    temperatureInput->resize( 1 );
    std::map<std::string, AMP::shared_ptr<std::vector<double>>> temperatureArgs1;
    temperatureArgs1.insert( std::make_pair( "enthalpy", enthalpyInput ) );
    std::map<std::string, AMP::shared_ptr<std::vector<double>>> temperatureArgs0;
    std::map<std::string, AMP::shared_ptr<std::vector<double>>> liquidEnthalpyArgs0;
    std::map<std::string, AMP::shared_ptr<std::vector<double>>> volumeArgs1;
    volumeArgs1.insert( std::make_pair( "enthalpy", enthalpyInput ) );
    std::map<std::string, AMP::shared_ptr<std::vector<double>>> volumeArgs0;
    std::map<std::string, AMP::shared_ptr<std::vector<double>>> conductivityArgs1;
    conductivityArgs1.insert( std::make_pair( "temperature", temperatureInput ) );
    std::map<std::string, AMP::shared_ptr<std::vector<double>>> conductivityArgs0;
    std::map<std::string, AMP::shared_ptr<std::vector<double>>> viscosityArgs1;
    viscosityArgs1.insert( std::make_pair( "temperature", temperatureInput ) );
    std::map<std::string, AMP::shared_ptr<std::vector<double>>> viscosityArgs0;
    std::map<std::string, AMP::shared_ptr<std::vector<double>>> enthalpyArgs1;
    enthalpyArgs1.insert( std::make_pair( "temperature", temperatureInput ) );
    std::map<std::string, AMP::shared_ptr<std::vector<double>>> enthalpyArgs0;

    // evaluate property functions
    std::vector<double> temperatureResult1Arg( 1 );
    subchannelPhysicsModel->getProperty( "Temperature", temperatureResult1Arg, temperatureArgs1 );
    std::vector<double> temperatureResult0Arg( 1 );
    subchannelPhysicsModel->getProperty( "Temperature", temperatureResult0Arg, temperatureArgs0 );
    std::vector<double> liquidEnthalpyResult0Arg( 1 );
    subchannelPhysicsModel->getProperty(
        "SaturatedLiquidEnthalpy", liquidEnthalpyResult0Arg, liquidEnthalpyArgs0 );
    std::vector<double> volumeResult1Arg( 1 );
    subchannelPhysicsModel->getProperty( "SpecificVolume", volumeResult1Arg, volumeArgs1 );
    std::vector<double> volumeResult0Arg( 1 );
    subchannelPhysicsModel->getProperty( "SpecificVolume", volumeResult0Arg, volumeArgs0 );
    std::vector<double> conductivityResult1Arg( 1 );
    subchannelPhysicsModel->getProperty(
        "ThermalConductivity", conductivityResult1Arg, conductivityArgs1 );
    std::vector<double> conductivityResult0Arg( 1 );
    subchannelPhysicsModel->getProperty(
        "ThermalConductivity", conductivityResult0Arg, conductivityArgs0 );
    std::vector<double> viscosityResult1Arg( 1 );
    subchannelPhysicsModel->getProperty( "DynamicViscosity", viscosityResult1Arg, viscosityArgs1 );
    std::vector<double> viscosityResult0Arg( 1 );
    subchannelPhysicsModel->getProperty( "DynamicViscosity", viscosityResult0Arg, viscosityArgs0 );
    std::vector<double> enthalpyResult1Arg( 1 );
    subchannelPhysicsModel->getProperty( "Enthalpy", enthalpyResult1Arg, enthalpyArgs1 );
    std::vector<double> enthalpyResult0Arg( 1 );
    subchannelPhysicsModel->getProperty( "Enthalpy", enthalpyResult0Arg, enthalpyArgs0 );

    // compare property function results to known values
    passed                      = true;
    double temperatureKnown1Arg = 392.224;
    passed =
        passed and areEqual( temperatureResult1Arg[0], temperatureKnown1Arg, "temperature, 1 arg" );
    double temperatureKnown0Arg = 320.835;
    passed =
        passed and areEqual( temperatureResult0Arg[0], temperatureKnown0Arg, "temperature, 0 arg" );
    double liquidEnthalpyKnown0Arg = 640.185e3;
    passed =
        passed and
        areEqual( liquidEnthalpyResult0Arg[0], liquidEnthalpyKnown0Arg, "liquidEnthalpy, 0 arg" );
    double volumeKnown1Arg = 0.00105962;
    passed = passed and areEqual( volumeResult1Arg[0], volumeKnown1Arg, "volume, 1 arg" );
    double volumeKnown0Arg = 0.00101083;
    passed = passed and areEqual( volumeResult0Arg[0], volumeKnown0Arg, "volume, 0 arg" );
    double conductivityKnown1Arg = 0.731;
    passed                       = passed and
             areEqual( conductivityResult1Arg[0], conductivityKnown1Arg, "conductivity, 1 arg" );
    double conductivityKnown0Arg = 0.668247;
    passed                       = passed and
             areEqual( conductivityResult0Arg[0], conductivityKnown0Arg, "conductivity, 0 arg" );
    double viscosityKnown1Arg = 0.000239;
    passed = passed and areEqual( viscosityResult1Arg[0], viscosityKnown1Arg, "viscosity, 1 arg" );
    double viscosityKnown0Arg = 0.000368895;
    passed = passed and areEqual( viscosityResult0Arg[0], viscosityKnown0Arg, "viscosity, 0 arg" );
    double enthalpyKnown1Arg = 533.121e3;
    passed = passed and areEqual( enthalpyResult1Arg[0], enthalpyKnown1Arg, "enthalpy, 1 arg" );
    double enthalpyKnown0Arg = 322.100e3;
    passed = passed and areEqual( enthalpyResult0Arg[0], enthalpyKnown0Arg, "enthalpy, 0 arg" );

    // determine if test passed or failed
    if ( passed )
        ut->passes( exeName + ": defaults evaluations" );
    else
        ut->failure( exeName + ": defaults evaluations" );
}

int main( int argc, char *argv[] )
{
    AMP::AMPManagerProperties startup_properties;
    startup_properties.use_MPI_Abort = false;
    AMP::AMPManager::startup( argc, argv, startup_properties );

    AMP::UnitTest ut;

    const int NUMFILES          = 1;
    std::string files[NUMFILES] = { "testSubchannelPhysicsModel" };

    for ( auto &file : files ) {
        try {
            nonlinearTest( &ut, file );
        } catch ( std::exception &err ) {
            std::cout << "ERROR: While testing " << argv[0] << err.what() << std::endl;
            ut.failure( "ERROR: While testing: " + file );
        } catch ( ... ) {
            std::cout << "ERROR: While testing " << argv[0] << "An unknown exception was thrown."
                      << std::endl;
            ut.failure( "ERROR: While testing: " + file );
        }
    }

    ut.report();

    int num_failed = ut.NumFailGlobal();
    AMP::AMPManager::shutdown();
    return num_failed;
}
