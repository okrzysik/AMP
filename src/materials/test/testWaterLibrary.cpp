/*
 * \file materials/test/testWaterLibrary.cc
 * \brief test of WaterLibrary.cc class
 */


#include <iostream>
#include <string>
#include <valarray>
#include <vector>

#include "AMP/materials/Material.h"
#include "AMP/utils/AMPManager.h"
#include "AMP/utils/UnitTest.h"
#include "AMP/utils/Utilities.h"


using materialMap = std::map<std::string, std::shared_ptr<std::vector<double>>>;
static inline auto make_shared_vector( size_t N, double v = 0 )
{
    return std::make_shared<std::vector<double>>( N, v );
}


void checkConsistency( double h, double p, double T, bool &allCorrect, bool &allConsistent )
{
    auto mat                 = AMP::Materials::getMaterial( "WaterLibrary" ); // get water library
    auto temperatureProperty = mat->property( "Temperature" ); // temperature property
    auto enthalpyProperty    = mat->property( "Enthalpy" );    // enthalpy property

    double tempOutput = temperatureProperty->eval( {}, "enthalpy", h, "pressure", p );
    // check that answer is correct
    if ( !AMP::Utilities::approx_equal( tempOutput, T, 0.01 ) ) {
        AMP::pout << "Incorrect value: Calculated T: " << tempOutput << ", Actual T: " << T
                  << std::endl;
        allCorrect = false;
    }
    // check that enthalpy function resturns original enthalpy
    double hOutput = enthalpyProperty->eval( {}, "temperature", tempOutput, "pressure", p );
    if ( !AMP::Utilities::approx_equal( hOutput, h, 0.01 ) )
        allConsistent = false;
}

int main( int argc, char **argv )
{
    AMP::AMPManager::startup( argc, argv );
    AMP::UnitTest ut;

    bool good = true;

    // test constructors for temperature
    auto mat                 = AMP::Materials::getMaterial( "WaterLibrary" ); // get water library
    auto temperatureProperty = mat->property( "Temperature" );                // temperature
    auto liquidEnthalpyProperty =
        mat->property( "SaturatedLiquidEnthalpy" ); // saturated liquid enthalpy
    auto vaporEnthalpyProperty =
        mat->property( "SaturatedVaporEnthalpy" );                      // saturated vapor enthalpy
    auto volumeProperty       = mat->property( "SpecificVolume" );      // specific volume
    auto conductivityProperty = mat->property( "ThermalConductivity" ); // thermal conductivity
    auto viscosityProperty    = mat->property( "DynamicViscosity" );    // dynamic viscosity
    auto enthalpyProperty     = mat->property( "Enthalpy" );            // enthalpy

    // test property accessors for temperature
    auto tcname = temperatureProperty->get_name();
    auto tcsorc = temperatureProperty->get_source();
    AMP::pout << "\n";
    good = good && tcname == "WaterLibrary_Temperature";
    AMP::pout << "Temperature name is " << tcname << "\n";
    AMP::pout << "Temperature source is " << tcsorc << "\n";
    auto args = temperatureProperty->get_arguments();
    good      = good && args[0] == "enthalpy";
    good      = good && args[1] == "pressure";
    AMP::pout << "Temperature property arguments are " << args[0] << " and " << args[1] << "\n\n";
    unsigned int nargs = temperatureProperty->get_number_arguments();
    good               = good && nargs == 2;

    // test material accessors, all arguments present
    const size_t n =
        3; // size of input and output arrays for comparison with known thermodynamic values
    auto enthalpyInput    = make_shared_vector( n ); // enthalpy input
    auto pressureInput    = make_shared_vector( n ); // pressure input
    auto temperatureInput = make_shared_vector( n ); // temperature input
    auto temperatureInputEnthalpy =
        make_shared_vector( n );                 // temperature input for enthalpy function
    auto densityInput = make_shared_vector( n ); // density input
    auto temperatureIdenticalInput =
        make_shared_vector( n ); // temperature input array with identical values
    auto enthalpyIdenticalInput =
        make_shared_vector( n ); // enthalpy input array with identical values
    auto pressureIdenticalInput =
        make_shared_vector( n );                   // pressure input array with identical values
    std::vector<double> temperatureOutput( n );    // temperature output
    std::vector<double> liquidEnthalpyOutput( n ); // saturated liquid enthalpy output
    std::vector<double> volumeOutput( n );         // specific volume output
    std::vector<double> conductivityOutput( n );   // thermal conductivity output
    std::vector<double> viscosityOutput( n );      // dynamic viscosity output
    std::vector<double> enthalpyOutput( n );       // enthalpy output
    std::vector<double> temperatureIdenticalOutput(
        n ); // temperature output array with identical values
    for ( size_t i = 0; i < n; i++ ) {
        ( *temperatureIdenticalInput )[i] = 400.0;   // temperature: 400 K
        ( *enthalpyIdenticalInput )[i]    = 500.0e3; // enthalpy: 500 kJ/kg
        ( *pressureIdenticalInput )[i]    = 1.0e6;   // pressure: 1 MPa
    }
    ( *enthalpyInput )[0] = 500.0e3;
    ( *enthalpyInput )[1] = 1.0e6;
    ( *enthalpyInput )[2] = 100.0e3;

    ( *pressureInput )[0] = 1.0e6;
    ( *pressureInput )[1] = 15.0e6;
    ( *pressureInput )[2] = 30.0e3;

    ( *temperatureInput )[0] = 400.0;
    ( *temperatureInput )[1] = 600.0;
    ( *temperatureInput )[2] = 300.0;

    ( *temperatureInputEnthalpy )[0] = 392.140;
    ( *temperatureInputEnthalpy )[1] = 504.658;
    ( *temperatureInputEnthalpy )[2] = 297.004;

    ( *densityInput )[0] = 937.871;
    ( *densityInput )[1] = 659.388;
    ( *densityInput )[2] = 996.526;

    // Block for temporary variables
    {
        // argument maps for each property function
        // temperature
        materialMap temperatureArgMap;
        temperatureArgMap.insert( std::make_pair( "enthalpy", enthalpyInput ) );
        temperatureArgMap.insert( std::make_pair( "pressure", pressureInput ) );
        // saturated liquid enthalpy
        materialMap liquidEnthalpyArgMap;
        liquidEnthalpyArgMap.insert( std::make_pair( "pressure", pressureInput ) );
        // specific volume
        materialMap volumeArgMap;
        volumeArgMap.insert( std::make_pair( "enthalpy", enthalpyInput ) );
        volumeArgMap.insert( std::make_pair( "pressure", pressureInput ) );
        // thermal conductivity
        materialMap conductivityArgMap;
        conductivityArgMap.insert( std::make_pair( "temperature", temperatureInput ) );
        conductivityArgMap.insert( std::make_pair( "density", densityInput ) );
        // dynamic viscosity
        materialMap viscosityArgMap;
        viscosityArgMap.insert( std::make_pair( "temperature", temperatureInput ) );
        viscosityArgMap.insert( std::make_pair( "density", densityInput ) );
        // enthalpy
        materialMap enthalpyArgMap;
        enthalpyArgMap.insert( std::make_pair( "temperature", temperatureInputEnthalpy ) );
        enthalpyArgMap.insert( std::make_pair( "pressure", pressureInput ) );
        // temperature identical values case
        materialMap temperatureIdenticalArgMap;
        temperatureIdenticalArgMap.insert( std::make_pair( "enthalpy", enthalpyIdenticalInput ) );
        temperatureIdenticalArgMap.insert( std::make_pair( "pressure", pressureIdenticalInput ) );

        // evaluate properties
        // temperature
        temperatureProperty->evalv( temperatureOutput, {}, temperatureArgMap );
        // saturated liquid enthalpy
        liquidEnthalpyProperty->evalv( liquidEnthalpyOutput, {}, liquidEnthalpyArgMap );
        // specific volume
        volumeProperty->evalv( volumeOutput, {}, volumeArgMap );
        // thermal conductivity
        conductivityProperty->evalv( conductivityOutput, {}, conductivityArgMap );
        // dynamic viscosity
        viscosityProperty->evalv( viscosityOutput, {}, viscosityArgMap );
        // enthalpy
        enthalpyProperty->evalv( enthalpyOutput, {}, enthalpyArgMap );
        // temperature identical values case
        std::vector<double> temperatureOutput_mat( temperatureIdenticalOutput );
        mat->property( "Temperature" )
            ->evalv( temperatureOutput_mat, {}, temperatureIdenticalArgMap );
        temperatureProperty->evalv( temperatureIdenticalOutput, {}, temperatureIdenticalArgMap );

        // known values for testing
        double temperatureKnown[n]    = { 392.140, 504.658, 297.004 }; // temperature
        double liquidEnthalpyKnown[n] = { 762.683e3,
                                          1610.15e3,
                                          289.229e3 }; // saturated liquid enthalpy
        double volumeKnown[n]         = { 0.00105925, 0.00119519, 0.00100259 }; // specific volume
        double conductivityKnown[n]   = { 0.684097, 0.503998, 0.610291 };    // thermal conductivity
        double viscosityKnown[n] = { 0.000218794, 7.72970e-5, 0.000853838 }; // dynamic viscosity
        double enthalpyKnown[n]  = { 500.0e3, 1.0e6, 100.0e3 };              // enthalpy

        // test property functions against known values
        for ( size_t i = 0; i < n; i++ ) {
            AMP::pout << "\nValue test " << i << ":\n=====================\n";
            // temperature
            if ( !AMP::Utilities::approx_equal(
                     temperatureOutput[i], temperatureKnown[i], 0.01 ) ) {
                ut.failure( "The answer is wrong." );
                AMP::pout << "The calculated temperature was " << temperatureOutput[i]
                          << " K and actual is ";
                AMP::pout << temperatureKnown[i] << " K\n";
            } else
                AMP::pout << "temperature value is approximately equal to known value.\n";
            // saturated liquid enthalpy
            if ( !AMP::Utilities::approx_equal(
                     liquidEnthalpyOutput[i], liquidEnthalpyKnown[i], 0.01 ) ) {
                ut.failure( "The answer is wrong." );
                AMP::pout << "The calculated saturated liquid enthalpy was "
                          << liquidEnthalpyOutput[i] << " J/kg and actual is ";
                AMP::pout << liquidEnthalpyKnown[i] << " J/kg\n";
            } else
                AMP::pout << "saturated liquid enthalpy value is approximately equal to known "
                             "value.\n";
            // specific volume
            if ( !AMP::Utilities::approx_equal( volumeOutput[i], volumeKnown[i], 0.01 ) ) {
                ut.failure( "The answer is wrong." );
                AMP::pout << "The calculated specific volume was " << volumeOutput[i]
                          << " m3/kg and actual is ";
                AMP::pout << volumeKnown[i] << " m3/kg\n";
            } else
                AMP::pout << "specific volume value is approximately equal to known value.\n";
            // thermal conductivity
            if ( !AMP::Utilities::approx_equal(
                     conductivityOutput[i], conductivityKnown[i], 0.01 ) ) {
                ut.failure( "The answer is wrong." );
                AMP::pout << "The calculated thermal conductivity was " << conductivityOutput[i]
                          << " W/m-K and actual is ";
                AMP::pout << conductivityKnown[i] << " W/m-K\n";
            } else
                AMP::pout << "thermal conductivity value is approximately equal to known value.\n";
            // dynamic viscosity
            if ( !AMP::Utilities::approx_equal( viscosityOutput[i], viscosityKnown[i], 0.01 ) ) {
                ut.failure( "The answer is wrong." );
                AMP::pout << "The calculated dynamic viscosity was " << viscosityOutput[i]
                          << " Pa-s and actual is ";
                AMP::pout << viscosityKnown[i] << " Pa-s\n";
            } else
                AMP::pout << "dynamic viscosity value is approximately equal to known value.\n";
            // enthalpy
            if ( !AMP::Utilities::approx_equal( enthalpyOutput[i], enthalpyKnown[i], 0.01 ) ) {
                ut.failure( "The answer is wrong." );
                AMP::pout << "The calculated enthalpy was " << enthalpyOutput[i]
                          << " J/kg and actual is ";
                AMP::pout << enthalpyKnown[i] << " J/kg\n";
            } else
                AMP::pout << "enthalpy value is approximately equal to known value.\n";
        }

        // identical values test: compare temperature values against each other
        for ( size_t i = 0; i < n; i++ ) {
            if ( !AMP::Utilities::approx_equal( temperatureIdenticalOutput[i],
                                                temperatureOutput_mat[i] ) ) {
                AMP::pout << "Identical values temperature test 1 failed: 1st value: "
                          << temperatureIdenticalOutput[i];
                AMP::pout << " and 2nd value: " << temperatureOutput_mat[i] << "\n";
                good = false;
            }
            if ( !AMP::Utilities::approx_equal( temperatureIdenticalOutput[0],
                                                temperatureIdenticalOutput[i] ) ) {
                AMP::pout << "Identical values temperature test 2 failed: 1st value: "
                          << temperatureIdenticalOutput[0];
                AMP::pout << " and 2nd value: " << temperatureIdenticalOutput[i] << "\n";
                good = false;
            }
        }
    }

    AMP::pout << "\nDefault arguments tests:\n============================\n";
    // set defaults
    // temperature
    std::vector<double> temperatureDefaults( 2 );
    temperatureDefaults[0] = 200e3; // enthalpy: 200 kJ/kg
    temperatureDefaults[1] = 0.5e6; // pressure: 0.5 MPa
    temperatureProperty->set_defaults( temperatureDefaults );
    // saturated liquid enthalpy
    std::vector<double> liquidEnthalpyDefaults( 1 );
    liquidEnthalpyDefaults[0] = 0.5e6; // pressure: 0.5 MPa
    liquidEnthalpyProperty->set_defaults( liquidEnthalpyDefaults );
    // specific volume
    std::vector<double> volumeDefaults( 2 );
    volumeDefaults[0] = 200e3; // enthalpy: 200 kJ/kg
    volumeDefaults[1] = 0.5e6; // pressure: 0.5 MPa
    volumeProperty->set_defaults( volumeDefaults );
    // thermal conductivity
    std::vector<double> conductivityDefaults( 2 );
    conductivityDefaults[0] = 350;     // temperature: 350 K
    conductivityDefaults[1] = 973.919; // density: 973.919 kg/m3
    conductivityProperty->set_defaults( conductivityDefaults );
    // dynamic viscosity
    std::vector<double> viscosityDefaults( 2 );
    viscosityDefaults[0] = 350;     // temperature: 350 K
    viscosityDefaults[1] = 973.919; // density: 973.919 kg/m3
    viscosityProperty->set_defaults( viscosityDefaults );
    // enthalpy
    std::vector<double> enthalpyDefaults( 2 );
    enthalpyDefaults[0] = 320.835; // temperature: 320.835
    enthalpyDefaults[1] = 0.5e6;   // pressure: 0.5 MPa
    enthalpyProperty->set_defaults( enthalpyDefaults );

    // default testing
    // temperature, one argument: enthalpy
    {
        double knownSolution = 392.224; // T [K]	@ {0.5 MPa, 500 kJ/kg}
        materialMap argMap;
        argMap.insert( std::make_pair( "enthalpy", enthalpyIdenticalInput ) );
        std::vector<double> temperatureOutput_def( temperatureOutput );
        temperatureProperty->evalv( temperatureOutput_def, {}, argMap );
        if ( !AMP::Utilities::approx_equal( temperatureOutput_def[0], knownSolution, 0.01 ) ) {
            AMP::pout << "Temperature w/ 1 arg incorrect: ";
            AMP::pout << temperatureOutput_def[0] << " vs " << knownSolution << "\n";
            good = false;
        }
    }
    // temperature, no argument
    {
        double knownSolution = 320.835; // T [K]	@ {0.5 MPa, 200 kJ/kg}
        materialMap argMap;
        std::vector<double> temperatureOutput_def( temperatureOutput );
        temperatureProperty->evalv( temperatureOutput_def, {}, argMap );
        if ( !AMP::Utilities::approx_equal( temperatureOutput_def[0], knownSolution, 0.01 ) ) {
            AMP::pout << "Temperature w/ 0 arg incorrect: ";
            AMP::pout << temperatureOutput_def[0] << " vs " << knownSolution << "\n";
            good = false;
        }
    }
    // saturated liquid enthalpy, no argument
    {
        double knownSolution = 640.185e3; // Hf [J/kg]	@ {0.5 MPa}
        materialMap argMap;
        std::vector<double> liquidEnthalpyOutput_def( liquidEnthalpyOutput );
        liquidEnthalpyProperty->evalv( liquidEnthalpyOutput_def, {}, argMap );
        if ( !AMP::Utilities::approx_equal( liquidEnthalpyOutput_def[0], knownSolution, 0.01 ) ) {
            AMP::pout << "Saturated liquid enthalpy w/ 0 arg incorrect: ";
            AMP::pout << liquidEnthalpyOutput_def[0] << " vs " << knownSolution << "\n";
            good = false;
        }
    }
    // specific volume, one argument: enthalpy
    {
        double knownSolution = 0.00105962; // v [m3/kg]	@ {0.5 MPa, 500 kJ/kg}
        materialMap argMap;
        argMap.insert( std::make_pair( "enthalpy", enthalpyIdenticalInput ) );
        std::vector<double> volumeOutput_def( volumeOutput );
        volumeProperty->evalv( volumeOutput_def, {}, argMap );
        if ( !AMP::Utilities::approx_equal( volumeOutput_def[0], knownSolution, 0.01 ) ) {
            AMP::pout << "Specific volume w/ 1 arg incorrect: ";
            AMP::pout << volumeOutput_def[0] << " vs " << knownSolution << "\n";
            good = false;
        }
    }
    // specific volume, no argument
    {
        double knownSolution = 0.00101083; // v [m3/kg]	@ {0.5 MPa, 200 kJ/kg}
        materialMap argMap;
        std::vector<double> volumeOutput_def( volumeOutput );
        volumeProperty->evalv( volumeOutput_def, {}, argMap );
        if ( !AMP::Utilities::approx_equal( volumeOutput_def[0], knownSolution, 0.01 ) ) {
            AMP::pout << "Specific volume w/ 0 arg incorrect: ";
            AMP::pout << volumeOutput_def[0] << " vs " << knownSolution << "\n";
            good = false;
        }
    }
    // thermal conductivity, one argument: enthalpy
    {
        double knownSolution = 0.731; // k [W/m-K]	@ {400 K, 973.919 kg/m3}
        materialMap argMap;
        argMap.insert( std::make_pair( "temperature", temperatureIdenticalInput ) );
        std::vector<double> conductivityOutput_def( conductivityOutput );
        conductivityProperty->evalv( conductivityOutput_def, {}, argMap );
        if ( !AMP::Utilities::approx_equal( conductivityOutput_def[0], knownSolution, 0.01 ) ) {
            AMP::pout << "Thermal conductivity w/ 1 arg incorrect: ";
            AMP::pout << conductivityOutput_def[0] << " vs " << knownSolution << "\n";
            good = false;
        }
    }
    // thermal conductivity, no argument
    {
        double knownSolution = 0.668247; // k [W/m-K]	@ {350 K, 973.919 kg/m3}
        materialMap argMap;
        std::vector<double> conductivityOutput_def( conductivityOutput );
        conductivityProperty->evalv( conductivityOutput_def, {}, argMap );
        if ( !AMP::Utilities::approx_equal( conductivityOutput_def[0], knownSolution, 0.01 ) ) {
            AMP::pout << "Thermal conductivity w/ 0 arg incorrect: ";
            AMP::pout << conductivityOutput_def[0] << " vs " << knownSolution << "\n";
            good = false;
        }
    }
    // dynamic viscosity, one argument: enthalpy
    {
        double knownSolution = 0.000239; // u [Pa-s]	@ {400 K, 973.919 kg/m3}
        materialMap argMap;
        argMap.insert( std::make_pair( "temperature", temperatureIdenticalInput ) );
        std::vector<double> viscosityOutput_def( viscosityOutput );
        viscosityProperty->evalv( viscosityOutput_def, {}, argMap );
        if ( !AMP::Utilities::approx_equal( viscosityOutput_def[0], knownSolution, 0.01 ) ) {
            AMP::pout << "Dynamic viscosity w/ 1 arg incorrect: ";
            AMP::pout << viscosityOutput_def[0] << " vs " << knownSolution << "\n";
            good = false;
        }
    }
    // dynamic viscosity, no argument
    {
        double knownSolution = 0.000368895; // u [Pa-s]	@ {350 K, 973.919 kg/m3}
        materialMap argMap;
        std::vector<double> viscosityOutput_def( viscosityOutput );
        viscosityProperty->evalv( viscosityOutput_def, {}, argMap );
        if ( !AMP::Utilities::approx_equal( viscosityOutput_def[0], knownSolution, 0.01 ) ) {
            AMP::pout << "Dynamic viscosity w/ 0 arg incorrect: ";
            AMP::pout << viscosityOutput_def[0] << " vs " << knownSolution << "\n";
            good = false;
        }
    }
    // enthalpy, one argument: temperature
    {
        double knownSolution = 533.121e3; // h [J/kg]	@ {400 K, 0.5 MPa}
        materialMap argMap;
        argMap.insert( std::make_pair( "temperature", temperatureIdenticalInput ) );
        std::vector<double> enthalpyOutput_def( enthalpyOutput );
        enthalpyProperty->evalv( enthalpyOutput_def, {}, argMap );
        if ( !AMP::Utilities::approx_equal( enthalpyOutput_def[0], knownSolution, 0.01 ) ) {
            AMP::pout << "Enthalpy w/ 1 arg incorrect: ";
            AMP::pout << enthalpyOutput_def[0] << " vs " << knownSolution << "\n";
            good = false;
        }
    }
    // enthalpy, no argument
    {
        double knownSolution = 200.0e3; // h [J/kg]	@ {320.835 K, 0.5 MPa}
        materialMap argMap;
        std::vector<double> enthalpyOutput_def( enthalpyOutput );
        enthalpyProperty->evalv( enthalpyOutput_def, {}, argMap );
        if ( !AMP::Utilities::approx_equal( enthalpyOutput_def[0], knownSolution, 0.01 ) ) {
            AMP::pout << "Enthalpy w/ 0 arg incorrect: ";
            AMP::pout << enthalpyOutput_def[0] << " vs " << knownSolution << "\n";
            good = false;
        }
    }

    if ( good )
        ut.passes( "basic tests of Material" );
    else
        ut.failure( "basic tests of Material" );

    // Test if thermodyamic properties are consistent
    bool allCorrect    = true;
    bool allConsistent = true;

    checkConsistency( 430.39e3, 15.0e6, 373.15, allCorrect, allConsistent );
    checkConsistency( 1334.4e3, 20.0e6, 573.15, allCorrect, allConsistent );
    checkConsistency( 176.37e3, 10.0e6, 313.15, allCorrect, allConsistent );
    checkConsistency( 507.19e3, 5.0e6, 393.15, allCorrect, allConsistent );
    checkConsistency( 684.01e3, 15.0e6, 433.15, allCorrect, allConsistent );

    if ( allCorrect )
        ut.passes( "Thermodynamic property value test" );
    else
        ut.failure( "Thermodynamic property value test" );
    if ( allConsistent )
        ut.passes( "Thermodynamic property consistency test" );
    else
        ut.failure( "Thermodynamic property consistency test" );

    // Extra tests of extended-range Temperature and Specific GeomType::Cell properties
    // Comparisons are to Matlab-computed values from Retran-3D source (see WaterLibrary.cc)
    // Matlab script used to generate values is located in data/waterlibrary.m
    AMP::pout << "\nExtended library tests:\n============================\n";
    bool extra_passed = true;
    auto pvec         = make_shared_vector( 1 );
    auto hvec         = make_shared_vector( 1 );
    double expected;

    // Test saturated liquid enthalpy
    materialMap pArgs;
    pArgs.insert( std::make_pair( "pressure", pvec ) );
    std::vector<double> h( 1 );
    double tol = 1e-6;

    // p=800 psi -- Eq. 5a
    ( *pvec )[0] = 5.51580582e6;
    expected     = 1.18586462e6;
    liquidEnthalpyProperty->evalv( h, {}, pArgs );
    if ( !AMP::Utilities::approx_equal( h[0], expected, tol ) ) {
        AMP::pout << "Saturated liquid enthalpy test 1 failed." << std::endl;
        AMP::pout << "Expected " << expected << ", computed " << h[0] << std::endl;
        extra_passed = false;
    }

    // p=2000 psi -- Eq. 5b
    ( *pvec )[0] = 1.37895146e7;
    expected     = 1.56325778e6;
    liquidEnthalpyProperty->evalv( h, {}, pArgs );
    if ( !AMP::Utilities::approx_equal( h[0], expected, tol ) ) {
        AMP::pout << "Saturated liquid enthalpy test 2 failed." << std::endl;
        AMP::pout << "Expected " << expected << ", computed " << h[0] << std::endl;
        extra_passed = false;
    }

    // p=3000 psi -- Eq. 5c
    ( *pvec )[0] = 2.06842718e7;
    expected     = 1.86516705e6;
    liquidEnthalpyProperty->evalv( h, {}, pArgs );
    if ( !AMP::Utilities::approx_equal( h[0], expected, tol ) ) {
        AMP::pout << "Saturated liquid enthalpy test 3 failed." << std::endl;
        AMP::pout << "Expected " << expected << ", computed " << h[0] << std::endl;
        extra_passed = false;
    }

    // Test saturated vapor enthalpy

    // p=800 psi -- Eq. 6a
    ( *pvec )[0] = 5.51580582e6;
    expected     = 2.78981197e6;
    vaporEnthalpyProperty->evalv( h, {}, pArgs );
    if ( !AMP::Utilities::approx_equal( h[0], expected, tol ) ) {
        AMP::pout << "Saturated vapor enthalpy test 1 failed." << std::endl;
        AMP::pout << "Expected " << expected << ", computed " << h[0] << std::endl;
        extra_passed = false;
    }

    // p = 2000 psi -- Eq. 6b
    ( *pvec )[0] = 1.37895146e7;
    expected     = 2.64769923e6;
    vaporEnthalpyProperty->evalv( h, {}, pArgs );
    if ( !AMP::Utilities::approx_equal( h[0], expected, tol ) ) {
        AMP::pout << "Saturated vapor enthalpy test 2 failed." << std::endl;
        AMP::pout << "Expected " << expected << ", computed " << h[0] << std::endl;
        extra_passed = false;
    }

    // p=3000 psi -- Eq. 5c
    ( *pvec )[0] = 2.06842718e7;
    expected     = 2.37307137e6;
    vaporEnthalpyProperty->evalv( h, {}, pArgs );
    if ( !AMP::Utilities::approx_equal( h[0], expected, tol ) ) {
        AMP::pout << "Saturated vapor enthalpy test 3 failed." << std::endl;
        AMP::pout << "Expected " << expected << ", computed " << h[0] << std::endl;
        extra_passed = false;
    }

    // Test TemperatureProp
    materialMap hpArgs;
    hpArgs.insert( std::make_pair( "enthalpy", hvec ) );
    hpArgs.insert( std::make_pair( "pressure", pvec ) );
    std::vector<double> T( 1 );

    // h = 400 Btu/lbm, p = 4000 psi -- Eq. 6b
    ( *hvec )[0] = 9.304e5;
    ( *pvec )[0] = 2.75790291e7;
    expected     = 488.207658;
    temperatureProperty->evalv( T, {}, hpArgs );
    if ( !AMP::Utilities::approx_equal( T[0], expected, tol ) ) {
        AMP::pout << "Extended temperature test 1 failed." << std::endl;
        AMP::pout << "Expected " << expected << ", computed " << T[0] << std::endl;
        extra_passed = false;
    }

    // h = 1400 Btu/lbm, p = 4000 psi -- Eq. 6d
    ( *hvec )[0] = 3.2564e6;
    ( *pvec )[0] = 2.75790291e7;
    expected     = 808.434684;
    temperatureProperty->evalv( T, {}, hpArgs );
    if ( !AMP::Utilities::approx_equal( T[0], expected, tol ) ) {
        AMP::pout << "Extended temperature test 2 failed." << std::endl;
        AMP::pout << "Expected " << expected << ", computed " << T[0] << std::endl;
        extra_passed = false;
    }

    // h = 400 Btu/lbm, p = 2000 psi -- Eq. 6a
    ( *hvec )[0] = 9.304e5;
    ( *pvec )[0] = 1.37895146e7;
    expected     = 489.651055;
    temperatureProperty->evalv( T, {}, hpArgs );
    if ( !AMP::Utilities::approx_equal( T[0], expected, tol ) ) {
        AMP::pout << "Extended temperature test 3 failed." << std::endl;
        AMP::pout << "Expected " << expected << ", computed " << T[0] << std::endl;
        extra_passed = false;
    }

    // h = 800 Btu/lbm, p = 2000 psi -- Eq. 6a at saturated liquid enthalpy
    ( *hvec )[0] = 1.8608e6;
    ( *pvec )[0] = 1.37895146e7;
    expected     = 609.133190;
    temperatureProperty->evalv( T, {}, hpArgs );
    if ( !AMP::Utilities::approx_equal( T[0], expected, tol ) ) {
        AMP::pout << "Extended temperature test 4 failed." << std::endl;
        AMP::pout << "Expected " << expected << ", computed " << T[0] << std::endl;
        extra_passed = false;
    }

    // h = 1400 Btu/lbm, p = 2000 psi -- Eq. 6c
    ( *hvec )[0] = 3.2564e6;
    ( *pvec )[0] = 1.37895146e7;
    expected     = 748.333926;
    temperatureProperty->evalv( T, {}, hpArgs );
    if ( !AMP::Utilities::approx_equal( T[0], expected, tol ) ) {
        AMP::pout << "Extended temperature test 5 failed." << std::endl;
        AMP::pout << "Expected " << expected << ", computed " << T[0] << std::endl;
        extra_passed = false;
    }

    // Test TemperatureProp
    std::vector<double> V( 1 );

    // h=200 Btu/lbm, p = 4000 psi -- Eq. 8a
    ( *hvec )[0] = 4.652e5;
    ( *pvec )[0] = 2.75790291e7;
    expected     = 1.03475993e-3;
    volumeProperty->evalv( V, {}, hpArgs );
    if ( !AMP::Utilities::approx_equal( V[0], expected, tol ) ) {
        AMP::pout << "Extended specific volume test 1 failed." << std::endl;
        AMP::pout << "Expected " << expected << ", computed " << V[0] << std::endl;
        extra_passed = false;
    }

    // h=600 Btu/lbm, p = 4000 psi -- Eq. 8b
    ( *hvec )[0] = 1.3956e6;
    ( *pvec )[0] = 2.75790291e7;
    expected     = 1.37847268e-3;
    volumeProperty->evalv( V, {}, hpArgs );
    if ( !AMP::Utilities::approx_equal( V[0], expected, tol ) ) {
        AMP::pout << "Extended specific volume test 2 failed." << std::endl;
        AMP::pout << "Expected " << expected << ", computed " << V[0] << std::endl;
        extra_passed = false;
    }

    // h=950 Btu/lbm, p = 4000 psi -- Eq. 8d
    ( *hvec )[0] = 2.2097e6;
    ( *pvec )[0] = 2.75790291e7;
    expected     = 3.18909350e-3;
    volumeProperty->evalv( V, {}, hpArgs );
    if ( !AMP::Utilities::approx_equal( V[0], expected, tol ) ) {
        AMP::pout << "Extended specific volume test 3 failed." << std::endl;
        AMP::pout << "Expected " << expected << ", computed " << V[0] << std::endl;
        extra_passed = false;
    }

    // h=1400 Btu/lbm, p = 4000 psi -- Eq. 8c
    ( *hvec )[0] = 3.2564e6;
    ( *pvec )[0] = 2.75790291e7;
    expected     = 1.07502487e-2;
    volumeProperty->evalv( V, {}, hpArgs );
    if ( !AMP::Utilities::approx_equal( V[0], expected, tol ) ) {
        AMP::pout << "Extended specific volume test 4 failed." << std::endl;
        AMP::pout << "Expected " << expected << ", computed " << V[0] << std::endl;
        extra_passed = false;
    }

    // h=200 Btu/lbm, p = 2000 psi -- Eq. 8a
    ( *hvec )[0] = 4.652e5;
    ( *pvec )[0] = 1.37895146e7;
    expected     = 1.04365756e-3;
    volumeProperty->evalv( V, {}, hpArgs );
    if ( !AMP::Utilities::approx_equal( V[0], expected, tol ) ) {
        AMP::pout << "Extended specific volume test 5 failed." << std::endl;
        AMP::pout << "Expected " << expected << ", computed " << V[0] << std::endl;
        extra_passed = false;
    }

    // h=500 Btu/lbm, p = 2000 psi -- Eq. 8b
    ( *hvec )[0] = 1.163e6;
    ( *pvec )[0] = 1.37895146e7;
    expected     = 1.267956e-3;
    volumeProperty->evalv( V, {}, hpArgs );
    if ( !AMP::Utilities::approx_equal( V[0], expected, tol ) ) {
        AMP::pout << "Extended specific volume test 6 failed." << std::endl;
        AMP::pout << "Expected " << expected << ", computed " << V[0] << std::endl;
        extra_passed = false;
    }

    // h=800 Btu/lbm, p = 2000 psi -- Eq. 8b and 8c at saturated values
    ( *hvec )[0] = 1.8608e6;
    ( *pvec )[0] = 1.37895146e7;
    expected     = 4.45266612e-3;
    volumeProperty->evalv( V, {}, hpArgs );
    if ( !AMP::Utilities::approx_equal( V[0], expected, tol ) ) {
        AMP::pout << "Extended specific volume test 7 failed." << std::endl;
        AMP::pout << "Expected " << expected << ", computed " << V[0] << std::endl;
        extra_passed = false;
    }

    // h=1400 Btu/lbm, p = 2000 psi -- Eq. 8c
    ( *hvec )[0] = 3.2564e6;
    ( *pvec )[0] = 1.37895146e7;
    expected     = 2.16299275e-2;
    volumeProperty->evalv( V, {}, hpArgs );
    if ( !AMP::Utilities::approx_equal( V[0], expected, tol ) ) {
        AMP::pout << "Extended specific volume test 8 failed." << std::endl;
        AMP::pout << "Expected " << expected << ", computed " << V[0] << std::endl;
        extra_passed = false;
    }

    // h=1040 Btu/lbm, p = 3000 psi -- Eq. 8d
    ( *hvec )[0] = 2.41904e6;
    ( *pvec )[0] = 2.06842718e7;
    expected     = 5.7290921e-3;
    volumeProperty->evalv( V, {}, hpArgs );
    if ( !AMP::Utilities::approx_equal( V[0], expected, tol ) ) {
        AMP::pout << "Extended specific volume test 9 failed." << std::endl;
        AMP::pout << "Expected " << expected << ", computed " << V[0] << std::endl;
        extra_passed = false;
    }

    if ( extra_passed )
        ut.passes( "Extended water library tests." );
    else
        ut.failure( "Extended water library tests." );

    ut.report();

    int num_failed = ut.NumFailGlobal();
    AMP::AMPManager::shutdown();
    return num_failed;
}
