/*
 * test_U02_MSRZC_09.cc
 *
 *  Created on: Feb 8, 2010
 *      Author: gad
 */


#include <iostream>
#include <string>
#include <valarray>
#include <vector>

#include "AMP/materials/Material.h"
#include "AMP/utils/AMPManager.h"
#include "AMP/utils/UnitTest.h"
#include "AMP/utils/Utilities.h"

using namespace std;

int main( int argc, char **argv )
{
    AMP::AMPManager::startup( argc, argv );
    AMP::UnitTest ut;

    using namespace AMP::Materials;

    bool good = true;

    // test constructors
    AMP::shared_ptr<AMP::Materials::Material> mat =
        AMP::voodoo::Factory<AMP::Materials::Material>::instance().create( "UO2_MSRZC_09" );
    auto prop = mat->property( "ThermalConductivity" );

    // test property accessors
    string tcname = prop->get_name();
    string tcsorc = prop->get_source();
    good          = good && tcname == string( "UO2_MSRZC_09_ThermalConductivity" );
    std::cout << "thermal conductivity name is " << tcname << "\n";
    std::cout << "thermal conductivity source is " << tcsorc << "\n";
    vector<string> args = prop->get_arguments();
    good                = good && args[0] == "temperature";
    good                = good && args[1] == "concentration";
    std::cout << "arguments are " << args[0] << " " << args[1] << "\n";
    unsigned int nargs = prop->get_number_arguments();
    good               = good && nargs == 2;

    // test material accessors, all arguments present
    size_t n = 10;
    AMP::shared_ptr<std::vector<double>> tv( new std::vector<double>( n ) );
    AMP::shared_ptr<std::vector<double>> uv( new std::vector<double>( n ) );
    vector<double> tcv( n );
    for ( size_t i = 0; i < n; i++ ) {
        ( *tv )[i] = 563.4;
        ( *uv )[i] = .05;
    }

    // Block for temporary variables
    {
        std::map<std::string, AMP::shared_ptr<std::vector<double>>> argMap;
        argMap.insert( std::make_pair( "temperature", tv ) );
        argMap.insert( std::make_pair( "concentration", uv ) );

        std::vector<double> tcv_mat( tcv );
        mat->property( "ThermalConductivity" )->evalv( tcv_mat, argMap );
        prop->evalv( tcv, argMap );
        for ( size_t i = 0; i < n; i++ ) {
            good = good && AMP::Utilities::approx_equal( tcv[i], tcv_mat[i] );
        }
        for ( size_t i = 0; i < n; i++ ) {
            good = good && AMP::Utilities::approx_equal( tcv[0], tcv[i] );
        }
    }

    // test material accessors, one argument present
    std::vector<double> defaults( 2 );
    defaults[0] = 563.4;
    defaults[1] = 0.05;
    prop->set_defaults( defaults );


    // Block for temporary variables
    {
        std::map<std::string, AMP::shared_ptr<std::vector<double>>> argMap;
        argMap.insert( std::make_pair( "temperature", tv ) );
        std::vector<double> tcv_def( tcv );
        prop->evalv( tcv_def, argMap );
        for ( size_t i = 0; i < n; i++ ) {
            good = good && AMP::Utilities::approx_equal( tcv[i], tcv_def[i] );
        }
    }

    // test material accessors, no arguments present
    {
        std::map<std::string, AMP::shared_ptr<std::vector<double>>> argMap;
        std::vector<double> tcv_def( tcv );
        prop->evalv( tcv_def, argMap );
        for ( size_t i = 0; i < n; i++ ) {
            good = good && AMP::Utilities::approx_equal( tcv[i], tcv_def[i] );
        }
    }

    if ( good )
        ut.passes( "basic tests of Material" );
    else
        ut.failure( "basic tests of Material" );

    ut.report();

    int num_failed = ut.NumFailGlobal();
    AMP::AMPManager::shutdown();
    return num_failed;
}
