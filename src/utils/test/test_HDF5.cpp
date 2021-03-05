// This file tests the HDF5 interfaces
#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <iostream>
#include <random>
#include <string>

#include "AMP/utils/AMPManager.h"
#include "AMP/utils/HDF5_IO.h"
#include "AMP/utils/HDF5_IO.hpp"
#include "AMP/utils/UnitTest.h"
#include "AMP/utils/Utilities.h"


// Test writing/reading a variable
template<class T>
void test_variable( const T &x, AMP::UnitTest &ut )
{
    remove( "test_HDF5.hdf5" );
    auto fid = AMP::openHDF5( "test_HDF5.hdf5", "w" );
    AMP::writeHDF5( fid, "test", x );
    AMP::closeHDF5( fid );
    fid = AMP::openHDF5( "test_HDF5.hdf5", "r" );
    T y;
    AMP::readHDF5( fid, "test", y );
    AMP::closeHDF5( fid );
    if ( x == y )
        ut.passes( typeid( x ).name() );
    else
        ut.failure( typeid( x ).name() );
}


// Create a random string
std::string random_string( size_t length )
{
    static const char alphanums[] = "0123456789"
                                    "abcdefghijklmnopqrstuvwxyz"
                                    "ABCDEFGHIJKLMNOPQRSTUVWXYZ";
    static std::mt19937 rg( std::chrono::system_clock::now().time_since_epoch().count() );
    static std::uniform_int_distribution<> pick( 0, sizeof( alphanums ) - 1 );
    std::string s;
    s.reserve( length );
    while ( length-- )
        s += alphanums[pick( rg )];
    return s;
}


// Main
int main( int argc, char *argv[] )
{
    AMP::AMPManager::startup( argc, argv );
    AMP::UnitTest ut;

    // Test catching an hdf5 error
    try {
        AMP::openHDF5( "", "w" );
        ut.failure( "Failed to catch HDF5 error" );
    } catch ( ... ) {
        ut.passes( "Caught HDF5 error" );
    }

    // Test some simple scalar types
    test_variable<char>( 'g', ut );
    test_variable<int>( 3, ut );
    test_variable<float>( 3.1459, ut );
    test_variable<double>( 3.141592653589793, ut );
    test_variable<std::string>( random_string( 25 ), ut );
    test_variable<std::complex<double>>( std::complex<double>( 2, 3 ), ut );

    // Test std::vector
    test_variable<std::vector<double>>( { 1, 2, 3 }, ut );
    test_variable<std::vector<std::string>>( { "a", "ab", "abc" }, ut );

    // Test std::array
    std::array<double, 3> x2               = { { 1, 2, 3 } };
    std::array<std::complex<double>, 2> x3 = { { std::complex<double>( 2, 3 ),
                                                 std::complex<double>( 4, 5 ) } };
    test_variable( x2, ut );
    test_variable( x3, ut );

    // Return
    int N_errors = ut.NumFailGlobal();
    ut.report();
    ut.reset();
    AMP::AMPManager::shutdown();
    return N_errors;
}
