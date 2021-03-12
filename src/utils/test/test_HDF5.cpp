// This file tests the HDF5 interfaces
#include <algorithm>
#include <array>
#include <cmath>
#include <complex>
#include <cstdio>
#include <cstring>
#include <iostream>
#include <string>
#include <vector>

#include "AMP/utils/AMPManager.h"
#include "AMP/utils/Array.h"
#include "AMP/utils/HDF5_Class.h"
#include "AMP/utils/HDF5_IO.h"
#include "AMP/utils/HDF5_IO.hpp"
#include "AMP/utils/UnitTest.h"
#include "AMP/utils/Utilities.h"


// Structure to contain variable types to write to hdf5
struct data_struct {
public: // Functions
    data_struct() : c( 0 ), i( 0 ), f( 0 ), d( 0 ) {}
    void initialize()
    {
        c             = 'g';
        i             = 3;
        f             = 3.1459;
        d             = 3.141592653589793;
        s             = "string";
        complex       = std::complex<double>( 2, 3 );
        vec_i         = { 1, 2, 3 };
        vec_d         = { 1.1, 2.2, 3.3 };
        vec_s         = { "a", "ab", "abc" };
        vec_complex   = { std::complex<double>( 2, 3 ), std::complex<double>( 4, 5 ) };
        array_i       = { 1, 2, 3 };
        array_d       = { 1.1, 2.2, 3.3 };
        array_s       = { "a", "ab", "abc" };
        array_complex = { std::complex<double>( 0, 1 ),
                          std::complex<double>( 2, 3 ),
                          std::complex<double>( 4, 5 ) };
        Array_i       = AMP::Array<int>( { 1, 2, 3 } );
        Array_d       = AMP::Array<double>( { 1.1, 2.2, 3.3 } );
        Array_s       = AMP::Array<std::string>( { "a", "ab", "abc" } );
        Array_complex = AMP::Array<std::complex<double>>(
            { std::complex<double>( 2, 3 ), std::complex<double>( 4, 5 ) } );
    }
    void clear() { *this = data_struct(); }

public: // Data members
    char c;
    int i;
    float f;
    double d;
    std::string s;
    std::complex<double> complex;
    std::vector<int> vec_i;
    std::vector<double> vec_d;
    std::vector<std::string> vec_s;
    std::vector<std::complex<double>> vec_complex;
    std::array<int, 3> array_i;
    std::array<double, 3> array_d;
    std::array<std::string, 3> array_s;
    std::array<std::complex<double>, 3> array_complex;
    AMP::Array<int> Array_i;
    AMP::Array<double> Array_d;
    AMP::Array<std::string> Array_s;
    AMP::Array<std::complex<double>> Array_complex;
};


// Write the data to HDF5
void writeHDF5( AMP::hid_t fid, const data_struct &data )
{
    AMP::writeHDF5( fid, "char", data.c );
    AMP::writeHDF5( fid, "int", data.i );
    AMP::writeHDF5( fid, "float", data.f );
    AMP::writeHDF5( fid, "double", data.d );
    AMP::writeHDF5( fid, "string", data.s );
    AMP::writeHDF5( fid, "complex", data.complex );
    AMP::writeHDF5( fid, "vector<int>", data.vec_i );
    AMP::writeHDF5( fid, "vector<double>", data.vec_d );
    AMP::writeHDF5( fid, "vector<string>", data.vec_s );
    AMP::writeHDF5( fid, "vector<complex>", data.vec_complex );
    AMP::writeHDF5( fid, "array<int>", data.array_i );
    AMP::writeHDF5( fid, "array<double>", data.array_d );
    AMP::writeHDF5( fid, "array<string>", data.array_s );
    AMP::writeHDF5( fid, "array<complex>", data.array_complex );
    AMP::writeHDF5( fid, "Array<int>", data.Array_i );
    AMP::writeHDF5( fid, "Array<double>", data.Array_d );
    AMP::writeHDF5( fid, "Array<string>", data.Array_s );
    AMP::writeHDF5( fid, "Array<complex>", data.Array_complex );
}


// Read the data from HDF5
void readHDF5( AMP::hid_t fid, data_struct &data )
{
    AMP::readHDF5( fid, "char", data.c );
    AMP::readHDF5( fid, "int", data.i );
    AMP::readHDF5( fid, "float", data.f );
    AMP::readHDF5( fid, "double", data.d );
    AMP::readHDF5( fid, "string", data.s );
    AMP::readHDF5( fid, "complex", data.complex );
    AMP::readHDF5( fid, "vector<int>", data.vec_i );
    AMP::readHDF5( fid, "vector<double>", data.vec_d );
    AMP::readHDF5( fid, "vector<string>", data.vec_s );
    AMP::readHDF5( fid, "vector<complex>", data.vec_complex );
    AMP::readHDF5( fid, "array<int>", data.array_i );
    AMP::readHDF5( fid, "array<double>", data.array_d );
    AMP::readHDF5( fid, "array<string>", data.array_s );
    AMP::readHDF5( fid, "array<complex>", data.array_complex );
    AMP::readHDF5( fid, "Array<int>", data.Array_i );
    AMP::readHDF5( fid, "Array<double>", data.Array_d );
    AMP::readHDF5( fid, "Array<string>", data.Array_s );
    AMP::readHDF5( fid, "Array<complex>", data.Array_complex );
}


// Read the data from HDF5
template<class TYPE>
void readScalar( AMP::HDF5data &ptr, const std::string &name, TYPE &x )
{
    AMP::Array<TYPE> data;
    auto ptr2 = ptr.getData( 0, name );
    ptr2->getData( data );
    AMP_ASSERT( data.length() == 1 );
    x = data( 0 );
}
template<class TYPE>
void readVector( AMP::HDF5data &ptr, const std::string &name, std::vector<TYPE> &x )
{
    AMP::Array<TYPE> data;
    auto ptr2 = ptr.getData( 0, name );
    ptr2->getData( data );
    x = std::vector<TYPE>( data.begin(), data.end() );
}
template<class TYPE, std::size_t N>
void readArray( AMP::HDF5data &ptr, const std::string &name, std::array<TYPE, N> &x )
{
    AMP::Array<TYPE> data;
    auto ptr2 = ptr.getData( 0, name );
    ptr2->getData( data );
    AMP_ASSERT( data.length() == N );
    for ( size_t i = 0; i < N; i++ )
        x[i] = data( i );
}
template<class TYPE>
void readArray( AMP::HDF5data &ptr, const std::string &name, AMP::Array<TYPE> &x )
{
    auto ptr2 = ptr.getData( 0, name );
    ptr2->getData( x );
}
void readHDF52( AMP::hid_t fid, data_struct &data )
{
    auto ptr = AMP::readHDF5( fid, "/" );
    readScalar( *ptr, "char", data.c );
    readScalar( *ptr, "int", data.i );
    readScalar( *ptr, "float", data.f );
    readScalar( *ptr, "double", data.d );
    readScalar( *ptr, "string", data.s );
    readVector( *ptr, "vector<int>", data.vec_i );
    readVector( *ptr, "vector<double>", data.vec_d );
    readVector( *ptr, "vector<string>", data.vec_s );
    readArray( *ptr, "array<int>", data.array_i );
    readArray( *ptr, "array<double>", data.array_d );
    readArray( *ptr, "array<string>", data.array_s );
    readArray( *ptr, "Array<int>", data.Array_i );
    readArray( *ptr, "Array<double>", data.Array_d );
    readArray( *ptr, "Array<string>", data.Array_s );
    // Complex is not well supported yet
    AMP::readHDF5( fid, "complex", data.complex );
    AMP::readHDF5( fid, "vector<complex>", data.vec_complex );
    AMP::readHDF5( fid, "array<complex>", data.array_complex );
    AMP::readHDF5( fid, "Array<complex>", data.Array_complex );
}


// Check the results
template<class TYPE>
void check( const std::string &name, const TYPE &x, const TYPE &y, AMP::UnitTest &ut )
{
    if ( x == y )
        ut.passes( name );
    else
        ut.failure( name );
}
void checkResults( const std::string &str,
                   const data_struct &data,
                   const data_struct &data2,
                   AMP::UnitTest &ut )
{
    check( str + "char", data.c, data2.c, ut );
    check( str + "int", data.i, data2.i, ut );
    check( str + "float", data.f, data2.f, ut );
    check( str + "double", data.d, data2.d, ut );
    check( str + "string", data.s, data2.s, ut );
    check( str + "complex", data.complex, data2.complex, ut );
    check( str + "vector<int>", data.vec_i, data2.vec_i, ut );
    check( str + "vector<double>", data.vec_d, data2.vec_d, ut );
    check( str + "vector<string>", data.vec_s, data2.vec_s, ut );
    check( str + "vector<complex>", data.vec_complex, data2.vec_complex, ut );
    check( str + "array<int>", data.array_i, data2.array_i, ut );
    check( str + "array<double>", data.array_d, data2.array_d, ut );
    check( str + "array<string>", data.array_s, data2.array_s, ut );
    check( str + "array<complex>", data.array_complex, data2.array_complex, ut );
    check( str + "Array<int>", data.Array_i, data2.Array_i, ut );
    check( str + "Array<double>", data.Array_d, data2.Array_d, ut );
    check( str + "Array<string>", data.Array_s, data2.Array_s, ut );
    check( str + "Array<complex>", data.Array_complex, data2.Array_complex, ut );
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

    // Write variables to HDF5
    data_struct data;
    data.initialize();
    auto fid = AMP::openHDF5( "test_HDF5.hdf5", "w" );
    writeHDF5( fid, data );
    AMP::closeHDF5( fid );

    // Read the variables from HDF5
    data_struct data2;
    fid = AMP::openHDF5( "test_HDF5.hdf5", "r" );
    readHDF5( fid, data2 );
    AMP::closeHDF5( fid );
    checkResults( "readHDF5 (1): ", data, data2, ut );

    // Read the data using the class interface
    data_struct data3;
    fid = AMP::openHDF5( "test_HDF5.hdf5", "r" );
    readHDF52( fid, data3 );
    AMP::closeHDF5( fid );
    checkResults( "readHDF5 (2): ", data, data3, ut );

    // Return
    data.clear();
    data2.clear();
    data3.clear();
    int N_errors = ut.NumFailGlobal();
    ut.report();
    ut.reset();
    AMP::AMPManager::shutdown();
    return N_errors;
}
