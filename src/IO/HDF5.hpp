// This file contains helper functions and interfaces for reading/writing HDF5

#ifndef included_AMP_HDF5_hpp
#define included_AMP_HDF5_hpp

#include "AMP/AMP_TPLs.h"
#include "AMP/IO/HDF5.h"
#include "AMP/IO/HDF5_Class.h"
#include "AMP/utils/Array.h"
#include "AMP/utils/TypeTraits.h"
#include "AMP/utils/Utilities.h"
#include "AMP/utils/typeid.h"

#include <array>
#include <complex>
#include <memory>
#include <type_traits>
#include <vector>

#ifdef AMP_USE_HDF5


namespace AMP {


/******************************************************************
 * Default implementation                                          *
 ******************************************************************/
template<class T>
void writeHDF5Array( hid_t, const std::string_view &, const AMP::Array<T> & );
template<class T>
void readHDF5Array( hid_t, const std::string_view &, AMP::Array<T> & );
template<class T>
void writeHDF5Scalar( hid_t, const std::string_view &, const T & );
template<class T>
void readHDF5Scalar( hid_t, const std::string_view &, T & );
template<class TYPE>
void writeHDF5( hid_t fid, const std::string_view &name, const TYPE &x )
{
    NULL_USE( fid );
    if constexpr ( AMP::is_shared_ptr<TYPE>::value ) {
        // We are dealing with a std::shared_ptr
        writeHDF5( fid, name, *x );
    } else if constexpr ( AMP::is_vector<TYPE>::value ) {
        // We are dealing with a std::vector
        typedef decltype( *x.begin() ) TYPE2;
        typedef typename std::remove_reference<TYPE2>::type TYPE3;
        typedef typename std::remove_cv<TYPE3>::type TYPE4;
        AMP::Array<TYPE4> y;
        y.viewRaw( { x.size() }, const_cast<TYPE4 *>( x.data() ) );
        writeHDF5Array( fid, name, y );
    } else if constexpr ( std::is_array<TYPE>::value ) {
        // We are dealing with a std::array
        typedef decltype( *x.begin() ) TYPE2;
        typedef typename std::remove_reference<TYPE2>::type TYPE3;
        typedef typename std::remove_cv<TYPE3>::type TYPE4;
        AMP::Array<TYPE4> y;
        y.viewRaw( { x.size() }, const_cast<TYPE4 *>( x.data() ) );
        writeHDF5Array( fid, name, y );
    } else if constexpr ( AMP::is_Array<TYPE>::value ) {
        // We are dealing with an Array
        writeHDF5Array( fid, name, x );
    } else if constexpr ( std::is_same_v<TYPE, std::string> ) {
        // We are dealing with a std::string
        writeHDF5Scalar( fid, name, x );
    } else if constexpr ( std::is_same<TYPE, std::string_view>::value ||
                          std::is_same<TYPE, char *>::value ||
                          std::is_same<TYPE, const char *>::value ) {
        // We are dealing with a string or char array
        writeHDF5( fid, name, std::string( x ) );
    } else if constexpr ( AMP::has_size<TYPE>::value ) {
        // We are dealing with a container
        typedef decltype( *x.begin() ) TYPE2;
        typedef typename std::remove_reference<TYPE2>::type TYPE3;
        typedef typename std::remove_cv<TYPE3>::type TYPE4;
        std::vector<TYPE4> x2( x.begin(), x.end() );
        writeHDF5<std::vector<TYPE4>>( fid, name, x2 );
    } else {
        writeHDF5Scalar( fid, name, x );
    }
}
template<class TYPE>
void readHDF5( hid_t fid, const std::string_view &name, TYPE &x )
{
    NULL_USE( fid );
    if constexpr ( AMP::is_shared_ptr<TYPE>::value ) {
        // We are dealing with a std::shared_ptr
        readHDF5( fid, name, *x );
    } else if constexpr ( AMP::is_vector<TYPE>::value ) {
        // We are dealing with a std::vector
        typedef typename std::remove_reference<decltype( *x.begin() )>::type TYPE2;
        AMP::Array<TYPE2> y;
        readHDF5Array( fid, name, y );
        x.resize( y.length() );
        // Swap the elements in the arrays to use the move operator
        for ( size_t i = 0; i < x.size(); i++ )
            std::swap( x[i], y( i ) );
    } else if constexpr ( std::is_array<TYPE>::value ) {
        // We are dealing with a std::array
        typedef typename std::remove_reference<decltype( *x.begin() )>::type TYPE2;
        AMP::Array<TYPE2> y;
        readHDF5Array( fid, name, y );
        AMP_ASSERT( y.length() == x.size() );
        // Swap the elements in the arrays to use the move operator
        for ( size_t i = 0; i < x.size(); i++ )
            std::swap( x[i], y( i ) );
    } else if constexpr ( AMP::is_Array<TYPE>::value ) {
        // We are dealing with an Array
        readHDF5Array( fid, name, x );
    } else if constexpr ( std::is_same_v<TYPE, std::string> ) {
        // We are dealing with a std::string
        readHDF5Scalar( fid, name, x );
    } else if constexpr ( std::is_same<TYPE, std::string_view>::value ||
                          std::is_same<TYPE, char *>::value ||
                          std::is_same<TYPE, const char *>::value ) {
        // We are dealing with a string or char array
        throw std::logic_error(
            "Reading data into a string_view, char*, const char* is not supported" );
    } else if constexpr ( AMP::has_size<TYPE>::value ) {
        // We are dealing with a container
        typedef typename std::remove_reference<decltype( *x.begin() )>::type TYPE2;
        AMP::Array<TYPE2> y;
        readHDF5Array( fid, name, y );
        if ( x.size() == y.length() ) {
            auto it = x.begin();
            for ( size_t i = 0; i < y.length(); i++, ++it )
                *it = y( i );
        } else {
            throw std::logic_error( "Reading data into an arbitrary container is not finished" );
        }
    } else {
        readHDF5Scalar( fid, name, x );
    }
}


/************************************************************************
 * Helper function to get the size of an Array                           *
 * Note that HDF5 uses C ordered arrays so we need to flip the dimensions*
 ************************************************************************/
inline std::vector<hsize_t> arraySize( const AMP::ArraySize &s1 )
{
    int N = s1.ndim();
    std::vector<hsize_t> s2( std::max( N, 1 ), 0 );
    for ( int i = 0; i < N; i++ )
        s2[N - i - 1] = static_cast<hsize_t>( s1[i] );
    return s2;
}
template<class T>
inline std::vector<hsize_t> arraySize( const AMP::Array<T> &x )
{
    return arraySize( x.size() );
}
inline std::vector<size_t> convertSize( int N, const hsize_t *dims )
{
    if ( N == 0 )
        return std::vector<size_t>( 1, 1 );
    std::vector<size_t> size( N, 0 );
    for ( int i = 0; i < N; i++ )
        size[N - i - 1] = static_cast<size_t>( dims[i] );
    return size;
}


} // namespace AMP


#endif
#endif
