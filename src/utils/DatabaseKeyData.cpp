#include "AMP/utils/Array.h"
#include "AMP/utils/Array.hpp"
#include "AMP/utils/Database.h"
#include "AMP/utils/Utilities.h"


#include <cmath>
#include <complex>
#include <cstdlib>
#include <cstring>
#include <iomanip>
#include <limits>
#include <sstream>
#include <string>


namespace AMP {


/********************************************************************
 * Helper function to perform data conversion                        *
 ********************************************************************/
template<class TYPE>
static constexpr TYPE getTol()
{
    if constexpr ( std::is_same<TYPE, float>::value ) {
        return 1e-6;
    } else if constexpr ( std::is_same<TYPE, double>::value ) {
        return 1e-12;
    } else {
        return 0;
    }
}
template<class TYPE>
static constexpr TYPE abs( const TYPE &x )
{
    if constexpr ( std::is_signed<TYPE>::value ) {
        return x < 0 ? -x : x;
    } else {
        return x;
    }
}
template<class TYPE1, class TYPE2>
Array<TYPE2> convert( const Array<TYPE1> &x )
{
    using AMP::Utilities::type_name;
    if constexpr ( std::is_same<TYPE1, TYPE2>::value ) {
        return x;
    } else if constexpr ( std::is_arithmetic<TYPE1>::value && std::is_arithmetic<TYPE2>::value ) {
        Array<TYPE2> y( x.size() );
        y.fill( 0 );
        bool pass = true;
        TYPE1 tol = getTol<TYPE1>();
        for ( size_t i = 0; i < x.length(); i++ ) {
            y( i ) = static_cast<TYPE2>( x( i ) );
            pass   = pass && abs( static_cast<TYPE1>( y( i ) - x( i ) ) ) <= tol * x( i );
        }
        if ( !pass ) {
            std::string msg = "Converting " + std::string( type_name<TYPE1>() ) + "-" +
                              std::string( type_name<TYPE2>() ) + " results in loss of precision";
            AMP_WARNING( msg );
        }
        return y;
    } else {
        std::string msg = "Invalid conversion: " + std::string( type_name<TYPE1>() ) + "-" +
                          std::string( type_name<TYPE2>() ) + " results in loss of precision";
        throw std::logic_error( msg );
    }
}


/********************************************************************
 * Scale data                                                        *
 ********************************************************************/
template<class TYPE>
void scaleData( TYPE &data, double factor )
{
    if constexpr ( std::is_same<TYPE, bool>::value ) {
        throw std::logic_error( "Unable to scale bool" );
    } else if constexpr ( std::is_same<TYPE, std::complex<float>>::value ) {
        data = static_cast<float>( factor ) * data;
    } else if constexpr ( std::is_same<TYPE, std::complex<double>>::value ) {
        data = factor * data;
    } else if constexpr ( std::is_arithmetic<TYPE>::value ) {
        data = static_cast<TYPE>( factor ) * data;
    } else {
        NULL_USE( factor );
        using AMP::Utilities::type_name;
        throw std::logic_error( "Unable to scale " + std::string( type_name<TYPE>() ) );
    }
}
template<class TYPE>
void scaleData( Array<TYPE> &data, double factor )
{
    if constexpr ( std::is_same<TYPE, bool>::value ) {
        throw std::logic_error( "Unable to scale bool" );
    } else if constexpr ( std::is_arithmetic<TYPE>::value ) {
        data.scale( factor );
    } else {
        NULL_USE( factor );
        using AMP::Utilities::type_name;
        throw std::logic_error( "Unable to scale " + std::string( type_name<TYPE>() ) );
    }
}


/********************************************************************
 * KeyDataScalar::operator==                                         *
 ********************************************************************/
template<class TYPE>
static inline bool compare( const TYPE &x, const TYPE &y );
template<>
inline bool compare( const double &x, const double &y )
{
    bool test = x == y;
    test      = test || fabs( x - y ) <= 1e-12 * fabs( x + y );
    test      = test || ( ( x != x ) && ( y != y ) );
    return test;
}
template<>
inline bool compare( const float &x, const float &y )
{
    bool test = x == y;
    test      = test || fabs( x - y ) <= 1e-7 * fabs( x + y );
    test      = test || ( ( x != x ) && ( y != y ) );
    return test;
}
template<>
inline bool compare( const std::complex<double> &x, const std::complex<double> &y )
{
    bool test = x == y;
    test      = test || std::abs( x - y ) <= 1e-12 * std::abs( x + y );
    test      = test || ( ( x != x ) && ( y != y ) );
    return test;
}
template<>
inline bool compare( const std::complex<float> &x, const std::complex<float> &y )
{
    bool test = x == y;
    test      = test || std::abs( x - y ) <= 1e-7 * std::abs( x + y );
    test      = test || ( ( x != x ) && ( y != y ) );
    return test;
}
template<class TYPE>
static inline bool compare( const TYPE &x, const TYPE &y )
{
    return x == y;
}
template<class TYPE>
bool KeyDataScalar<TYPE>::operator==( const KeyData &rhs ) const
{
    auto tmp1 = dynamic_cast<const KeyDataScalar<TYPE> *>( &rhs );
    auto tmp2 = dynamic_cast<const KeyDataArray<TYPE> *>( &rhs );
    if ( tmp1 ) {
        return compare( d_data, tmp1->d_data );
    } else if ( tmp2 ) {
        if ( tmp2->get().size() != 1 )
            return false;
        return compare( d_data, tmp2->get()( 0 ) );
    } else if ( ( is_floating_point() || is_integral() ) &&
                ( rhs.is_floating_point() || rhs.is_integral() ) ) {
        auto data1 = convertToDouble();
        auto data2 = rhs.convertToDouble();
        if ( data1.size() != data2.size() )
            return false;
        bool test = true;
        for ( size_t i = 0; i < data1.length(); i++ )
            test = test && compare( data1( i ), data2( i ) );
        return test;
    }
    return false;
}
template<class TYPE>
bool KeyDataArray<TYPE>::operator==( const KeyData &rhs ) const
{
    auto tmp1 = dynamic_cast<const KeyDataScalar<TYPE> *>( &rhs );
    auto tmp2 = dynamic_cast<const KeyDataArray<TYPE> *>( &rhs );
    if ( tmp1 ) {
        if ( d_data.size() != 1 )
            return false;
        return compare( d_data( 0 ), tmp1->get() );
    } else if ( tmp2 ) {
        return compare( d_data, tmp2->get() );
    } else if ( ( is_floating_point() || is_integral() ) &&
                ( rhs.is_floating_point() || rhs.is_integral() ) ) {
        auto data1 = convertToDouble();
        auto data2 = rhs.convertToDouble();
        if ( data1.size() != data2.size() )
            return false;
        bool test = true;
        for ( size_t i = 0; i < data1.length(); i++ )
            test = test && compare( data1( i ), data2( i ) );
        return test;
    }
    return false;
}


/********************************************************************
 * DatabaseBox                                                       *
 ********************************************************************/
DatabaseBox::DatabaseBox() : d_dim( 0 )
{
    d_lower.fill( 0 );
    d_upper.fill( 0 );
}
DatabaseBox::DatabaseBox( int dim, const int *lower, const int *upper ) : d_dim( dim )
{
    AMP_ASSERT( dim <= 5 );
    d_lower.fill( 0 );
    d_upper.fill( 0 );
    for ( int d = 0; d < dim; d++ ) {
        d_lower[d] = lower[d];
        d_upper[d] = upper[d];
    }
}
DatabaseBox::DatabaseBox( const std::string_view &str ) : d_dim( 0 )
{
    d_lower.fill( 0 );
    d_upper.fill( 0 );
    size_t j     = 0;
    size_t k     = 0;
    char tmp[10] = { 0 };
    for ( size_t i = 0; i < str.length(); i++ ) {
        if ( str[i] == ' ' || str[i] == '(' || str[i] == '[' || str[i] == ']' ) {
            continue;
        } else if ( str[i] == ',' || str[i] == ')' ) {
            if ( str[i] == ')' && d_dim == 0 ) {
                d_dim = k + 1;
                while ( str[i] != ',' && i + 1 < str.length() ) {
                    ++i;
                }
            }
            tmp[j] = 0;
            j      = 0;
            int v  = atoi( tmp );
            if ( k < d_dim || d_dim == 0 )
                d_lower[k++] = v;
            else
                d_upper[( k++ ) - d_dim] = v;
        } else {
            tmp[j++] = str[i];
        }
    }
}
uint8_t &DatabaseBox::dim() { return d_dim; }
uint8_t DatabaseBox::dim() const { return d_dim; }
bool DatabaseBox::empty() const
{
    if ( d_dim == 0 )
        return true;
    for ( int d = 0; d < d_dim; d++ ) {
        if ( d_upper[d] < d_lower[d] )
            return true;
    }
    return false;
}
int &DatabaseBox::lower( uint8_t d )
{
    AMP_ASSERT( d_dim <= 5 && d < d_dim );
    return d_lower[d];
}
int DatabaseBox::lower( uint8_t d ) const
{
    AMP_ASSERT( d_dim <= 5 && d < d_dim );
    return d_lower[d];
}
int &DatabaseBox::upper( uint8_t d )
{
    AMP_ASSERT( d_dim <= 5 && d < d_dim );
    return d_upper[d];
}
int DatabaseBox::upper( uint8_t d ) const
{
    AMP_ASSERT( d_dim <= 5 && d < d_dim );
    return d_upper[d];
}
bool DatabaseBox::operator==( const DatabaseBox &box ) const
{
    bool equal = d_dim == box.d_dim;
    for ( int d = 0; d < d_dim; d++ ) {
        equal = equal && d_lower[d] == box.d_lower[d];
        equal = equal && d_upper[d] == box.d_upper[d];
    }
    return equal;
}
std::ostream &operator<<( std::ostream &out, const DatabaseBox &box )
{
    if ( box.empty() )
        out << "[(),()]";
    out << "[(" << box.lower( 0 );
    for ( int d = 1; d < box.dim(); d++ )
        out << "," << box.lower( d );
    out << "),(" << box.upper( 0 );
    for ( int d = 1; d < box.dim(); d++ )
        out << "," << box.upper( d );
    out << ")";
    return out;
}


/********************************************************************
 * Explicit instantiations                                           *
 ********************************************************************/
#define instantiate( FUN )       \
    FUN( bool );                 \
    FUN( char );                 \
    FUN( int8_t );               \
    FUN( int16_t );              \
    FUN( int32_t );              \
    FUN( int64_t );              \
    FUN( uint8_t );              \
    FUN( uint16_t );             \
    FUN( uint32_t );             \
    FUN( uint64_t );             \
    FUN( float );                \
    FUN( double );               \
    FUN( long double );          \
    FUN( std::complex<float> );  \
    FUN( std::complex<double> ); \
    FUN( std::string );          \
    FUN( std::_Bit_reference );  \
    FUN( AMP::DatabaseBox )
#define instantiateConvert( TYPE )                                             \
    template AMP::Array<TYPE> convert( const Array<bool> &x );                 \
    template AMP::Array<TYPE> convert( const Array<char> &x );                 \
    template AMP::Array<TYPE> convert( const Array<int8_t> &x );               \
    template AMP::Array<TYPE> convert( const Array<int16_t> &x );              \
    template AMP::Array<TYPE> convert( const Array<int32_t> &x );              \
    template AMP::Array<TYPE> convert( const Array<int64_t> &x );              \
    template AMP::Array<TYPE> convert( const Array<uint8_t> &x );              \
    template AMP::Array<TYPE> convert( const Array<uint16_t> &x );             \
    template AMP::Array<TYPE> convert( const Array<uint32_t> &x );             \
    template AMP::Array<TYPE> convert( const Array<uint64_t> &x );             \
    template AMP::Array<TYPE> convert( const Array<float> &x );                \
    template AMP::Array<TYPE> convert( const Array<double> &x );               \
    template AMP::Array<TYPE> convert( const Array<long double> &x );          \
    template AMP::Array<TYPE> convert( const Array<std::complex<float>> &x );  \
    template AMP::Array<TYPE> convert( const Array<std::complex<double>> &x ); \
    template AMP::Array<TYPE> convert( const Array<std::string> &x );          \
    template AMP::Array<TYPE> convert( const Array<std::_Bit_reference> &x );  \
    template AMP::Array<TYPE> convert( const Array<AMP::DatabaseBox> &x )
#define instantiateScaleData( TYPE )                             \
    template void scaleData<TYPE>( TYPE & data, double factor ); \
    template void scaleData<TYPE>( AMP::Array<TYPE> & data, double factor )
#define instantiateKeyDataScalar( TYPE ) template class KeyDataScalar<TYPE>
#define instantiateKeyDataArray( TYPE ) template class KeyDataArray<TYPE>
instantiate( instantiateConvert );       // convert
instantiate( instantiateScaleData );     // scaleData
instantiate( instantiateKeyDataScalar ); // KeyDataScalar
instantiate( instantiateKeyDataArray );  // KeyDataArray


/********************************************************
 *  Explicit instantiations of Array<DatabaseBox>        *
 ********************************************************/
template Array<DatabaseBox, FunctionTable>::Array();
template Array<DatabaseBox, FunctionTable>::~Array();
template Array<DatabaseBox, FunctionTable>::Array( size_t );
template Array<DatabaseBox, FunctionTable>::Array( size_t, size_t );
template Array<DatabaseBox, FunctionTable>::Array( size_t, size_t, size_t );
template Array<DatabaseBox, FunctionTable>::Array( size_t, size_t, size_t, size_t );
template Array<DatabaseBox, FunctionTable>::Array( size_t, size_t, size_t, size_t, size_t );
template Array<DatabaseBox, FunctionTable>::Array( const Array<DatabaseBox, FunctionTable> & );
template Array<DatabaseBox, FunctionTable>::Array( Array<DatabaseBox, FunctionTable> && );
template Array<DatabaseBox, FunctionTable> &
Array<DatabaseBox, FunctionTable>::operator=( const Array<DatabaseBox, FunctionTable> & );
template Array<DatabaseBox, FunctionTable> &
Array<DatabaseBox, FunctionTable>::operator=( Array<DatabaseBox, FunctionTable> && );
template Array<DatabaseBox, FunctionTable> &
Array<DatabaseBox, FunctionTable>::operator=( const std::vector<DatabaseBox> & );
template void Array<DatabaseBox, FunctionTable>::clear();
template void
Array<DatabaseBox, FunctionTable>::viewRaw( ArraySize const &, DatabaseBox *, bool, bool );
template bool
Array<DatabaseBox, FunctionTable>::operator==( Array<DatabaseBox, FunctionTable> const & ) const;
template void Array<DatabaseBox, FunctionTable>::resize( ArraySize const & );


} // namespace AMP
