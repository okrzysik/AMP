#include "AMP/utils/Array.h"
#include "AMP/utils/Array.hpp"
#include "AMP/utils/Database.h"
#include "AMP/utils/Utilities.h"


#include <complex>
#include <cstdlib>
#include <cstring>
#include <iomanip>
#include <limits>
#include <sstream>
#include <string>


// clang-format off
namespace AMP {


/********************************************************************
 * Scale data                                                        *
 ********************************************************************/
#define scaleDataValid( TYPE )                            \
    template<>                                            \
    void scaleData( Array<TYPE>& data, double factor )    \
    {                                                     \
        data.scale( factor );                             \
    }                                                     \
    template<>                                            \
    void scaleData( TYPE& data, double factor )           \
    {                                                     \
        data = static_cast<TYPE>( factor * data );        \
    }
#define scaleDataInvalid( TYPE )                          \
    template<>                                            \
    void scaleData( Array<TYPE>&, double )                \
    {                                                     \
        throw std::logic_error( "Unable to scale data" ); \
    }                                                     \
    template<>                                            \
    void scaleData( TYPE&, double )                       \
    {                                                     \
        throw std::logic_error( "Unable to scale data" ); \
    }
scaleDataValid( int8_t )
scaleDataValid( int16_t )
scaleDataValid( int32_t )
scaleDataValid( int64_t )
scaleDataValid( uint8_t )
scaleDataValid( uint16_t )
scaleDataValid( uint32_t )
scaleDataValid( uint64_t )
scaleDataValid( float )
scaleDataValid( double )
scaleDataValid( long double )
scaleDataValid( std::complex<double> )
scaleDataInvalid( bool )
scaleDataInvalid( char )
scaleDataInvalid( std::string )
scaleDataInvalid( DatabaseBox )
template<>
void scaleData( Array<std::complex<float>>& data, double factor )
{
    data.scale( factor );
}
template<>
void scaleData( std::complex<float>& data, double factor )
{
    data = static_cast<float>( factor ) * data;
}


/********************************************************************
 * ConvertData                                                       *
 ********************************************************************/
#define convertDataValid( TYPE )                                           \
    template<>                                                             \
    Array<double> KeyDataScalar<TYPE>::convertToDouble() const             \
    {                                                                      \
        Array<double> data( 1 );                                           \
        data(0) = d_data;                                                  \
        return data;                                                       \
    }                                                                      \
    template<>                                                             \
    Array<int64_t> KeyDataScalar<TYPE>::convertToInt64() const             \
    {                                                                      \
        Array<int64_t> data( 1 );                                          \
        data(0) = d_data;                                                  \
        return data;                                                       \
    }                                                                      \
    template<>                                                             \
    Array<double> KeyDataArray<TYPE>::convertToDouble() const              \
    {                                                                      \
        Array<double> data( d_data.size() );                               \
        data.copy( d_data );                                               \
        return data;                                                       \
    }                                                                      \
    template<>                                                             \
    Array<int64_t> KeyDataArray<TYPE>::convertToInt64() const              \
    {                                                                      \
        Array<int64_t> data( d_data.size() );                              \
        data.copy( d_data );                                               \
        return data;                                                       \
    }                                                                      \
    template<>                                                             \
    Array<TYPE> convertFromDouble( const Array<double>& data )             \
    {                                                                      \
        Array<TYPE> data2( data.size() );                                  \
        data2.copy( data );                                                \
        return data2;                                                      \
    }
#define convertDataInvalid( TYPE )                                        \
    template<>                                                            \
    Array<double> KeyDataScalar<TYPE>::convertToDouble() const            \
    {                                                                     \
        throw std::logic_error( "Invalid conversion: " #TYPE "-double" ); \
    }                                                                     \
    template<>                                                            \
    Array<int64_t> KeyDataScalar<TYPE>::convertToInt64() const            \
    {                                                                     \
        throw std::logic_error( "Invalid conversion: " #TYPE"-int" );     \
    }                                                                     \
    template<>                                                            \
    Array<double> KeyDataArray<TYPE>::convertToDouble() const             \
    {                                                                     \
        throw std::logic_error( "Invalid conversion: " #TYPE "-double" ); \
    }                                                                     \
    template<>                                                            \
    Array<int64_t> KeyDataArray<TYPE>::convertToInt64() const             \
    {                                                                     \
        throw std::logic_error( "Invalid conversion: " #TYPE "-int" );    \
    }                                                                     \
    template<>                                                            \
    Array<TYPE> convertFromDouble( const Array<double>& )                 \
    {                                                                     \
        throw std::logic_error( "Invalid conversion: double-" #TYPE );    \
    }
convertDataValid( char )
convertDataValid( int8_t )
convertDataValid( int16_t )
convertDataValid( int32_t )
convertDataValid( int64_t )
convertDataValid( uint8_t )
convertDataValid( uint16_t )
convertDataValid( uint32_t )
convertDataValid( uint64_t )
convertDataValid( float )
convertDataValid( double )
convertDataValid( long double )
convertDataInvalid( std::complex<double> )
convertDataInvalid( std::complex<float> )
convertDataInvalid( bool )
convertDataInvalid( std::_Bit_reference )
convertDataInvalid( std::string )
convertDataInvalid( DatabaseBox )


/********************************************************************
 * KeyDataScalar::operator==                                         *
 ********************************************************************/
template<class TYPE>
static inline bool compare( const TYPE& x, const TYPE& y );
template<> inline bool compare( const double& x, const double& y )
{
    bool test = x == y;
    test = test || fabs( x - y ) <= 1e-12 * fabs( x + y );
    test = test || ( ( x != x ) && ( y != y ) );
    return test;
}
template<> inline bool compare( const float& x, const float& y )
{
    bool test = x == y;
    test = test || fabs( x - y ) <= 1e-7 * fabs( x + y );
    test = test || ( ( x != x ) && ( y != y ) );
    return test;
}
template<> inline bool compare( const std::complex<double>& x, const std::complex<double>& y )
{
    bool test = x == y;
    test = test || std::abs( x - y ) <= 1e-12 * std::abs( x + y );
    test = test || ( ( x != x ) && ( y != y ) );
    return test;
}
template<> inline bool compare( const std::complex<float>& x, const std::complex<float>& y )
{
    bool test = x == y;
    test = test || std::abs( x - y ) <= 1e-7 * std::abs( x + y );
    test = test || ( ( x != x ) && ( y != y ) );
    return test;
}
template<class TYPE>
static inline bool compare( const TYPE& x, const TYPE& y )
{
    return x == y;
}
#define compareKeyData( TYPE )                                        \
    template<>                                                        \
    bool KeyDataScalar<TYPE>::operator==( const KeyData& rhs ) const  \
    {                                                                 \
        auto tmp1 = dynamic_cast<const KeyDataScalar<TYPE>*>( &rhs ); \
        auto tmp2 = dynamic_cast<const KeyDataArray<TYPE>*>( &rhs );  \
        if ( tmp1 ) {                                                 \
            return compare( d_data, tmp1->d_data );                   \
        } else if ( tmp2 ) {                                          \
            if ( tmp2->get().size() != 1 )                            \
                return false;                                         \
            return compare( d_data, tmp2->get()(0) );                 \
        } else if ( ( is_floating_point() || is_integral() ) &&       \
                 ( rhs.is_floating_point() || rhs.is_integral() ) ) { \
            auto data1 = convertToDouble();                           \
            auto data2 = rhs.convertToDouble();                       \
            if ( data1.size() != data2.size() )                       \
                return false;                                         \
            bool test = true;                                         \
            for ( size_t i=0; i<data1.length(); i++)                  \
                test = test && compare( data1(i), data2(i) );         \
            return test;                                              \
        }                                                             \
        return false;                                                 \
    }                                                                 \
    template<>                                                        \
    bool KeyDataArray<TYPE>::operator==( const KeyData& rhs ) const   \
    {                                                                 \
        auto tmp1 = dynamic_cast<const KeyDataScalar<TYPE>*>( &rhs ); \
        auto tmp2 = dynamic_cast<const KeyDataArray<TYPE>*>( &rhs ); \
        if ( tmp1 ) {                                                 \
            if ( d_data.size() != 1 )                                 \
                return false;                                         \
            return compare( d_data(0), tmp1->get() );                 \
        } else if ( tmp2 ) {                                          \
            return compare( d_data, tmp2->get() );                    \
        } else if ( ( is_floating_point() || is_integral() ) &&       \
                 ( rhs.is_floating_point() || rhs.is_integral() ) ) { \
            auto data1 = convertToDouble();                           \
            auto data2 = rhs.convertToDouble();                       \
            if ( data1.size() != data2.size() )                       \
                return false;                                         \
            bool test = true;                                         \
            for ( size_t i=0; i<data1.length(); i++)                  \
                test = test && compare( data1(i), data2(i) );         \
            return test;                                              \
        }                                                             \
        return false;                                                 \
    }
compareKeyData( int8_t )
compareKeyData( int16_t )
compareKeyData( int32_t )
compareKeyData( int64_t )
compareKeyData( uint8_t )
compareKeyData( uint16_t )
compareKeyData( uint32_t )
compareKeyData( uint64_t )
compareKeyData( float )
compareKeyData( double )
compareKeyData( long double )
compareKeyData( std::complex<double> )
compareKeyData( std::complex<float> )
compareKeyData( bool )
compareKeyData( char )
compareKeyData( std::string )
compareKeyData( std::_Bit_reference )
compareKeyData( DatabaseBox )

    // clang-format on


    /********************************************************************
     * DatabaseBox                                                       *
     ********************************************************************/
    DatabaseBox::DatabaseBox()
    : d_dim( 0 )
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
template Array<DatabaseBox, FunctionTable> &Array<DatabaseBox, FunctionTable>::
operator=( const Array<DatabaseBox, FunctionTable> & );
template Array<DatabaseBox, FunctionTable> &Array<DatabaseBox, FunctionTable>::
operator=( Array<DatabaseBox, FunctionTable> && );
template Array<DatabaseBox, FunctionTable> &Array<DatabaseBox, FunctionTable>::
operator=( const std::vector<DatabaseBox> & );
template void Array<DatabaseBox, FunctionTable>::clear();
template void
Array<DatabaseBox, FunctionTable>::viewRaw( ArraySize const &, DatabaseBox *, bool, bool );
template bool Array<DatabaseBox, FunctionTable>::
operator==( Array<DatabaseBox, FunctionTable> const & ) const;
template void Array<DatabaseBox, FunctionTable>::resize( ArraySize const & );


} // namespace AMP
