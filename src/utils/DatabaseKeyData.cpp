#include "AMP/utils/Database.h"

#include <complex>
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
#define scaleDataValid( TYPE )                               \
    template<>                                               \
    void scaleData( std::vector<TYPE>& data, double factor ) \
    {                                                        \
        for ( size_t i = 0; i < data.size(); i++ )           \
            data[i] = static_cast<TYPE>( factor * data[i] ); \
    }                                                        \
    template<>                                               \
    void scaleData( TYPE& data, double factor )              \
    {                                                        \
        data = static_cast<TYPE>( factor * data );           \
    }
#define scaleDataInvalid( TYPE )                          \
    template<>                                            \
    void scaleData( std::vector<TYPE>&, double )          \
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
scaleDataValid( std::complex<double> )
scaleDataInvalid( bool )
scaleDataInvalid( char )
scaleDataInvalid( std::string )
template<>
void scaleData( std::vector<std::complex<float>>& data, double factor )
{
    for ( size_t i = 0; i < data.size(); i++ )
        data[i] = static_cast<float>( factor ) * data[i];
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
    std::vector<double> KeyDataScalar<TYPE>::convertToDouble() const       \
    {                                                                      \
        return std::vector<double>( 1, static_cast<double>( d_data ) );    \
    }                                                                      \
    template<>                                                             \
    std::vector<int64_t> KeyDataScalar<TYPE>::convertToInt64() const       \
    {                                                                      \
        return std::vector<int64_t>( 1, static_cast<int64_t>( d_data ) );  \
    }                                                                      \
    template<>                                                             \
    std::vector<double> KeyDataVector<TYPE>::convertToDouble() const       \
    {                                                                      \
        std::vector<double> data( d_data.size() );                         \
        for ( size_t i = 0; i < data.size(); i++ )                         \
            data[i] = static_cast<double>( d_data[i] );                    \
        return data;                                                       \
    }                                                                      \
    template<>                                                             \
    std::vector<int64_t> KeyDataVector<TYPE>::convertToInt64() const       \
    {                                                                      \
        std::vector<int64_t> data( d_data.size() );                        \
        for ( size_t i = 0; i < data.size(); i++ )                         \
            data[i] = static_cast<int64_t>( d_data[i] );                   \
        return data;                                                       \
    }                                                                      \
    template<>                                                             \
    std::vector<TYPE> convertFromDouble( const std::vector<double>& data ) \
    {                                                                      \
        std::vector<TYPE> data2( data.size() );                            \
        for ( size_t i = 0; i < data.size(); i++ )                         \
            data2[i] = static_cast<TYPE>( data[i] );                       \
        return data2;                                                      \
    }
#define convertDataInvalid( TYPE )                                    \
    template<>                                                        \
    std::vector<double> KeyDataScalar<TYPE>::convertToDouble() const  \
    {                                                                 \
        throw std::logic_error( "Invalid conversion" );               \
    }                                                                 \
    template<>                                                        \
    std::vector<int64_t> KeyDataScalar<TYPE>::convertToInt64() const  \
    {                                                                 \
        throw std::logic_error( "Invalid conversion" );               \
    }                                                                 \
    template<>                                                        \
    std::vector<double> KeyDataVector<TYPE>::convertToDouble() const  \
    {                                                                 \
        throw std::logic_error( "Invalid conversion" );               \
    }                                                                 \
    template<>                                                        \
    std::vector<int64_t> KeyDataVector<TYPE>::convertToInt64() const  \
    {                                                                 \
        throw std::logic_error( "Invalid conversion" );               \
    }                                                                 \
    template<>                                                        \
    std::vector<TYPE> convertFromDouble( const std::vector<double>& ) \
    {                                                                 \
        throw std::logic_error( "Invalid conversion" );               \
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
convertDataInvalid( std::complex<double> )
convertDataInvalid( std::complex<float> )
convertDataInvalid( bool )
convertDataInvalid( std::_Bit_reference )
convertDataInvalid( std::string )


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
        auto tmp2 = dynamic_cast<const KeyDataVector<TYPE>*>( &rhs ); \
        if ( tmp1 ) {                                                 \
            return compare( d_data, tmp1->d_data );                   \
        } else if ( tmp2 ) {                                          \
            if ( tmp2->get().size() != 1 )                            \
                return false;                                         \
            return compare( d_data, tmp2->get()[0] );                 \
        } else if ( ( is_floating_point() || is_integral() ) &&       \
                 ( rhs.is_floating_point() || rhs.is_integral() ) ) { \
            auto data1 = convertToDouble();                           \
            auto data2 = rhs.convertToDouble();                       \
            if ( data1.size() != data2.size() )                       \
                return false;                                         \
            bool test = true;                                         \
            for ( size_t i=0; i<data1.size(); i++)                    \
                test = test && compare( data1[i], data2[i] );         \
            return test;                                              \
        }                                                             \
        return false;                                                 \
    }                                                                 \
    template<>                                                        \
    bool KeyDataVector<TYPE>::operator==( const KeyData& rhs ) const  \
    {                                                                 \
        auto tmp1 = dynamic_cast<const KeyDataScalar<TYPE>*>( &rhs ); \
        auto tmp2 = dynamic_cast<const KeyDataVector<TYPE>*>( &rhs ); \
        if ( tmp1 ) {                                                 \
            if ( d_data.size() != 1 )                                 \
                return false;                                         \
            return compare( d_data[0], tmp1->get() );                 \
        } else if ( tmp2 ) {                                          \
            return compare( d_data, tmp2->get() );                    \
        } else if ( ( is_floating_point() || is_integral() ) &&       \
                 ( rhs.is_floating_point() || rhs.is_integral() ) ) { \
            auto data1 = convertToDouble();                           \
            auto data2 = rhs.convertToDouble();                       \
            if ( data1.size() != data2.size() )                       \
                return false;                                         \
            bool test = true;                                         \
            for ( size_t i=0; i<data1.size(); i++)                    \
                test = test && compare( data1[i], data2[i] );         \
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
compareKeyData( std::complex<double> )
compareKeyData( std::complex<float> )
compareKeyData( bool )
compareKeyData( char )
compareKeyData( std::string )
compareKeyData( std::_Bit_reference )

}

// clang-format on
