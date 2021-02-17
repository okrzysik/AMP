#ifndef included_extended_int_hpp
#define included_extended_int_hpp

#include <algorithm>
#include <limits>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <string>
#include <vector>


namespace AMP::extended {


/********************************************************************
 * Constructors                                                      *
 ********************************************************************/
template<unsigned char N>
constexpr int64N<N>::int64N() : data( { 0 } )
{
}
template<unsigned char N>
template<unsigned char N2>
constexpr int64N<N>::int64N( const int64N<N2> &rhs ) : data( { 0 } )
{
    uint64_t fill = ( rhs.data.u64[N2 - 1] >> 63 ) ? ( ~( (uint64_t) 0 ) ) : 0;
    for ( size_t i = 0; i < N; i++ )
        data.u64[i] = fill;
    for ( size_t i = 0; i < std::min( N, N2 ); i++ )
        data.u64[i] = rhs.data.u64[i];
}
template<unsigned char N>
constexpr int64N<N>::int64N( const int64_t &rhs ) : data( { 0 } )
{
    uint64_t fill = ( rhs < 0 ) ? ( ~( (uint64_t) 0 ) ) : 0;
    for ( size_t i = 0; i < N; i++ )
        data.u64[i] = fill;
    data.u64[0] = static_cast<uint64_t>( rhs );
}
template<unsigned char N>
constexpr int64N<N>::int64N( const int &rhs ) : data( { 0 } )
{
    uint64_t fill = ( rhs < 0 ) ? ( ~( (uint64_t) 0 ) ) : 0;
    for ( size_t i = 0; i < N; i++ )
        data.u64[i] = fill;
    data.u64[0] = static_cast<uint64_t>( static_cast<int64_t>( rhs ) );
}
template<unsigned char N>
constexpr int64N<N>::int64N( const char *rhs ) : data( { 0 } )
{
    for ( size_t i = 0; i < N; i++ )
        data.u64[i] = 0;
    int length = 0;
    while ( rhs[length] != 0 )
        length++;
    if ( rhs[0] == '0' && rhs[1] == 'b' ) {
        // Read base 2 string
        throw std::logic_error( "Not finished" );
    } else if ( rhs[0] == '0' && rhs[1] == 'x' ) {
        // Reading a base 16 string
        int Nb = std::min<int>( ( length + 7 - 2 ) / 8, 2 * N );
        uint64_t blocks[2 * N];
        memset( blocks, 0, sizeof( blocks ) );
        for ( int i = 0; i < Nb; i++ ) {
            char tmp[9];
            memset( tmp, 0, sizeof( tmp ) );
            int N2 = std::min<int>( length - 2 - 8 * i, 8 );
            memcpy( tmp, &rhs[length - i * 8 - N2], N2 );
            blocks[i] = static_cast<uint64_t>( strtoul( tmp, NULL, 16 ) );
        }
        for ( int i = 0; i < N; i++ )
            data.u64[i] = blocks[2 * i] + ( blocks[2 * i + 1] << 32 );
    } else {
        // Read base 10 string
        /*for (int i=0; i<length; i++) {
            if ( rhs[i]<48 || rhs[i]>59 ) {
                if ( rhs[i]!='+' && rhs[i]!='-' && rhs[i]!='e' )
                    throw std::logic_error("Invalid initialization");
            }
        }
        int i=0;
        int sign = 1;
        if ( rhs[0]=='-' || rhs[0]=='+' ) {
            i = 1;
            sign = (rhs[0]=='-') ? -1:1;
        }*/
        throw std::logic_error( "Not finished" );
    }
}


/********************************************************************
 * Some simple functions                                             *
 ********************************************************************/
template<unsigned char N>
constexpr int64_t int64N<N>::get_int64() const
{
    return static_cast<int64_t>( data.u64[0] );
}
template<unsigned char N>
constexpr int int64N<N>::get_int() const
{
    return static_cast<int>( data.u64[0] );
}
template<unsigned char N>
constexpr double int64N<N>::get_double() const
{
    uint64_t data2[N];
    for ( size_t i = 0; i < N; i++ )
        data2[i] = data.u64[i];
    int s = sign();
    if ( s == -1 ) {
        for ( size_t i = 0; i < N; i++ )
            data2[i] = ~data2[i];
        data2[0]++;
        int i = 1;
        while ( data2[i - 1] == 0 && i < N ) {
            data2[i]++;
            i++;
        }
    }
    const double scale = 1.84467440737095516e19; // 2^64
    double result      = 0.0;
    double tmp         = 1.0;
    for ( size_t i = 0; i < N; i++ ) {
        if ( data2[i] != 0 )
            result += tmp * static_cast<double>( data2[i] );
        tmp *= scale;
    }
    return s * result;
}
template<unsigned char N>
std::string int64N<N>::get_hex() const
{
    char hex[16 * N + 3] = "0x";
    for ( size_t i = 0; i < 2 * N; i++ ) {
        char tmp[64];
        sprintf( tmp, "%08x", static_cast<unsigned int>( data.u32[i] ) );
        memcpy( &hex[8 * ( 2 * N - i - 1 ) + 2], tmp, 8 );
    }
    hex[16 * N + 2] = 0;
    return hex;
}
template<unsigned char N>
constexpr int int64N<N>::sign() const
{
    constexpr uint64_t mask = 0x8000000000000000;
    return ( data.u64[N - 1] & mask ) == 0 ? 1 : -1;
}
template<unsigned char N>
constexpr int64N<N> int64N<N>::operator!() const
{
    int64N<N> y( *this );
    for ( size_t i = 0; i < N; i++ )
        y.data.u64[i] = ~y.data.u64[i];
    return y;
}
template<unsigned char N>
constexpr void int64N<N>::compliment()
{
    for ( size_t i = 0; i < N; i++ )
        data.u64[i] = ~data.u64[i];
    this->operator+=( static_cast<uint64_t>( 1u ) );
}


/********************************************************************
 * Arithmetic functions                                              *
 * Note: this and x may be the same object                           *
 ********************************************************************/
template<unsigned char N>
constexpr int64N<N> &int64N<N>::operator+=( const int64N<N> &x )
{
    bool carry = false;
    for ( size_t i = 0; i < N; i++ ) {
        uint64_t data0 = data.u64[i];
        data.u64[i] += x.data.u64[i] + carry;
        carry = ( ( data.u64[i] < data0 ) || ( data.u64[i] < x.data.u64[i] ) );
    }
    return *this;
}
template<unsigned char N>
constexpr int64N<N> &int64N<N>::operator+=( const int64_t x )
{
    return operator+=( static_cast<uint64_t>( x ) );
}
template<unsigned char N>
constexpr int64N<N> &int64N<N>::operator+=( const uint64_t x )
{
    data.u64[0] = data.u64[0] + x;
    bool carry  = data.u64[0] < x;
    size_t i    = 1;
    while ( carry && i < N ) {
        ++data.u64[i];
        carry = data.u64[i] == 0;
        ++i;
    }
    return *this;
}
template<unsigned char N>
constexpr int64N<N> &int64N<N>::operator-=( const int64N<N> &x )
{
    uint64_t carry = 1;
    for ( size_t i = 0; i < N; i++ ) {
        uint64_t data0 = data.u64[i];
        uint64_t data1 = ~( x.data.u64[i] );
        data.u64[i] += data1 + carry;
        carry = ( ( data.u64[i] < data0 ) || ( data.u64[i] < data1 ) );
    }
    return *this;
}
template<unsigned char N>
constexpr int64N<N> &int64N<N>::operator*=( const int64N<N> &x )
{
    // Break the numbers into blocks of 32 bits
    // Note:  only works for little-endian
    uint64_t a[2 * N], b[2 * N];
    for ( size_t i = 0; i < 2 * N; i++ ) {
        a[i] = data.u32[i];
        b[i] = x.data.u32[i];
    }
    // Multiply the blocks
    uint64_t c[2 * N + 1], carry0[2 * N + 2], carry[N + 1];
    for ( size_t i = 0; i < 2 * N + 2; i++ )
        carry0[i] = 0;
    for ( size_t j = 0; j < 2 * N; j++ )
        c[j] = a[0] * b[j];
    for ( size_t i = 1; i < 2 * N; i++ ) {
        for ( size_t j = 0, j2 = i; j < 2 * N - i; ++j, ++j2 ) {
            uint64_t tmp = c[j2];
            c[i + j] += a[i] * b[j];
            carry0[i + j + 2] += c[i + j] < tmp;
        }
    }
    for ( size_t i = 0; i < N; i++ )
        carry[i] = carry0[2 * i] + ( carry0[2 * i + 1] << 32 );
    // Combine the blocks
    for ( size_t i = 0; i < N; i++ ) {
        data.u64[i] = c[2 * i] + ( c[2 * i + 1] << 32 ) + carry[i];
        carry[i + 1] += c[2 * i + 1] >> 32;
        carry[i + 1] += data.u64[i] < c[2 * i];
    }
    return *this;
}
template<unsigned char N>
constexpr int64N<N> &int64N<N>::operator/=( const int64N<N> & )
{
    throw std::logic_error( "Not finished" );
    return *this;
}


/********************************************************************
 * Shift operator overloading                                        *
 ********************************************************************/
template<unsigned char N>
constexpr int64N<N> &int64N<N>::operator>>=( int shift )
{
    return this->operator<<=( -shift );
}
template<unsigned char N>
constexpr int64N<N> &int64N<N>::operator<<=( int shift )
{
    if ( shift > 0 ) {
        // Left shift
        while ( shift >= 64 ) {
            for ( int i = N - 1; i > 0; i-- )
                data.u64[i] = data.u64[i - 1];
            data.u64[0] = 0;
            shift -= 64;
        }
        for ( int i = N - 1; i > 0; i-- )
            data.u64[i] = ( data.u64[i] << shift ) + ( data.u64[i - 1] >> ( 64 - shift ) );
        data.u64[0] <<= shift;
    } else {
        // Right shift
        while ( -shift >= 64 ) {
            for ( int i = 0; i < N - 1; i++ )
                data.u64[i] = data.u64[i + 1];
            data.u64[N - 1] = 0;
            shift += 64;
        }
        for ( int i = 0; i < N - 1; i++ )
            data.u64[i] = ( data.u64[i] >> -shift ) + ( data.u64[i + 1] << ( 64 + shift ) );
        data.u64[N - 1] >>= -shift;
    }
    return *this;
}


/********************************************************************
 * Comparison operator overloading                                   *
 ********************************************************************/
template<unsigned char N>
constexpr bool int64N<N>::operator==( const int64N<N> &rhs ) const
{
    bool equal = true;
    for ( int i = 0; i < N; i++ )
        equal &= data.u64[i] == rhs.data.u64[i];
    return equal;
}
template<unsigned char N>
constexpr bool int64N<N>::operator!=( const int64N<N> &rhs ) const
{
    return !operator==( rhs );
}
template<unsigned char N>
constexpr bool int64N<N>::operator>( const int64N<N> &rhs ) const
{
    if ( this->sign() != rhs.sign() )
        return this->sign() > rhs.sign();
    if ( this->sign() > rhs.sign() )
        for ( int i = N - 1; i > 0; i-- ) {
            if ( data.u64[i] < rhs.data.u64[i] )
                return false;
            if ( data.u64[i] > rhs.data.u64[i] )
                return true;
        }
    return data.u64[0] > rhs.data.u64[0];
}
template<unsigned char N>
constexpr bool int64N<N>::operator<( const int64N<N> &rhs ) const
{
    if ( this->sign() != rhs.sign() )
        return this->sign() < rhs.sign();
    for ( int i = N - 1; i > 0; i-- ) {
        if ( data.u64[i] > rhs.data.u64[i] )
            return false;
        if ( data.u64[i] < rhs.data.u64[i] )
            return true;
    }
    return data.u64[0] < rhs.data.u64[0];
}
template<unsigned char N>
constexpr bool int64N<N>::operator>=( const int64N<N> &rhs ) const
{
    if ( this->sign() != rhs.sign() )
        return this->sign() > rhs.sign();
    for ( int i = N - 1; i > 0; i-- ) {
        if ( data.u64[i] < rhs.data.u64[i] )
            return false;
        if ( data.u64[i] > rhs.data.u64[i] )
            return true;
    }
    return data.u64[0] >= rhs.data.u64[0];
}
template<unsigned char N>
constexpr bool int64N<N>::operator<=( const int64N<N> &rhs ) const
{
    if ( this->sign() != rhs.sign() )
        return this->sign() < rhs.sign();
    for ( int i = N - 1; i > 0; i-- ) {
        if ( data.u64[i] > rhs.data.u64[i] )
            return false;
        if ( data.u64[i] < rhs.data.u64[i] )
            return true;
    }
    return data.u64[0] <= rhs.data.u64[0];
}
template<unsigned char N>
constexpr bool int64N<N>::operator==( int rhs ) const
{
    return operator==( int64N<N>( rhs ) );
}
template<unsigned char N>
constexpr bool int64N<N>::operator!=( int rhs ) const
{
    return operator!=( int64N<N>( rhs ) );
}
template<unsigned char N>
constexpr bool int64N<N>::operator>( int rhs ) const
{
    return operator>( int64N<N>( rhs ) );
}
template<unsigned char N>
constexpr bool int64N<N>::operator<( int rhs ) const
{
    return operator<( int64N<N>( rhs ) );
}
template<unsigned char N>
constexpr bool int64N<N>::operator>=( int rhs ) const
{
    return operator>=( int64N<N>( rhs ) );
}
template<unsigned char N>
constexpr bool int64N<N>::operator<=( int rhs ) const
{
    return operator<=( int64N<N>( rhs ) );
}


/********************************************************************
 * Arithmetic operator overloading                                   *
 ********************************************************************/
template<unsigned char N>
constexpr int64N<N> operator+( const int64N<N> &x, const int64N<N> &y )
{
    int64N<N> result( x );
    result += y;
    return result;
}
template<unsigned char N>
constexpr int64N<N> operator-( const int64N<N> &x, const int64N<N> &y )
{
    int64N<N> result( x );
    result -= y;
    return result;
}
template<unsigned char N>
constexpr int64N<N> operator-( const int64_t &x, const int64N<N> &y )
{
    int64N<N> result( x );
    result -= y;
    return result;
}
template<unsigned char N>
constexpr int64N<N> operator-( const int64N<N> &x )
{
    int64N<N> result( x );
    result.compliment();
    return result;
}
template<unsigned char N>
constexpr int64N<N> operator*( const int64N<N> &x, const int64N<N> &y )
{
    int64N<N> result( x );
    result *= y;
    return result;
}
template<unsigned char N>
constexpr int64N<N> operator/( const int64N<N> &x, const int64N<N> &y )
{
    int64N<N> result( x );
    result /= y;
    return result;
}


/********************************************************************
 * ostream and shift operator overloading                            *
 ********************************************************************/
template<unsigned char N>
std::ostream &operator<<( std::ostream &out, const int64N<N> &x )
{
    out << x.get_hex();
    return out;
}
template<unsigned char N>
constexpr int64N<N> operator<<( const int64N<N> &rhs, int shift )
{
    int64N<N> result( rhs );
    result <<= shift;
    return result;
}
template<unsigned char N>
constexpr int64N<N> operator>>( const int64N<N> &rhs, int shift )
{
    int64N<N> result( rhs );
    result >>= shift;
    return result;
}


} // namespace AMP::extended


#endif
