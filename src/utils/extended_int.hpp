#ifndef included_extended_int_hpp
#define included_extended_int_hpp

#include <algorithm>
#include <limits>
#include <stdexcept>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <string>
#include <vector>


namespace AMP::extended {


/********************************************************************
 * Constructors                                                      *
 ********************************************************************/
template<uint8_t N>
constexpr int64N<N>::int64N() : data{ 0 }
{
}
template<uint8_t N>
template<uint8_t N2>
constexpr int64N<N>::int64N( const int64N<N2> &rhs ) : data{ 0 }
{
    uint64_t fill = ( rhs.data[N2 - 1] >> 63 ) ? ( ~( (uint64_t) 0 ) ) : 0;
    for ( size_t i = 0; i < N; i++ )
        data[i] = fill;
    for ( size_t i = 0; i < std::min( N, N2 ); i++ )
        data[i] = rhs.data[i];
}
template<uint8_t N>
constexpr int64N<N>::int64N( int64_t rhs ) : data{ 0 }
{
    uint64_t fill = ( rhs < 0 ) ? ( ~( (uint64_t) 0 ) ) : 0;
    for ( size_t i = 0; i < N; i++ )
        data[i] = fill;
    data[0] = static_cast<uint64_t>( rhs );
}
template<uint8_t N>
constexpr int64N<N>::int64N( int rhs ) : data{ 0 }
{
    uint64_t fill = ( rhs < 0 ) ? ( ~( (uint64_t) 0 ) ) : 0;
    for ( size_t i = 0; i < N; i++ )
        data[i] = fill;
    data[0] = static_cast<uint64_t>( static_cast<int64_t>( rhs ) );
}
template<uint8_t N>
constexpr int64N<N>::int64N( const char *str ) : data{ 0 }
{
    for ( size_t i = 0; i < N; i++ )
        data[i] = 0;
    int length = 0;
    while ( str[length] != 0 )
        length++;
    if ( str[0] == '0' && str[1] == 'b' ) {
        // Read base 2 string
        throw std::logic_error( "Not finished" );
    } else if ( str[0] == '0' && str[1] == 'x' ) {
        // Reading a base 16 string
        if ( length > 16 * N + 3 )
            throw "hex does not fit in integer";
        for ( int i = 2; i < length; i++ ) {
            uint8_t v = 0;
            if ( str[i] >= 48 && str[i] <= 57 )
                v = str[i] - 48;
            else if ( str[i] >= 65 && str[i] <= 70 )
                v = str[i] - 65 + 10;
            else if ( str[i] >= 97 && str[i] <= 102 )
                v = str[i] - 97 + 10;
            else
                throw std::logic_error( "Invalid char in hex string" );
            operator<<=( 4 );
            data[0] += v;
        }
    } else {
        // Read base 10 string
        /*for (int i=0; i<length; i++) {
            if ( str[i]<48 || str[i]>59 ) {
                if ( str[i]!='+' && str[i]!='-' && str[i]!='e' )
                    throw std::logic_error("Invalid initialization");
            }
        }
        int i=0;
        int sign = 1;
        if ( str[0]=='-' || str[0]=='+' ) {
            i = 1;
            sign = (str[0]=='-') ? -1:1;
        }*/
        throw std::logic_error( "Not finished" );
    }
}


/********************************************************************
 * Convert to a hex string                                           *
 ********************************************************************/
template<uint8_t N>
constexpr std::array<char, 16 * N + 3> int64N<N>::hex( bool fixedWidth ) const
{
    constexpr char hexmap[]         = "0123456789abcdef";
    std::array<char, 16 *N + 3> hex = { 0 };
    hex[0]                          = '0';
    hex[1]                          = 'x';
    size_t k                        = 2;
    for ( int i = N - 1; i >= 0; i-- ) {
        for ( int j = 60; j >= 0; j -= 4 ) {
            char c = hexmap[( data[i] >> j ) & 0x0F];
            if ( k > 2 || c != '0' || fixedWidth )
                hex[k++] = c;
        }
    }
    return hex;
}


/********************************************************************
 * Some simple functions                                             *
 ********************************************************************/
template<uint8_t N>
constexpr int64N<N>::operator int64_t() const
{
    return static_cast<int64_t>( data[0] );
}
template<uint8_t N>
constexpr int64N<N>::operator int() const
{
    return static_cast<int>( static_cast<int64_t>( data[0] ) );
}
template<uint8_t N>
constexpr int64N<N>::operator double() const
{
    // Split the data into sets of unsigned 64-bit numbers
    uint64_t data2[N] = { 0 };
    for ( size_t i = 0; i < N; i++ )
        data2[i] = data[i];
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
    // Check for overflow
    uint8_t N2 = N;
    if constexpr ( N > 16 ) {
        N2            = 16;
        bool overflow = false;
        for ( size_t i = 16; i < N; i++ )
            overflow = overflow || data2[i] != 0;
        if ( overflow ) {
            if ( s >= 0 )
                return std::numeric_limits<double>::infinity();
            else
                return -std::numeric_limits<double>::infinity();
        }
    }
    // Convert to double
    // Note: using long double in not support in constexpr on IBM
    constexpr double scale[16] = { 1.0,
                                   1.8446744073709551616e19,
                                   3.4028236692093846346e38,
                                   6.2771017353866807638e57,
                                   1.1579208923731619542e77,
                                   2.1359870359209100823e96,
                                   3.9402006196394479212e115,
                                   7.2683872429560689054e134,
                                   1.3407807929942597099e154,
                                   2.4733040147310453406e173,
                                   4.5624406176221952186e192,
                                   8.4162174424773976115e211,
                                   1.5525180923007089351e231,
                                   2.8638903918474961204e250,
                                   5.2829453113566524635e269,
                                   9.7453140113999990803e288 };
    double result              = 0.0;
    for ( size_t i = 0; i < N2; i++ ) {
        if ( data2[i] != 0 )
            result += scale[i] * static_cast<double>( data2[i] );
    }
    return s * result;
}
template<uint8_t N>
constexpr int int64N<N>::sign() const
{
    constexpr uint64_t mask = 0x8000000000000000;
    return ( data[N - 1] & mask ) == 0 ? 1 : -1;
}
template<uint8_t N>
constexpr int64N<N> int64N<N>::operator!() const
{
    int64N<N> y( *this );
    for ( size_t i = 0; i < N; i++ )
        y.data[i] = ~y.data[i];
    return y;
}
template<uint8_t N>
constexpr void int64N<N>::compliment()
{
    for ( size_t i = 0; i < N; i++ )
        data[i] = ~data[i];
    this->operator+=( static_cast<uint64_t>( 1u ) );
}


/********************************************************************
 * Arithmetic functions                                              *
 * Note: this and x may be the same object                           *
 ********************************************************************/
template<uint8_t N>
constexpr int64N<N> &int64N<N>::operator+=( const int64N<N> &x )
{
    bool carry = false;
    for ( size_t i = 0; i < N; i++ ) {
        uint64_t data0 = data[i];
        data[i] += x.data[i] + ( carry ? 1 : 0 );
        carry = ( ( data[i] < data0 ) || ( data[i] < x.data[i] ) );
    }
    return *this;
}
template<uint8_t N>
constexpr int64N<N> &int64N<N>::operator+=( const int64_t x )
{
    return operator+=( static_cast<uint64_t>( x ) );
}
template<uint8_t N>
constexpr int64N<N> &int64N<N>::operator+=( const uint64_t x )
{
    data[0]    = data[0] + x;
    bool carry = data[0] < x;
    size_t i   = 1;
    while ( carry && i < N ) {
        ++data[i];
        carry = data[i] == 0;
        ++i;
    }
    return *this;
}
template<uint8_t N>
constexpr int64N<N> &int64N<N>::operator-=( const int64N<N> &x )
{
    uint64_t carry = 1;
    for ( size_t i = 0; i < N; i++ ) {
        uint64_t data0 = data[i];
        uint64_t data1 = ~( x.data[i] );
        data[i] += data1 + carry;
        carry = ( ( data[i] < data0 ) || ( data[i] < data1 ) );
    }
    return *this;
}
template<uint8_t N>
constexpr int64N<N> &int64N<N>::operator*=( const int64N<N> &x )
{
    // Break the numbers into blocks of 32 bits
    // Note:  only works for little-endian
    uint64_t a[2 * N] = { 0 }, b[2 * N] = { 0 };
    for ( size_t i = 0; i < N; i++ ) {
        a[2 * i + 0] = data[i] & 0xFFFFFFFF;
        b[2 * i + 0] = x.data[i] & 0xFFFFFFFF;
        a[2 * i + 1] = data[i] >> 32;
        b[2 * i + 1] = x.data[i] >> 32;
    }
    // Multiply the blocks
    uint64_t c[2 * N + 1] = { 0 }, carry0[2 * N + 2] = { 0 }, carry[N + 1] = { 0 };
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
        data[i] = c[2 * i] + ( c[2 * i + 1] << 32 ) + carry[i];
        carry[i + 1] += c[2 * i + 1] >> 32;
        carry[i + 1] += data[i] < c[2 * i];
    }
    return *this;
}
template<uint8_t N>
constexpr int64N<N> &int64N<N>::operator/=( const int64N<N> & )
{
    throw std::logic_error( "Not finished" );
    return *this;
}


/********************************************************************
 * Shift operator overloading                                        *
 ********************************************************************/
template<uint8_t N>
constexpr int64N<N> &int64N<N>::operator>>=( unsigned shift )
{
    // Right shift
    while ( shift >= 64 ) {
        for ( int i = 0; i < N - 1; i++ )
            data[i] = data[i + 1];
        data[N - 1] = 0;
        shift -= 64;
    }
    for ( int i = 0; i < N - 1; i++ )
        data[i] = ( data[i] >> shift ) + ( data[i + 1] << ( 64 - shift ) );
    data[N - 1] >>= shift;
    return *this;
}
template<uint8_t N>
constexpr int64N<N> &int64N<N>::operator<<=( unsigned shift )
{
    while ( shift >= 64 ) {
        for ( int i = N - 1; i > 0; i-- )
            data[i] = data[i - 1];
        data[0] = 0;
        shift -= 64;
    }
    for ( int i = N - 1; i > 0; i-- )
        data[i] = ( data[i] << shift ) + ( data[i - 1] >> ( 64 - shift ) );
    data[0] <<= shift;
    return *this;
}


/********************************************************************
 * Comparison operator overloading                                   *
 ********************************************************************/
template<uint8_t N>
constexpr bool int64N<N>::operator==( const int64N<N> &rhs ) const
{
    bool equal = true;
    for ( int i = 0; i < N; i++ )
        equal &= data[i] == rhs.data[i];
    return equal;
}
template<uint8_t N>
constexpr bool int64N<N>::operator!=( const int64N<N> &rhs ) const
{
    return !operator==( rhs );
}
template<uint8_t N>
constexpr bool int64N<N>::operator>( const int64N<N> &rhs ) const
{
    if ( this->sign() != rhs.sign() )
        return this->sign() > rhs.sign();
    if ( this->sign() > rhs.sign() )
        for ( int i = N - 1; i > 0; i-- ) {
            if ( data[i] < rhs.data[i] )
                return false;
            if ( data[i] > rhs.data[i] )
                return true;
        }
    return data[0] > rhs.data[0];
}
template<uint8_t N>
constexpr bool int64N<N>::operator<( const int64N<N> &rhs ) const
{
    if ( this->sign() != rhs.sign() )
        return this->sign() < rhs.sign();
    for ( int i = N - 1; i > 0; i-- ) {
        if ( data[i] > rhs.data[i] )
            return false;
        if ( data[i] < rhs.data[i] )
            return true;
    }
    return data[0] < rhs.data[0];
}
template<uint8_t N>
constexpr bool int64N<N>::operator>=( const int64N<N> &rhs ) const
{
    if ( this->sign() != rhs.sign() )
        return this->sign() > rhs.sign();
    for ( int i = N - 1; i > 0; i-- ) {
        if ( data[i] < rhs.data[i] )
            return false;
        if ( data[i] > rhs.data[i] )
            return true;
    }
    return data[0] >= rhs.data[0];
}
template<uint8_t N>
constexpr bool int64N<N>::operator<=( const int64N<N> &rhs ) const
{
    if ( this->sign() != rhs.sign() )
        return this->sign() < rhs.sign();
    for ( int i = N - 1; i > 0; i-- ) {
        if ( data[i] > rhs.data[i] )
            return false;
        if ( data[i] < rhs.data[i] )
            return true;
    }
    return data[0] <= rhs.data[0];
}
template<uint8_t N>
constexpr bool int64N<N>::operator==( int rhs ) const
{
    return operator==( int64N<N>( rhs ) );
}
template<uint8_t N>
constexpr bool int64N<N>::operator!=( int rhs ) const
{
    return operator!=( int64N<N>( rhs ) );
}
template<uint8_t N>
constexpr bool int64N<N>::operator>( int rhs ) const
{
    return operator>( int64N<N>( rhs ) );
}
template<uint8_t N>
constexpr bool int64N<N>::operator<( int rhs ) const
{
    return operator<( int64N<N>( rhs ) );
}
template<uint8_t N>
constexpr bool int64N<N>::operator>=( int rhs ) const
{
    return operator>=( int64N<N>( rhs ) );
}
template<uint8_t N>
constexpr bool int64N<N>::operator<=( int rhs ) const
{
    return operator<=( int64N<N>( rhs ) );
}


/********************************************************************
 * Arithmetic operator overloading                                   *
 ********************************************************************/
template<uint8_t N>
constexpr int64N<N> operator+( const int64N<N> &x, const int64N<N> &y )
{
    int64N<N> result( x );
    result += y;
    return result;
}
template<uint8_t N>
constexpr int64N<N> operator-( const int64N<N> &x, const int64N<N> &y )
{
    int64N<N> result( x );
    result -= y;
    return result;
}
template<uint8_t N>
constexpr int64N<N> operator-( const int64_t &x, const int64N<N> &y )
{
    int64N<N> result( x );
    result -= y;
    return result;
}
template<uint8_t N>
constexpr int64N<N> operator-( const int64N<N> &x )
{
    int64N<N> result( x );
    result.compliment();
    return result;
}
template<uint8_t N>
constexpr int64N<N> operator*( const int64N<N> &x, const int64N<N> &y )
{
    int64N<N> result( x );
    result *= y;
    return result;
}
template<uint8_t N>
constexpr int64N<N> operator/( const int64N<N> &x, const int64N<N> &y )
{
    int64N<N> result( x );
    result /= y;
    return result;
}


/********************************************************************
 * ostream and shift operator overloading                            *
 ********************************************************************/
template<uint8_t N>
std::ostream &operator<<( std::ostream &out, const int64N<N> &x )
{
    out << x.hex().data();
    return out;
}
template<uint8_t N>
constexpr int64N<N> operator<<( const int64N<N> &rhs, int shift )
{
    int64N<N> result( rhs );
    result <<= shift;
    return result;
}
template<uint8_t N>
constexpr int64N<N> operator>>( const int64N<N> &rhs, int shift )
{
    int64N<N> result( rhs );
    result >>= shift;
    return result;
}


} // namespace AMP::extended


#endif
