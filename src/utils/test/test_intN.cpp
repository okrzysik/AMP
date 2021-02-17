#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <ctime>
#include <iostream>
#include <random>
#include <stdexcept>
#include <string>
#include <sys/stat.h>
#include <vector>

#include "AMP/utils/Utilities.h"
#include "AMP/utils/extended_int.h"

inline double time() { return AMP::Utilities::time(); }

#define printp printf

using eint64   = AMP::extended::int64N<1>;
using eint128  = AMP::extended::int128_t;
using eint256  = AMP::extended::int256_t;
using eint512  = AMP::extended::int512_t;
using eint1024 = AMP::extended::int1024_t;
using eint2048 = AMP::extended::int2048_t;


// Check the sizes of the different types
static_assert( sizeof( eint64 ) == 8, " eint64 is not 8 bytes" );
static_assert( sizeof( eint128 ) == 16, " eint128 is not 16 bytes" );
static_assert( sizeof( eint256 ) == 32, " eint256 is not 32 bytes" );
static_assert( sizeof( eint512 ) == 64, " eint512 is not 64 bytes" );
static_assert( sizeof( eint1024 ) == 128, " eint1024 is not 128 bytes" );
static_assert( sizeof( eint2048 ) == 256, " eint2048 is not 256 bytes" );


// Test numeric limits (can be done at compile time)
template<class TYPE>
constexpr bool test_numeric_limits()
{
    constexpr TYPE zero( 0 );
    constexpr size_t bytes    = sizeof( TYPE );
    constexpr size_t digits   = 8 * bytes - 1;
    constexpr size_t digits10 = digits * 0.301029995663981;
    static_assert( bytes % 8 == 0 );
    static_assert( std::numeric_limits<TYPE>::is_specialized );
    static_assert( std::numeric_limits<TYPE>::is_specialized );
    static_assert( std::numeric_limits<TYPE>::is_signed );
    static_assert( std::numeric_limits<TYPE>::is_integer );
    static_assert( std::numeric_limits<TYPE>::is_exact );
    static_assert( !std::numeric_limits<TYPE>::has_infinity );
    static_assert( !std::numeric_limits<TYPE>::has_quiet_NaN );
    static_assert( !std::numeric_limits<TYPE>::has_signaling_NaN );
    static_assert( std::numeric_limits<TYPE>::has_denorm == std::denorm_absent );
    static_assert( !std::numeric_limits<TYPE>::has_denorm_loss );
    static_assert( std::numeric_limits<TYPE>::round_style == std::round_toward_zero );
    static_assert( !std::numeric_limits<TYPE>::is_iec559 );
    static_assert( std::numeric_limits<TYPE>::is_bounded );
    static_assert( !std::numeric_limits<TYPE>::is_modulo );
    static_assert( std::numeric_limits<TYPE>::digits == digits );
    static_assert( std::numeric_limits<TYPE>::digits10 == digits10 );
    static_assert( std::numeric_limits<TYPE>::max_digits10 == 0 );
    static_assert( std::numeric_limits<TYPE>::radix == 2 );
    static_assert( std::numeric_limits<TYPE>::min_exponent == 0 );
    static_assert( std::numeric_limits<TYPE>::min_exponent10 == 0 );
    static_assert( std::numeric_limits<TYPE>::max_exponent == 0 );
    static_assert( std::numeric_limits<TYPE>::max_exponent10 == 0 );
    static_assert( std::numeric_limits<TYPE>::traps );
    static_assert( !std::numeric_limits<TYPE>::tinyness_before );
    static_assert( std::numeric_limits<TYPE>::min() == std::numeric_limits<TYPE>::lowest() );
    static_assert( std::numeric_limits<TYPE>::max() > zero );
    static_assert( zero < std::numeric_limits<TYPE>::max() );
    static_assert( std::numeric_limits<TYPE>::min() < zero );
    static_assert( zero > std::numeric_limits<TYPE>::min() );
    return true;
}


// The main function
int main( int, char *[] )
{
    int N_errors = 0;

    // Print the size of the largest integer
    if ( sizeof( intmax_t ) != 8 )
        printp( "sizeof(intmax_t) = %u, optimization possible\n", (uint8_t) sizeof( intmax_t ) );

    // Test numeric_limits
    static_assert( test_numeric_limits<int64_t>() );
    static_assert( test_numeric_limits<eint64>() );
    static_assert( test_numeric_limits<eint128>() );
    static_assert( test_numeric_limits<eint256>() );
    static_assert( test_numeric_limits<eint512>() );
    static_assert( test_numeric_limits<eint1024>() );

    std::cout << eint64( 0 ).sign() << std::endl;
    std::cout << std::numeric_limits<eint64>::min().sign() << std::endl;

    // Create some basic variables
    std::cout << eint64( 1 ) << std::endl;
    std::cout << eint128( 1 ) << std::endl;
    std::cout << eint256( 1 ) << std::endl;
    std::cout << eint128( -1 ) << std::endl;
    std::cout << ( eint128( -1 ) * eint128( -1 ) ) << std::endl << std::endl;
    eint256 a;
    eint256 b( 12 );
    eint256 c( -13 );
    eint256 d( b );
    // eint256 e = c;
    // eint256 f("-15");

    // Check some very basic math
    if ( a.get_int() != 0 || b.get_int() != 12 || c.get_int() != -13 || d.get_int() != 12 ||
         ( -b ).get_int() != -12 || ( a + b ).get_int() != 12 || ( b + b ).get_int() != 24 ||
         ( b + c ).get_int() != -1 || ( c + c ).get_int() != -26 || ( b - c ).get_int() != 25 ||
         ( c - b ).get_int() != -25 ) {
        printp( "Failed basic math:\n" );
        printp( "   %i %i %i %i %i\n",
                a.get_int(),
                b.get_int(),
                c.get_int(),
                d.get_int(),
                ( -b ).get_int() );
        printp( "   %i %i %i\n", ( a + b ).get_int(), ( b + b ).get_int(), ( b + c ).get_int() );
        printp( "   %i %i %i\n", ( c + c ).get_int(), ( b - c ).get_int(), ( c - b ).get_int() );
        N_errors++;
    }

    // Check some shifts
    {
        eint128 tmp( 1 );
        tmp <<= 6;
        if ( tmp.get_int() != 64 || ( tmp << -1 ).get_int() != 32 ) {
            printp( "Failed leftshift\n" );
            N_errors++;
        }
        tmp = eint128( 512 );
        tmp >>= 6;
        if ( tmp.get_int() != 8 || ( tmp >> -1 ).get_int() != 16 || ( tmp >> 12 ).get_int() != 0 ) {
            printp( "Failed rightshift\n" );
            N_errors++;
        }
        eint256 tmp2     = eint256( 1 ) << 250;
        const double ans = 1.809251394333066e+75;
        if ( fabs( tmp2.get_double() - ans ) / ans > 1e-14 || ( tmp2 >> 245 ).get_int() != 32 ) {
            printp( "Failed << 250\n" );
            std::cout << "   " << tmp2 << std::endl;
            std::cout << "   " << tmp2.get_double() << std::endl;
            std::cout << "   " << ans << std::endl;
            std::cout << "   " << ( tmp2 >> 245 ).get_int() << std::endl;
            N_errors++;
        }
    }

    // Test some more math
    {
        AMP::extended::int64N<3> tmp1( "0x10000000000000000" ); // 2^64
        AMP::extended::int64N<3> tmp2 = -tmp1;
        AMP::extended::int64N<3> tmp3 = tmp1;
        tmp3.compliment();
        AMP::extended::int64N<3> tmp4 = 7 - tmp1;
        if ( tmp1.get_hex() != "0x000000000000000000000000000000010000000000000000" ||
             tmp2.get_hex() != "0xffffffffffffffffffffffffffffffff0000000000000000" ||
             tmp3.get_hex() != "0xffffffffffffffffffffffffffffffff0000000000000000" ||
             tmp4.get_hex() != "0xffffffffffffffffffffffffffffffff0000000000000007" ) {
            printp( "Failed test: 2^64, -2^64, -2^64, 7-2^64\n" );
            std::cout << "   " << tmp1 << std::endl;
            std::cout << "   " << tmp2 << std::endl;
            std::cout << "   " << tmp3 << std::endl;
            std::cout << "   " << tmp4 << std::endl << std::endl;
            N_errors++;
        }
    }


    // Test multiple by repeated multiplication of -1e8
    eint2048 x( -100000000 ), x2( 1 );
    printp( "%24.16e\n", x2.get_double() );
    double ans = 1.0;
    for ( int i = 0; i < 42; i++ ) {
        if ( i <= 38 && fabs( ( x2.get_double() - ans ) / ans ) > 1e-12 )
            N_errors++;
        x2 = x2 * x;
        ans *= -1e8;
        printp( "%24.16e\n", x2.get_double() );
    }
    std::cout << std::endl;

    // Test the performance of add
    int N_it = 10000;
    int N    = 1000;
    std::default_random_engine gen;
    std::uniform_int_distribution<int> dist( 0, 1000 );
    double t1 = time();
    for ( volatile int j = 0; j < N_it; j++ ) {
        int64_t y1( dist( gen ) );
        int64_t y2( dist( gen ) );
        for ( int i = 0; i < N; i++ )
            y1 += y2;
        NULL_USE( j );
        NULL_USE( y1 );
    }
    double t2 = time();
    double t3 = time();
    for ( volatile int j = 0; j < N_it; j++ ) {
        eint256 y5( dist( gen ) );
        eint256 y6( dist( gen ) );
        for ( int i = 0; i < N; i++ )
            y5 += y6;
        NULL_USE( j );
        NULL_USE( y5 );
    }
    double t4 = time();
    printp( "\nPerforamnce (add):\n" );
    printp( "   int64_t: %0.3f ns\n", 1e9 * ( t2 - t1 ) / ( N * N_it ) );
    printp( "   eint256: %0.3f ns\n", 1e9 * ( t4 - t3 ) / ( N * N_it ) );

    // Test the performance of multiply
    N  = 1000;
    t1 = time();
    for ( int j = 0; j < N_it; j++ ) {
        int64_t y1( dist( gen ) % 3 + 1 );
        int64_t y2( dist( gen ) % 3 + 1 );
        for ( int i = 0; i < N; i++ )
            y1 *= y2;
        NULL_USE( y1 );
    }
    t2 = time();
    t3 = time();
    for ( int j = 0; j < N_it; j++ ) {
        eint256 y5( dist( gen ) % 3 + 1 );
        eint256 y6( dist( gen ) % 3 + 1 );
        for ( int i = 0; i < N; i++ )
            y5 *= y6;
        NULL_USE( y5 );
    }
    t4 = time();
    printp( "Perforamnce (multiply):\n" );
    printp( "   int64_t: %0.3f ns\n", 1e9 * ( t2 - t1 ) / ( N * N_it ) );
    printp( "   eint256: %0.3f ns\n", 1e9 * ( t4 - t3 ) / ( N * N_it ) );
    printp( "\n" );

    // Run some known problems
    {
        eint128 tmp1( "0x9040ac62a509a917168" );
        eint128 tmp2( "0x152CDCE1C94F7300127" );
        if ( eint256( tmp1 ).get_double() != tmp1.get_double() ) {
            printp( "Failed conversionn\n" );
            N_errors++;
        }
        eint256 tmp3   = eint256( tmp1 ) * eint256( tmp2 );
        std::string s3 = tmp3.get_hex();
        if ( s3 != "0x000000000000000000000000000bee95b886e9fd2b0d774eecc702b2b919aed8" ) {
            printp( "Failed known multiplication\n" );
            N_errors++;
        }
    }

    if ( N_errors == 0 )
        printp( "All tests passed\n" );
    else
        printp( "Some tests failed\n" );
    return N_errors;
}
