#ifndef included_AMP_UtilityHelpers
#define included_AMP_UtilityHelpers

#include "AMP/utils/UnitTest.h"
#include "AMP/utils/Utilities.h"
#include "AMP/utils/Utilities.hpp"
#include "AMP/utils/enable_shared_from_this.h"
#include "AMP/utils/typeid.h"

#include "StackTrace/StackTrace.h"

#include <cmath>
#include <cstdio>
#include <limits>
#include <random>
#include <type_traits>


// Helper function to record pass/failure
#define PASS_FAIL( test, MSG ) \
    do {                       \
        if ( test )            \
            ut.passes( MSG );  \
        else                   \
            ut.failure( MSG ); \
    } while ( 0 )


// Subtract two size_t numbers, returning the absolute value
size_t abs_diff( size_t a, size_t b ) { return ( a >= b ) ? a - b : b - a; }


/****************************************************************
 * Checks approx_equal                                           *
 ****************************************************************/
template<class T>
void testApproxEqualInt( AMP::UnitTest &ut )
{
    using namespace AMP::Utilities;
    std::string type = typeid( T ).name();
    if ( approx_equal( 100000, 100000 ) && approx_equal_abs( 100000, 100000 ) &&
         !approx_equal( 100000, 100001 ) && !approx_equal( 100000, 100001 ) )
        ut.passes( "Integer (" + type + ") passes simple check." );
    else
        ut.failure( "Integer (" + type + ") passes simple check." );

    if ( approx_equal_abs( 100001, 100000, 1 ) && !approx_equal_abs( 100002, 100000, 1 ) )
        ut.passes( "Integer (" + type + ") passes close simple check." );
    else
        ut.failure( "Integer (" + type + ") passes close simple check." );
}
template<class T>
void testApproxEqual( AMP::UnitTest &ut )
{
    using namespace AMP::Utilities;
    std::string type = typeid( T ).name();
    constexpr T one  = static_cast<T>( 1.0 );
    const T closeTol = std::pow( std::numeric_limits<T>::epsilon(), (T) 0.8 );
    const T wrongTol = std::pow( std::numeric_limits<T>::epsilon(), (T) 0.7 );

    T mine      = one;
    T close_rel = mine * static_cast<T>( one + closeTol );
    T wrong_rel = mine * static_cast<T>( one + wrongTol );
    T close_abs = mine + closeTol;
    T wrong_abs = mine + wrongTol;
    if ( approx_equal( mine, close_rel ) && approx_equal_abs( mine, close_abs ) &&
         !approx_equal( mine, wrong_rel ) && !approx_equal_abs( mine, wrong_abs ) )
        ut.passes( type + " passes simple check near 1" );
    else
        ut.failure( type + " passes simple check near 1" );

    mine      = static_cast<T>( 1e-6 );
    close_rel = mine * static_cast<T>( one + closeTol );
    wrong_rel = mine * static_cast<T>( one + wrongTol );
    close_abs = mine + closeTol;
    wrong_abs = mine + wrongTol;
    if ( approx_equal( mine, close_rel ) && approx_equal_abs( mine, close_abs ) &&
         !approx_equal( mine, wrong_rel ) && !approx_equal_abs( mine, wrong_abs ) )
        ut.passes( type + " passes simple check near 1e-6" );
    else
        ut.failure( type + " passes simple check near 1e-6" );

    mine      = static_cast<T>( -1e-32 );
    close_rel = mine * static_cast<T>( one + closeTol );
    wrong_rel = mine * static_cast<T>( one + wrongTol );
    close_abs = mine + closeTol;
    wrong_abs = mine + wrongTol;
    if ( approx_equal( mine, close_rel ) && approx_equal_abs( mine, close_abs ) &&
         !approx_equal( mine, wrong_rel ) && !approx_equal_abs( mine, wrong_abs ) )
        ut.passes( type + " passes simple check near -1e-32" );
    else
        ut.failure( type + " passes simple check near -1e-32" );
}


/****************************************************************
 * Function to return the call stack                             *
 ****************************************************************/
std::vector<StackTrace::stack_info> get_call_stack() { return StackTrace::getCallStack(); }


/****************************************************************
 * Function to test the interpolants                             *
 ****************************************************************/
void test_interp( AMP::UnitTest &ut )
{
    const double a  = 1.0;
    const double bx = 1.0;
    const double by = -1.0;
    const double bz = 0.5;
    int Nx          = 20;
    int Ny          = 10;
    int Nz          = 5;
    std::vector<double> x( Nx, 0.0 );
    std::vector<double> y( Ny, 0.0 );
    std::vector<double> z( Nz, 0.0 );
    std::vector<double> f1( Nx, 0.0 );
    std::vector<double> f2( Nx * Ny, 0.0 );
    std::vector<double> f3( Nx * Ny * Nz, 0.0 );
    for ( int i = 0; i < Nx; i++ ) {
        x[i]  = ( (double) i ) / ( (double) ( Nx - 1 ) );
        f1[i] = a + bx * x[i];
        for ( int j = 0; j < Ny; j++ ) {
            y[j]           = ( (double) j ) / ( (double) ( Ny - 1 ) );
            f2[i + j * Nx] = a + bx * x[i] + by * y[j];
            for ( int k = 0; k < Nz; k++ ) {
                z[k]                         = ( (double) k ) / ( (double) ( Nz - 1 ) );
                f3[i + j * Nx + k * Nx * Ny] = a + bx * x[i] + by * y[j] + bz * z[k];
            }
        }
    }
    bool pass_linear    = true;
    bool pass_bilinear  = true;
    bool pass_trilinear = true;
    int Nix             = 100;
    int Niy             = 200;
    int Niz             = 50;
    for ( int i = 0; i < Nix; i++ ) {
        double xi = ( (double) i - 2 ) / ( (double) ( Nix - 5 ) );
        double fi = AMP::Utilities::linear( x, f1, xi );
        if ( !AMP::Utilities::approx_equal( fi, a + bx * xi, 1e-12 ) )
            pass_linear = false;
        for ( int j = 0; j < Niy; j++ ) {
            double yi = ( (double) j - 2 ) / ( (double) ( Niy - 5 ) );
            fi        = AMP::Utilities::bilinear( x, y, f2, xi, yi );
            if ( !AMP::Utilities::approx_equal( fi, a + bx * xi + by * yi, 1e-12 ) )
                pass_bilinear = false;
            for ( int k = 0; k < Niz; k++ ) {
                double zi = ( (double) k - 2 ) / ( (double) ( Niz - 5 ) );
                fi        = AMP::Utilities::trilinear( x, y, z, f3, xi, yi, zi );
                if ( !AMP::Utilities::approx_equal( fi, a + bx * xi + by * yi + bz * zi, 1e-12 ) )
                    pass_trilinear = false;
            }
        }
    }
    PASS_FAIL( pass_linear, "Linear interpolation" );
    PASS_FAIL( pass_bilinear, "Bi-linear interpolation" );
    PASS_FAIL( pass_trilinear, "Tri-linear interpolation" );
}


/****************************************************************
 * Test enable_shared_from_this                                  *
 ****************************************************************/
class dummy : public AMP::enable_shared_from_this<dummy>
{
public:
    std::shared_ptr<dummy> getPtr() { return shared_from_this(); }
};
static inline bool test_shared_from_this_pointer( const std::shared_ptr<dummy> &p1 )
{
    bool pass = p1.use_count() == 1;
    int count1, count2;
    auto p2 = p1->getPtr();
    count1  = p1.use_count();
    count2  = p2.use_count();
    pass    = pass && count1 == 2 && count2 == 2;
    std::shared_ptr<dummy> p3( p1 );
    count1 = p2.use_count();
    count2 = p3.use_count();
    pass   = pass && count1 == 3 && count2 == 3;
    p2.reset();
    count1  = p2.use_count();
    count2  = p3.use_count();
    pass    = pass && count1 == 0 && count2 == 2;
    auto p4 = p3->getPtr();
    count1  = p3.use_count();
    count2  = p4.use_count();
    pass    = pass && count1 == 3 && count2 == 3;
    std::shared_ptr<dummy> p5( p3.get(), []( void * ) {} );
    count1 = p3.use_count();
    count2 = p5.use_count();
    pass   = pass && count1 == 3 && count2 == 1;
    p5.reset();
    count1 = p3.use_count();
    count2 = p5.use_count();
    pass   = pass && p3.use_count() == 3 && p5.use_count() == 0;
    return pass;
}
void test_shared_from_this( AMP::UnitTest &ut )
{
    bool pass = true;
    try {
        auto ptr = std::make_shared<dummy>();
        pass     = test_shared_from_this_pointer( ptr );
    } catch ( ... ) {
        pass = false;
    }
    PASS_FAIL( pass, "shared_from_this 1" );
    try {
        auto *p1 = new dummy;
        auto p2  = p1->getPtr();
        pass     = test_shared_from_this_pointer( p2 );
    } catch ( ... ) {
        pass = false;
    }
    PASS_FAIL( pass, "shared_from_this 2" );
}


/****************************************************************
 *  Test quicksort/quickselect                                    *
 ****************************************************************/
void testQuickSort( AMP::UnitTest &ut, std::vector<int> &data1, const std::string &str )
{
    auto data2 = data1;
    auto data3 = data1;
    auto data4 = data1;
    double t1  = AMP::Utilities::time();
    AMP::Utilities::quicksort( data1 );
    double t2 = AMP::Utilities::time();
    std::sort( data2.begin(), data2.end() );
    double t3 = AMP::Utilities::time();
    std::sort( &data3[0], &data3[0] + data3.size() );
    double t4 = AMP::Utilities::time();
    bool pass = true;
    for ( size_t i = 0; i < data1.size(); i++ ) {
        if ( data1[i] != data2[i] )
            pass = false;
    }
    PASS_FAIL( pass, "quicksort sorts correctly: " + str );
    std::cout << "quicksort:" << str << " = " << t2 - t1 << ", std::sort = " << t3 - t2
              << ", std::sort(2) = " << t4 - t3 << std::endl;
}
void testQuickSort( AMP::UnitTest &ut )
{
    size_t N = 500000;
    std::vector<int> data( N, 31 );
    testQuickSort( ut, data, "identical" );
    for ( size_t i = 0; i < N; i++ )
        data[i] = i;
    testQuickSort( ut, data, "sorted" );
    std::random_device rd;
    std::mt19937 gen( rd() );
    std::uniform_int_distribution<int> dist( 1, 10000000 );
    for ( size_t i = 0; i < N; i++ )
        data[i] = dist( gen );
    testQuickSort( ut, data, "random" );
}
void testQuickSelect( AMP::UnitTest &ut )
{
    size_t N = 500000;
    std::random_device rd;
    std::mt19937 gen( rd() );
    std::uniform_int_distribution<int> dist( 1, 10000000 );
    std::vector<int> data( N, 31 );
    for ( size_t i = 0; i < N; i++ )
        data[i] = dist( gen );
    std::sort( data.begin(), data.end() );
    size_t N_it = 200;
    bool pass   = true;
    auto t1     = AMP::Utilities::time();
    for ( size_t i = 0; i < N_it; i++ ) {
        size_t k = dist( gen ) % N;
        auto v   = AMP::Utilities::quickselect( data.size(), data.data(), k );
        pass     = v == data[k];
    }
    double t = AMP::Utilities::time() - t1;
    if ( pass )
        ut.passes( "quickselect" );
    else
        ut.failure( "quickselect" );
    std::cout << "quickselect = " << t / N_it << std::endl;
}


/****************************************************************
 * Test precision                                                *
 ****************************************************************/
template<class TYPE>
constexpr int calculateDigits()
{
    int d = 40;
    while ( true ) {
        TYPE x  = 1;
        TYPE y  = std::pow( 10.0, -d );
        TYPE z  = std::pow( 10.0, d );
        TYPE x2 = x + y;
        TYPE x3 = x2 - x;
        if ( std::abs( ( x3 - y ) * z ) < 0.5 )
            break;
        d--;
    }
    return d;
}
template<class T>
void test_precision( [[maybe_unused]] AMP::UnitTest &ut )
{
    auto type = AMP::getTypeID<T>();
    if constexpr ( std::is_integral_v<T> ) {
        static_assert( std::numeric_limits<T>::digits >= sizeof( T ) - 1 );
        auto min = static_cast<double>( std::numeric_limits<T>::min() );
        auto max = static_cast<double>( std::numeric_limits<T>::max() );
        printf( "<%s>: Integer %e-%e\n", type.name, min, max );
    } else if constexpr ( std::is_floating_point_v<T> ) {
        int digits  = std::numeric_limits<T>::digits10;
        int digits2 = calculateDigits<T>();
        bool match  = std::abs( digits - digits2 ) <= 1;
        printf( "<%s>: Floating Point %i/%i\n", type.name, digits, digits2 );
        if constexpr ( std::is_same_v<T, float> ) {
            PASS_FAIL( digits == 6 && match, "<float> matches IEEE" );
        } else if constexpr ( std::is_same_v<T, double> ) {
            PASS_FAIL( digits == 15 && match, "<double> matches IEEE" );
        } else if constexpr ( std::is_same_v<T, long double> ) {
            PASS_FAIL( digits >= 18 && match, "<long double> matches expected" );
        } else {
            ut.failure( "Unknown type" );
        }
    } else {
        static_assert( !std::is_same_v<T, T>, "Invalid type" );
    }
}


/********************************************************
 *  Test typeid                                          *
 ********************************************************/
template<class T>
void check( const std::string &name, AMP::UnitTest &ut, bool expected )
{
    auto type = AMP::getTypeID<T>();
    bool pass =
        std::string_view( type.name ) == name && type.bytes == sizeof( T ) && type.hash != 0;
    if ( pass )
        ut.passes( "typeid: " + name );
    else if ( expected )
        ut.expected_failure( "typeid: " + name + " --- " + std::string( type.name ) );
    else
        ut.failure( "typeid: " + name + " --- " + std::string( type.name ) );
}
void testTypeID( AMP::UnitTest &ut )
{
    check<std::shared_ptr<double>>( "std::shared_ptr<double>", ut, false ); // Fails windows
    check<double *>( "double*", ut, true );                                 // Fails clang-16
    check<const double *>( "const double*", ut, true );                     // Fails clang-16
    check<double const *>( "const double*", ut, true );                     // Fails clang-16
}


#endif
