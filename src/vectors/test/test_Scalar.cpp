#include "AMP/vectors/Scalar.h"

#include <chrono>
#include <string_view>


// Test auto creation of a Scalar
bool fun( const AMP::Scalar &x ) { return x.get<double>() != 0.0; }
bool fun2( const std::any &x ) { return std::any_cast<size_t>( x ) != 0; }


// Test storing and getting a value (integer)
template<class TYPE>
bool testGet( TYPE x )
{
    // Store the scalar
    AMP::Scalar y = AMP::Scalar::create( x );
    bool pass     = true;
    // Test getting the scalar
    auto z = std::real( x );
    pass   = pass && y.get<int>() == static_cast<int>( z );
    pass   = pass && y.get<uint64_t>() == static_cast<uint64_t>( z );
    pass   = pass && y.get<float>() == static_cast<float>( z );
    pass   = pass && y.get<double>() == static_cast<double>( z );
    pass   = pass && y.get<std::complex<int>>() == std::complex<int>( z, 0.0 );
    pass   = pass && y.get<std::complex<float>>() == std::complex<float>( z, 0.0 );
    pass   = pass && y.get<std::complex<double>>() == std::complex<double>( z, 0.0 );
    return pass;
}


// Test basic arithmetic
template<class TYPE>
bool testArithmetic()
{
    TYPE a        = 17.2;
    TYPE b        = 3.7;
    AMP::Scalar x = AMP::Scalar::create( a );
    AMP::Scalar y = AMP::Scalar::create( b );
    bool pass     = true;
    // Test some different operations
    pass = pass && fabs( ( a - b ) - ( x - y ).get<TYPE>() ) < 1e-8;
    pass = pass && fabs( ( a + b ) - ( x + y ).get<TYPE>() ) < 1e-8;
    pass = pass && fabs( ( a * b ) - ( x * y ).get<TYPE>() ) < 1e-8;
    pass = pass && fabs( ( a / b ) - ( x / y ).get<TYPE>() ) < 1e-8;
    return pass;
}


// Test the performance
void testPerformance()
{
    size_t N = 1000000;

    auto start = std::chrono::system_clock::now();
    for ( size_t i = 1; i <= N; i++ ) {
        fun2( i );
    }
    auto stop  = std::chrono::system_clock::now();
    int64_t ns = std::chrono::duration_cast<std::chrono::nanoseconds>( stop - start ).count();
    printf( "Time to create std::any: %i ns\n", static_cast<int>( ns / N ) );

    start = std::chrono::system_clock::now();
    for ( size_t i = 1; i <= N; i++ )
        fun( i );
    stop = std::chrono::system_clock::now();
    ns   = std::chrono::duration_cast<std::chrono::nanoseconds>( stop - start ).count();
    printf( "Time to store/get value: %i ns\n", static_cast<int>( ns / N ) );
}


int main( int, char ** )
{
    bool pass = true;

    // Test some basic types
    pass = pass && testGet<char>( 'x' );
    pass = pass && testGet<int>( 1 );
    pass = pass && testGet<int64_t>( 2 );
    pass = pass && testGet<uint64_t>( 3 );
    pass = pass && testGet<float>( 4 );
    pass = pass && testGet<double>( 5 );
    pass = pass && testGet<long double>( 6 );
    pass = pass && testGet( std::complex<int>( 7.0, 0.0 ) );
    pass = pass && testGet( std::complex<float>( 8.0, 0.0 ) );
    pass = pass && testGet( std::complex<double>( 9.0, 0.0 ) );

    // Test passing values
    pass = pass && fun( 1 );
    pass = pass && fun( 2u );
    pass = pass && fun( 3.0f );
    pass = pass && fun( 4.0 );

    // Test complex
    auto c1 = AMP::Scalar::create( std::complex<float>( 3.0, 1.0 ) );
    pass    = pass && c1.get<std::complex<float>>() == std::complex<float>( 3.0, 1.0 );

    // Test copy
    auto c2 = c1;
    pass    = pass && c1.get<std::complex<float>>() == c2.get<std::complex<float>>();

    // Test arithmetic
    pass = pass && testArithmetic<double>();
    pass = pass && testArithmetic<int64_t>();
    pass = pass && testArithmetic<std::complex<double>>();

    // Test the performance
    testPerformance();

    if ( pass )
        std::cout << "All tests passed\n";
    else
        std::cout << "Some tests failed\n";
    return pass ? 0 : 1;
}
