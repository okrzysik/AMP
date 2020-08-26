#include "AMP/vectors/Scalar.h"


// Test auto creation of a Scalar
bool fun( const AMP::Scalar &x ) { return x.get<double>() != 0.0; }


// Test storing and getting a value (integer)
template<class TYPE>
bool testGet( TYPE x )
{
    // Store the scalar
    auto y    = AMP::Scalar( x );
    bool pass = true;
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
    auto c1 = AMP::Scalar( std::complex<float>( 3.0, 1.0 ) );
    pass    = pass && c1.get<std::complex<float>>() == std::complex<float>( 3.0, 1.0 );

    // Test copy
    auto c2 = c1;
    pass    = pass && c1.get<std::complex<float>>() == c2.get<std::complex<float>>();

    if ( pass )
        std::cout << "All tests passed\n";
    else
        std::cout << "Some tests failed\n";
    return pass ? 0 : 1;
}
