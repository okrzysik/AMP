#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <ctime>
#include <set>
#include <vector>

#include "AMP/utils/AMPManager.h"
#include "AMP/utils/Utilities.h"
#include "AMP/utils/kdtree.h"

#include "ProfilerApp.h"

// Get a random double value between 0 and 1
double getRand()
{
    return ( (double) rand() / (double) RAND_MAX ) + 1e-7 * ( (double) rand() / (double) RAND_MAX );
}


// Factor a number into it's prime factors
std::vector<size_t> factor( size_t number )
{
    if ( number <= 3 )
        return std::vector<size_t>( 1, number );
    size_t n, n_max;
    bool factor_found;
    // Compute the maximum number of factors
    int N_primes_max = 1;
    n                = number;
    while ( n >>= 1 )
        ++N_primes_max;
    // Initialize n, factors
    n = number;
    std::vector<size_t> factors;
    while ( true ) {
        // Check if n is a trivial prime number
        if ( n == 2 || n == 3 || n == 5 ) {
            factors.push_back( (int) n );
            break;
        }
        // Check if n is divisible by 2
        if ( n % 2 == 0 ) {
            factors.push_back( 2 );
            n /= 2;
            continue;
        }
        // Check each odd number until a factor is reached
        n_max        = (size_t) floor( sqrt( (double) n ) );
        factor_found = false;
        for ( size_t i = 3; i <= n_max; i += 2 ) {
            if ( n % i == 0 ) {
                factors.push_back( i );
                n /= i;
                factor_found = true;
                break;
            }
        }
        if ( factor_found )
            continue;
        // No factors were found, the number must be prime
        factors.push_back( (int) n );
        break;
    }
    return factors;
}


// Compute the average distance between N points in a NDIM hypercube
// This is done by dividing the volume of the hypercube by the volume of the hyperspheres
double compute_avg_dist( int DIM, int N )
{
    double C        = 1.0;
    const double pi = 3.14159265359;
    if ( DIM == 1 ) {
        C = 1.0;
    } else if ( DIM == 2 ) {
        C = 4.0 / pi;
    } else if ( DIM == 3 ) {
        C = 6.0 / pi;
    } else if ( DIM == 4 ) {
        C = 32 / ( pi * pi );
    } else if ( DIM == 5 ) {
        C = 60.0 / pi * pi;
    }
    return pow( C / ( (double) N ), 1.0 / ( (double) DIM ) );
}


// Run the kdtree tests
int run_kdtree_test( int DIM, size_t Nx, size_t Ns )
{
    // Initialize the seed
    srand( static_cast<unsigned int>( time( nullptr ) ) );

    // Create the coordinates
    std::vector<std::vector<double>> points( DIM, std::vector<double>( Nx, 0.0 ) );
    for ( int d = 0; d < DIM; d++ ) {
        for ( size_t i = 0; i < Nx; i++ ) {
            points[d][i] = getRand();
        }
    }

    // Create the search points
    std::vector<std::vector<double>> search( DIM, std::vector<double>( Ns, 0.0 ) );
    for ( int d = 0; d < DIM; d++ ) {
        for ( size_t i = 0; i < Ns; i++ ) {
            search[d][i] = getRand();
        }
    }

    // Create the kdtree
    PROFILE_START( "create_tree" );
    const auto x = new const double *[DIM];
    for ( int d = 0; d < DIM; d++ )
        x[d] = &points[d][0];
    AMP::kdtree tree( DIM, Nx, x );
    delete[] x;
    PROFILE_STOP( "create_tree" );

    // Search for the local points
    PROFILE_START( "search_local" );
    bool error = false;
    double dist, xs[100], pos[100];
    for ( size_t i = 0; i < Nx; i++ ) {
        for ( int d = 0; d < DIM; d++ )
            xs[d] = points[d][i];
        size_t j = tree.find_nearest( xs, &dist, pos );
        if ( j != i || dist != 0.0 )
            error = true;
        for ( int d = 0; d < DIM; d++ )
            if ( pos[d] != xs[d] )
                error = true;
    }
    PROFILE_STOP( "search_local" );
    if ( error ) {
        printf( "Unable to correctly find local point in kdtree (%i,%i,%i)\n",
                DIM,
                (int) Nx,
                (int) Ns );
        return -1;
    }

    // Search for the unknown points
    PROFILE_START( "search_unknown" );
    error           = false;
    double dist_max = sqrt( (double) DIM ); // Maximum possible distance between two points
    double dist_avg = compute_avg_dist( DIM, (int) Nx ); // Average distance between two points
    for ( size_t i = 0; i < Ns; i++ ) {
        for ( int d = 0; d < DIM; d++ )
            xs[d] = search[d][i];
        size_t j = tree.find_nearest( xs, &dist, pos );
        if ( j >= Nx || dist <= 0.0 || dist > dist_max || dist > 100.0 * dist_avg )
            error = true;
    }
    PROFILE_STOP( "search_unknown" );
    if ( error ) {
        printf( "Unable to find neighbor to unknown point in kdtree (%i,%i,%i)\n",
                DIM,
                (int) Nx,
                (int) Ns );
        return -2;
    }

    return 0;
}


// Main
int main( int argc, char *argv[] )
{
    AMP::AMPManager::startup( argc, argv );

    int error = 0;
    PROFILE_ENABLE( 3 );
    PROFILE_ENABLE_TRACE();

    // Run a 1D test
    PROFILE_START( "1D kdtree" );
    error += run_kdtree_test( 1, 10, 10 );
    error += run_kdtree_test( 1, 10000, 1000 );
    PROFILE_STOP( "1D kdtree" );

    // Run a 2D test
    PROFILE_START( "2D kdtree" );
    error += run_kdtree_test( 2, 10, 10 );
    error += run_kdtree_test( 2, 100000, 10000 );
    PROFILE_STOP( "2D kdtree" );

    // Run a 3D test
    PROFILE_START( "3D kdtree" );
    error += run_kdtree_test( 3, 10, 10 );
    error += run_kdtree_test( 3, 100000, 100000 );
    PROFILE_STOP( "3D kdtree" );

    // Run a 5D test
    PROFILE_START( "5D kdtree" );
    error += run_kdtree_test( 5, 10, 10 );
    error += run_kdtree_test( 5, 10000, 2000 );
    PROFILE_STOP( "5D kdtree" );

    // Run a 10D test
    //    PROFILE_START( "10D kdtree" );
    //    error += run_kdtree_test( 10, 10, 10 );
    //    error += run_kdtree_test( 10, 10000, 2000 );
    //    PROFILE_STOP( "10D kdtree" );

    // Save the results
    PROFILE_SAVE( "test_kdtree" );
    if ( error == 0 )
        printf( "All tests passed\n" );
    else
        printf( "Some tests failed\n" );
    AMP::AMPManager::shutdown();
    return error;
}
