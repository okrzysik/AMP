#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <ctime>
#include <random>
#include <set>
#include <vector>

#include "AMP/utils/AMPManager.h"
#include "AMP/utils/UnitTest.h"
#include "AMP/utils/Utilities.h"
#include "AMP/utils/kdtree.h"

#include "ProfilerApp.h"


// Add the result to the tests
void checkResult( AMP::UnitTest &ut, bool pass, const std::string &msg )
{
    if ( pass )
        ut.passes( msg );
    else
        ut.failure( msg );
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
void run_kdtree_test( AMP::UnitTest &ut, int DIM, size_t Nx, size_t Ns )
{
    auto prefix = AMP::Utilities::stringf( "kdtree<%i,%i,%i>::", DIM, Nx, Ns );

    // Initialize the random number
    static std::random_device rd;
    static std::mt19937 gen( rd() );
    static std::uniform_real_distribution<double> dis( 0, 1 );

    // Create the coordinates
    std::vector<std::vector<double>> points( DIM, std::vector<double>( Nx, 0.0 ) );
    for ( int d = 0; d < DIM; d++ ) {
        for ( size_t i = 0; i < Nx; i++ )
            points[d][i] = dis( gen );
    }

    // Create the search points
    std::vector<std::vector<double>> search( DIM, std::vector<double>( Ns, 0.0 ) );
    for ( int d = 0; d < DIM; d++ ) {
        for ( size_t i = 0; i < Ns; i++ )
            search[d][i] = dis( gen );
    }

    // Create the kdtree
    PROFILE_START( "create_tree" );
    const auto x = new const double *[DIM];
    for ( int d = 0; d < DIM; d++ )
        x[d] = &points[d][0];
    AMP::kdtree tree( DIM, Nx, x );
    delete[] x;
    PROFILE_STOP( "create_tree" );

    // Check the bounding box
    auto box  = tree.box();
    bool pass = (int) box.size() == 2 * DIM;
    for ( int d = 0; d < DIM; d++ ) {
        pass = pass && box[2 * d + 0] >= 0.0;
        pass = pass && box[2 * d + 1] <= 1.0;
    }
    checkResult( ut, pass, prefix + "box" );

    // Search for the local points
    PROFILE_START( "search_local" );
    pass = true;
    double dist, xs[100], pos[100];
    for ( size_t i = 0; i < Nx; i++ ) {
        for ( int d = 0; d < DIM; d++ )
            xs[d] = points[d][i];
        size_t j = tree.find_nearest( xs, &dist, pos );
        if ( j != i || dist != 0.0 )
            pass = false;
        for ( int d = 0; d < DIM; d++ )
            if ( pos[d] != xs[d] )
                pass = false;
    }
    PROFILE_STOP( "search_local" );
    checkResult( ut, pass, prefix + "find_nearest (local)" );

    // Search for the unknown points
    PROFILE_START( "search_unknown" );
    pass            = true;
    double dist_max = sqrt( (double) DIM ); // Maximum possible distance between two points
    double dist_avg = compute_avg_dist( DIM, (int) Nx ); // Average distance between two points
    for ( size_t i = 0; i < Ns; i++ ) {
        for ( int d = 0; d < DIM; d++ )
            xs[d] = search[d][i];
        size_t j = tree.find_nearest( xs, &dist, pos );
        if ( j >= Nx || dist <= 0.0 || dist > dist_max || dist > 100.0 * dist_avg )
            pass = false;
    }
    PROFILE_STOP( "search_unknown" );
    checkResult( ut, pass, prefix + "find_nearest (unknown)" );
}


// Main
int main( int argc, char *argv[] )
{
    AMP::AMPManager::startup( argc, argv );
    AMP::UnitTest ut;

    PROFILE_ENABLE( 3 );
    PROFILE_ENABLE_TRACE();

    // Run a 1D test
    PROFILE_START( "1D kdtree" );
    run_kdtree_test( ut, 1, 10, 10 );
    run_kdtree_test( ut, 1, 10000, 1000 );
    PROFILE_STOP( "1D kdtree" );

    // Run a 2D test
    PROFILE_START( "2D kdtree" );
    run_kdtree_test( ut, 2, 10, 10 );
    run_kdtree_test( ut, 2, 100000, 10000 );
    PROFILE_STOP( "2D kdtree" );

    // Run a 3D test
    PROFILE_START( "3D kdtree" );
    run_kdtree_test( ut, 3, 10, 10 );
    run_kdtree_test( ut, 3, 100000, 100000 );
    PROFILE_STOP( "3D kdtree" );

    // Run a 5D test
    PROFILE_START( "5D kdtree" );
    run_kdtree_test( ut, 5, 10, 10 );
    run_kdtree_test( ut, 5, 10000, 2000 );
    PROFILE_STOP( "5D kdtree" );

    // Run a 10D test
    // PROFILE_START( "10D kdtree" );
    // run_kdtree_test( ut, 10, 10, 10 );
    // run_kdtree_test( ut, 10, 10000, 2000 );
    // PROFILE_STOP( "10D kdtree" );

    // Save the results
    PROFILE_SAVE( "test_kdtree" );
    ut.report();
    int N_errors = ut.NumFailGlobal();
    ut.reset();
    AMP::AMPManager::shutdown();
    return N_errors;
}
