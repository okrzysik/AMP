#include "AMP/utils/AMPManager.h"
#include "AMP/utils/DelaunayInterpolation.h"
#include "AMP/utils/DelaunayTessellation.h"
#include "AMP/utils/NearestPairSearch.h"
#include "AMP/utils/UnitTest.h"

#include "test_DelaunayInterpolation.h"

#include "ProfilerApp.h"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <ctime>
#include <limits>
#include <memory>
#include <random>
#include <stdexcept>
#include <string>
#include <typeinfo>
#include <vector>


#define printp printf

#define NDIM_MAX 3 // The maximum number of dimensions supported (currently 3)

#define ASSERT( EXP )                                                     \
    do {                                                                  \
        if ( !( EXP ) ) {                                                 \
            char tmp[1000];                                               \
            sprintf( tmp,                                                 \
                     "Failed assertion: %s\n   File: %s\n:   Line: %i\n", \
                     #EXP,                                                \
                     __FILE__,                                            \
                     __LINE__ );                                          \
            throw std::logic_error( tmp );                                \
        }                                                                 \
                                                                          \
    } while ( 0 )


// Helper function to check if two numbers are approximately equal
inline bool approx_equal( double x, double y, double tol = 1e-8 )
{
    return 2.0 * fabs( x - y ) <= tol * fabs( x + y );
}


// Helper function to get the class name
template<class TYPE>
const char *getName();
template<>
const char *getName<int>()
{
    return "int";
}
template<>
const char *getName<double>()
{
    return "double";
}


// Helper function to write points to a file
template<class TYPE>
void writePoints( const char *filename, int ndim, const std::vector<TYPE> &x );
template<>
void writePoints<double>( const char *filename, int ndim, const std::vector<double> &x )
{
    FILE *fid;
    fid = fopen( filename, "wb" );
    fprintf(
        fid, "%i points in %iD in double precision\n", static_cast<int>( x.size() / ndim ), ndim );
    fwrite( &x[0], sizeof( double ), x.size(), fid );
    fclose( fid );
    printp( "  Failed points written to %s\n", filename );
}
template<>
void writePoints<int>( const char *filename, int ndim, const std::vector<int> &x )
{
    FILE *fid;
    fid = fopen( filename, "wb" );
    fprintf(
        fid, "%i points in %iD in int precision\n", static_cast<int>( x.size() / ndim ), ndim );
    fwrite( &x[0], sizeof( int ), x.size(), fid );
    fclose( fid );
    printp( "  Failed points written to %s\n", filename );
}


// Helper function to read points from a file
std::pair<int, std::vector<double>> readPoints( const char *filename )
{
    FILE *fid;
    fid = fopen( filename, "rb" );
    if ( fid == nullptr )
        throw std::logic_error( "Input file not found" );
    char tline[200];
    char *tmp1 = fgets( tline, 200, fid );
    ASSERT( tmp1 == tline );
    std::string header( tline );
    size_t i = header.find( " points in " );
    ASSERT( i != std::string::npos );
    std::string tmp = header.substr( 0, i );
    int N           = atoi( tmp.c_str() );
    i               = header.find( "D in " );
    ASSERT( i != std::string::npos );
    tmp   = header.substr( i - 2, 2 );
    int D = atoi( tmp.c_str() );
    ASSERT( N > 0 && D > 0 );
    std::vector<double> x( D * N, 0.0 );
    size_t tmp2 = fread( &x[0], sizeof( double ), N * D, fid );
    ASSERT( tmp2 == (size_t) N * D );
    fclose( fid );
    return std::pair<int, std::vector<double>>( D, x );
}


// Test for the nearest pair search
template<class TYPE>
void testPointSearch( AMP::UnitTest *ut, int ndim, const std::vector<TYPE> &x )
{
    // Test the nearest point search
    PROFILE_START( "Test point search 1" );
    double dist1 = 1e100;
    std::pair<int, int> index1( -1, -1 );
    for ( size_t i = 0; i < x.size() / ndim; i++ ) {
        for ( size_t j = i + 1; j < x.size() / ndim; j++ ) {
            double dist2 = 0.0;
            for ( int d = 0; d < ndim; d++ ) {
                auto tmp = static_cast<double>( x[d + i * ndim] - x[d + j * ndim] );
                dist2 += tmp * tmp;
            }
            if ( dist2 < dist1 ) {
                dist1  = dist2;
                index1 = std::pair<int, int>( (int) i, (int) j );
            }
        }
    }
    PROFILE_STOP( "Test point search 1" );
    PROFILE_START( "Test point search 2" );
    std::pair<int, int> index2;
    if ( ndim == 1 )
        index2 = AMP::find_min_dist<1, TYPE>( (int) x.size(), &x[0] );
    else if ( ndim == 2 )
        index2 = AMP::find_min_dist<2, TYPE>( (int) x.size() / 2, &x[0] );
    else if ( ndim == 3 )
        index2 = AMP::find_min_dist<3, TYPE>( (int) x.size() / 3, &x[0] );
    else if ( ndim == 4 )
        index2 = AMP::find_min_dist<4, TYPE>( (int) x.size() / 4, &x[0] );
    else if ( ndim == 5 )
        index2 = AMP::find_min_dist<5, TYPE>( (int) x.size() / 5, &x[0] );
    double dist2 = 0.0;
    for ( int d = 0; d < ndim; d++ ) {
        auto tmp = static_cast<double>( x[d + index2.first * ndim] - x[d + index2.second * ndim] );
        dist2 += tmp * tmp;
    }
    PROFILE_STOP( "Test point search 2" );
    char msg[100];
    sprintf( msg,
             "Test point search (%i,%i,%s)",
             ndim,
             static_cast<int>( x.size() ) / ndim,
             getName<TYPE>() );
    if ( approx_equal( dist1, dist2, 1e-12 ) ) {
        ut->passes( msg );
    } else {
        char msg2[200];
        sprintf( msg2,
                 "%s:  (%i,%i,%e) - (%i,%i,%e)",
                 msg,
                 index1.first,
                 index1.second,
                 dist1,
                 index2.first,
                 index2.second,
                 dist2 );
        ut->failure( msg2 );
        writePoints( "failed_points", ndim, x );
    }
}


// Function to create and test the construction of the tessellation
template<class TYPE>
std::shared_ptr<AMP::DelaunayInterpolation<TYPE>>
createAndTestDelaunayInterpolation( AMP::UnitTest *ut, int ndim, const std::vector<TYPE> &x )
{
    size_t N = x.size() / ndim;
    char tmp[32];
    sprintf( tmp, "(%i,%i,%s)", ndim, static_cast<int>( x.size() ) / ndim, getName<TYPE>() );
    std::string msg( tmp );

    // Create the tessellation
    PROFILE_START( "Create Tessellation" );
    auto data      = std::make_shared<AMP::DelaunayInterpolation<TYPE>>( ndim );
    int error_code = data->create_tessellation( (int) x.size() / ndim, &x[0] );
    size_t N_tri   = data->get_N_tri();
    PROFILE_STOP( "Create Tessellation" );
    if ( error_code == 0 && N_tri > 0 ) {
        ut->passes( "Created tessellation " + msg );
    } else {
        ut->failure( "Created tessellation " + msg );
        return nullptr;
    }

    // Check the triangle neighbors
    int *tri_nab = data->get_tri_nab( 0 );
    bool pass    = true;
    for ( size_t i = 0; i < N_tri; i++ ) {
        for ( int d = 0; d < ndim + 1; d++ ) {
            if ( tri_nab[d + i * ( ndim + 1 )] == (int) i || tri_nab[d + i * ( ndim + 1 )] < -1 ||
                 tri_nab[d + i * ( ndim + 1 )] >= (int) N_tri )
                pass = false;
        }
    }
    delete[] tri_nab;
    if ( !pass ) {
        ut->failure( "Triangle neighbors are invalid " + msg );
        return nullptr;
    }

    // Check the current memory usage
    PROFILE_START( "Check memory usage", 1 );
    size_t bytes  = data->memory_usage();
    size_t bytes1 = data->memory_usage( 1 );
    size_t bytes2 = data->memory_usage( 2 );
    size_t bytes3 = data->memory_usage( 3 );
    size_t bytes4 = data->memory_usage( 4 );
    if ( bytes1 !=
         sizeof( AMP::DelaunayInterpolation<TYPE> ) + ( ndim + 1 ) * N_tri * sizeof( int ) )
        ut->failure( "memory_usage(1) fails " + msg );
    if ( bytes2 != bytes1 + ndim * N * sizeof( TYPE ) )
        ut->failure( "memory_usage(2) fails " + msg );
    if ( bytes3 <= bytes1 || bytes3 > 10 * bytes1 )
        ut->failure( "memory_usage(3) fails " + msg );
    if ( bytes4 != bytes3 + bytes2 - bytes1 )
        ut->failure( "memory_usage(4) fails " + msg );
    if ( bytes < bytes1 || bytes > 1.1 * bytes4 )
        ut->failure( "memory_usage() fails " + msg );
    PROFILE_STOP( "Check memory usage", 1 );

    // Copy the tessellation
    PROFILE_START( "Copy tessellation", 1 );
    int N2, N_tri2;
    auto x2  = new TYPE[N * ndim];
    auto tri = new int[N_tri * ( ndim + 1 )];
    data->copy_tessellation( &N2, x2, 0, &N_tri2, tri );
    bool error = static_cast<size_t>( N2 ) != N || static_cast<size_t>( N_tri2 ) != N_tri;
    if ( error == 0 ) {
        for ( size_t i = 0; i < ndim * N; i++ ) {
            if ( x[i] != x2[i] )
                error = true;
        }
        for ( size_t i = 0; i < ( ndim + 1 ) * N_tri; i++ ) {
            if ( tri[i] < 0 || tri[i] >= N2 )
                error = true;
        }
    }
    AMP::DelaunayInterpolation<TYPE> data2( ndim );
    data2.create_tessellation( (int) N, &x[0], 0, (int) N_tri, tri );
    if ( data2.get_N_tri() != N_tri )
        error = true;
    if ( bytes1 != data2.memory_usage( 1 ) || bytes2 != data2.memory_usage( 2 ) ||
         bytes3 != data2.memory_usage( 3 ) || bytes4 != data2.memory_usage( 4 ) )
        error = true;
    delete[] x2;
    delete[] tri;
    data2.create_tessellation( 0, nullptr, 0, 0, nullptr );
    if ( data2.memory_usage() != sizeof( AMP::DelaunayInterpolation<TYPE> ) )
        error = true;
    if ( !error )
        ut->passes( "Copy of tessellation " + msg );
    else
        ut->failure( "Copy of tessellation " + msg );
    PROFILE_STOP( "Copy tessellation", 1 );

    // Check the behavior of get_circumsphere and test_in_circumsphere
    PROFILE_START( "Check circumsphere", 1 );
    {
        bool pass = true;
        int *tri  = data->get_tri( 0 );
        TYPE x1[NDIM_MAX * ( NDIM_MAX + 1 )];
        double R, c[NDIM_MAX], xi[NDIM_MAX], x2[NDIM_MAX * ( NDIM_MAX + 1 )];
        for ( size_t i = 0; i < N_tri; i++ ) {
            for ( int d1 = 0; d1 < ndim + 1; d1++ ) {
                int k = tri[d1 + i * ( ndim + 1 )];
                for ( int d2 = 0; d2 < ndim; d2++ ) {
                    x1[d2 + d1 * ndim] = x[d2 + k * ndim];
                    x2[d2 + d1 * ndim] = x[d2 + k * ndim];
                }
            }
            AMP::DelaunayTessellation::get_circumsphere( ndim, x1, R, c );
            if ( R < 0 ) {
                pass = false;
            }
            for ( int j = 0; j < ndim + 1; j++ ) {
                double dist = 0.0;
                for ( int d = 0; d < ndim; d++ )
                    dist += ( x2[d + j * ndim] - c[d] ) * ( x2[d + j * ndim] - c[d] );
                if ( fabs( sqrt( dist ) - R ) / R > 1e-7 )
                    pass = false;
            }
            for ( int j = 0; j < ndim; j++ ) {
                for ( int d = 0; d < ndim; d++ )
                    xi[d] = c[d];
                xi[j]     = c[j] - 10.0 * R;
                int test1 = AMP::DelaunayTessellation::test_in_circumsphere( ndim, x2, xi, 0 );
                xi[j]     = c[j] - 0.1 * R;
                int test2 = AMP::DelaunayTessellation::test_in_circumsphere( ndim, x2, xi, 0 );
                xi[j]     = c[j] + 10.0 * R;
                int test3 = AMP::DelaunayTessellation::test_in_circumsphere( ndim, x2, xi, 0 );
                xi[j]     = c[j] + 0.1 * R;
                int test4 = AMP::DelaunayTessellation::test_in_circumsphere( ndim, x2, xi, 0 );
                if ( test1 != -1 || test2 != 1 || test3 != -1 || test4 != 1 ) {
                    pass = false;
                }
            }
        }
        delete[] tri;
        if ( pass )
            ut->passes( "get_circumsphere and test_in_circumsphere " + msg );
        else
            ut->expected_failure( "get_circumsphere and test_in_circumsphere " + msg );
    }
    PROFILE_STOP( "Check circumsphere", 1 );

    // Check that the volume of each triangle is positive and > 0
    PROFILE_START( "Check volume", 1 );
    {
        int *tri       = data->get_tri( 0 );
        double vol_min = 1e100;
        double vol_max = 0.0;
        bool neg_vol   = false;
        for ( size_t i = 0; i < N_tri; i++ ) {
            int tri2[NDIM_MAX + 1];
            double x2[NDIM_MAX * ( NDIM_MAX + 1 )];
            for ( int d1 = 0; d1 < ndim + 1; d1++ ) {
                tri2[d1] = tri[d1 + i * ( ndim + 1 )];
                for ( int d2 = 0; d2 < ndim; d2++ )
                    x2[d2 + d1 * ndim] = x[d2 + tri2[d1] * ndim];
            }
            double vol = AMP::DelaunayTessellation::calc_volume( ndim, x2 );
            if ( vol < 0 ) {
                neg_vol = true;
                vol     = -vol;
            }
            vol_min = std::min( vol_min, vol );
            vol_max = std::max( vol_max, vol );
        }
        if ( !neg_vol && vol_min > 0.0 && vol_min >= 1e-10 * vol_max )
            ut->passes( "Tessellation volume is valid " + msg );
        else
            ut->failure( "Tessellation volume is invalid " + msg );
        delete[] tri;
        PROFILE_STOP( "Check volume", 1 );
    }

    // Perform a rigorous check of the tessellation by checking that no point in the
    // tessellation lies withing the circumsphere of any triangle
    // Note: This is an N^2 test and will only be run for "reasonable" problem sizes
    if ( N <= 5000 ) {
        PROFILE_START( "Check tessellation", 1 );
        int *tri               = data->get_tri( 0 );
        bool pass_circumsphere = true;
        for ( size_t i = 0; i < N_tri; i++ ) {
            int tri2[NDIM_MAX + 1];
            TYPE x2[NDIM_MAX * ( NDIM_MAX + 1 )], xi[NDIM_MAX];
            for ( int d1 = 0; d1 < ndim + 1; d1++ ) {
                tri2[d1] = tri[d1 + i * ( ndim + 1 )];
                for ( int d2 = 0; d2 < ndim; d2++ )
                    x2[d2 + d1 * ndim] = x[d2 + tri2[d1] * ndim];
            }
            for ( size_t j = 0; j < N; j++ ) {
                for ( int d = 0; d < ndim; d++ )
                    xi[d] = x[d + j * ndim];
                int test = AMP::DelaunayTessellation::test_in_circumsphere( ndim, x2, xi, 1e-8 );
                if ( test == 1 ) {
                    pass_circumsphere = false;
                    break;
                }
            }
        }
        delete[] tri;
        if ( pass_circumsphere )
            ut->passes( "Tessellation is valid " + msg );
        else
            ut->failure( "Tessellation is invalid " + msg );
        PROFILE_STOP( "Check tessellation", 1 );
    }

    // Set the storage level to 2 and check that the memory usage changed
    // Note: this does not apply to 1d
    PROFILE_START( "Change storage", 1 );
    int *tri_nab1 = data->get_tri_nab( 0 );
    data->set_storage_level( 2 );
    size_t bytes_new = data->memory_usage();
    if ( bytes_new == bytes2 && bytes_new != bytes )
        ut->passes( "Storage level changed " + msg );
    else
        ut->failure( "Storage level changed " + msg );
    data->set_storage_level( 4 );
    int *tri_nab2 = data->get_tri_nab( 0 );
    pass          = true;
    for ( size_t i = 0; i < N_tri * ( ndim + 1 ); i++ ) {
        if ( tri_nab1[i] != tri_nab2[i] )
            pass = false;
    }
    if ( !pass )
        ut->failure( "Triangle neighbors changed during call to change memory usage " + msg );
    delete[] tri_nab1;
    delete[] tri_nab2;
    PROFILE_STOP( "Change storage", 1 );

    return data;
}


// Create f and g for a given problem
template<class TYPE>
std::string initialize_problem( int p,
                                int ndim,
                                size_t N,
                                const TYPE *x0,
                                const double *xrange,
                                double *f,
                                double *g,
                                double &tol_grad,
                                double &tol_linear,
                                double &tol_cubic,
                                double &tol_cubic_grad,
                                double &tol_cubic_extrap1,
                                double &tol_cubic_extrap2 )
{
    std::string problem;
    tol_grad          = 1e-8;  // Tolerance for gradient calculation
    tol_linear        = 1e-8;  // Tolerance for linear interpolation
    tol_cubic         = 1e-10; // Tolerance for cubic interpolation (interior)
    tol_cubic_grad    = 1e-8;  // Tolerance for cubic interpolation (interior)
    tol_cubic_extrap1 = 1e-10; // Tolerance for linear extrapolation with cubic interpolation
    tol_cubic_extrap2 = 1e-10; // Tolerance for quadratic extrapolation with cubic interpolation
    // Set f and g for each problem
    double grad_max[10] = { 0 };
    if ( p == 0 ) {
        // Constant profile
        problem = std::string( "f(x) = 2.5" );
        for ( size_t i = 0; i < N; i++ )
            f[i] = 2.5;
        memset( g, 0, N * ndim * sizeof( double ) );
    } else if ( p == 1 ) {
        // Linear profile
        problem               = std::string( "f(x) = x-y+2*z" );
        const double df_dx[3] = { 1.0, -1.0, 2.0 };
        for ( size_t i = 0; i < N; i++ ) {
            f[i] = 0.0;
            for ( int d = 0; d < ndim; d++ ) {
                f[i] += df_dx[d] * x0[d + i * ndim] / xrange[d];
                g[d + i * ndim] = df_dx[d] / xrange[d];
            }
        }
        for ( int d = 0; d < ndim; d++ )
            grad_max[d] = df_dx[d] / xrange[d];
    } else if ( p == 2 ) {
        // Quadratic profile
        problem           = std::string( "f(x) = 1+x-y+2*z+0.1*x^2+0.2*y^2-0.3*z^2" );
        double dx2        = 1.0 / pow( 1.0 * N, 1.0 / ndim );
        tol_grad          = 2.0 * dx2;
        tol_linear        = std::min( 4 * pow( 3.0, ndim ) * dx2, 0.5 );
        tol_cubic         = 1e-8;
        tol_cubic_grad    = 1e-6;
        tol_cubic_extrap1 = tol_linear;
        tol_cubic_extrap2 = 1e-8;
        for ( size_t i = 0; i < N; i++ ) {
            double x2[3] = { 0, 0, 0 }, g2[3];
            for ( int d = 0; d < ndim; d++ )
                x2[d] = x0[d + i * ndim];
            double x = x2[0] / xrange[0];
            double y = x2[1] / xrange[1];
            double z = x2[2] / xrange[2];
            f[i]     = 1 + x - y + 2 * z + 0.1 * x * x + 0.2 * y * y - 0.3 * z * z;
            g2[0]    = 1 + 0.2 * x;
            g2[1]    = -1 + 0.4 * y;
            g2[2]    = 2 - 0.6 * z;
            for ( int d = 0; d < ndim; d++ )
                g[d + i * ndim] = g2[d] / xrange[d];
        }
        const double grad_max0[3] = { 1.2, 1, 2 };
        for ( int d = 0; d < ndim; d++ )
            grad_max[d] = grad_max0[d] / xrange[d];
        if ( N < 10 ) {
            tol_grad          = 1e100;
            tol_linear        = 1e100;
            tol_cubic         = 1e100;
            tol_cubic_grad    = 1e100;
            tol_cubic_extrap1 = 1e100;
            tol_cubic_extrap2 = 1e100;
        }
    } else {
        throw std::logic_error( "Unknown problem" );
    }
    tol_grad *=
        sqrt( grad_max[0] * grad_max[0] + grad_max[1] * grad_max[1] + grad_max[2] * grad_max[2] );
    tol_cubic_grad *=
        sqrt( grad_max[0] * grad_max[0] + grad_max[1] * grad_max[1] + grad_max[2] * grad_max[2] );
    tol_grad       = std::max( tol_grad, 1e-12 );
    tol_cubic_grad = std::max( tol_cubic_grad, 1e-12 );
    return problem;
}


template<class TYPE>
std::vector<double> convert_to_double( const std::vector<TYPE> &x0 );
template<>
std::vector<double> convert_to_double<double>( const std::vector<double> &x0 )
{
    return x0;
}
template<>
std::vector<double> convert_to_double<int>( const std::vector<int> &x0 )
{
    std::vector<double> x( x0.size(), 0.0 );
    for ( size_t i = 0; i < x0.size(); i++ )
        x[i] = static_cast<double>( x0[i] ) / static_cast<double>( R_INT );
    return x;
}
template<class TYPE>
void testInterpolation( AMP::UnitTest *ut,
                        int ndim,
                        const std::vector<TYPE> &x,
                        bool check_extrap = true )
{
    check_extrap = false; // Extrapolation is having major issues (disable for now)

    size_t N = x.size() / ndim;
    if ( N == 0 )
        return;

    // Check the points
    if ( x.size() / ndim < 10000 )
        testPointSearch( ut, ndim, x );

    // Create the tessellation
    auto data = createAndTestDelaunayInterpolation<TYPE>( ut, ndim, x );
    if ( data == nullptr )
        return;
    size_t N_tri = data->get_N_tri();

    // Create a grid for interpolation at the center of each triangle
    size_t N1  = N_tri;
    auto order = new size_t[N1];
    for ( size_t i = 0; i < N1; i++ )
        order[i] = i;
    std::shuffle( order, &order[N1], std::mt19937( std::random_device()() ) );
    auto x1        = new double[ndim * N1];
    int *triangles = data->get_tri( 0 );
    for ( size_t i = 0; i < N1; i++ ) {
        for ( int d = 0; d < ndim; d++ ) {
            x1[d + i * ndim] = 0.0;
            for ( int j = 0; j <= ndim; j++ ) {
                int k = triangles[j + order[i] * ( ndim + 1 )];
                x1[d + i * ndim] += x[d + k * ndim];
            }
            x1[d + i * ndim] /= ndim + 1;
        }
    }

    // Create a second grid for interpolation
    double dx[3]     = { 0, 0, 0 };
    double xrange[3] = { 1, 1, 1 };
    double tol       = std::min( 1e-2, 1.0 / pow( 1.0 * N, 1.0 / ndim ) );
    for ( int d = 0; d < ndim; d++ ) {
        double x_min = 1e100;
        double x_max = -1e100;
        for ( size_t i = 0; i < N; i++ ) {
            x_min = std::min<double>( x_min, x[d + i * ndim] );
            x_max = std::max<double>( x_max, x[d + i * ndim] );
        }
        dx[d]     = tol * ( x_max - x_min );
        xrange[d] = std::max( fabs( x_max ), fabs( x_min ) );
    }
    int Nn    = 3;
    size_t N2 = Nn * N;
    auto x2   = new TYPE[N2 * ndim];
    for ( size_t i = 0; i < N; i++ ) {
        for ( int j = 0; j < Nn; j++ ) {
            for ( int d = 0; d < ndim; d++ )
                x2[d + j * ndim + i * ndim * Nn] =
                    x[d + i * ndim] + 0.1 * dx[d] * ( 2 * rand_double() - 1 );
        }
    }

    // Test the nearest point search
    PROFILE_START( "nearest point search", 1 );
    auto nearest = new unsigned int[N];
    data->find_nearest( (int) N, &x[0], 0, nearest );
    bool pass_nearest = true;
    for ( size_t i = 0; i < N; i++ ) {
        if ( nearest[i] != i )
            pass_nearest = false;
    }
    if ( pass_nearest )
        ut->passes( "Found nearest point" );
    else
        ut->failure( "Found nearest point" );
    PROFILE_STOP( "nearest point search", 1 );

    // Test find_tri
    PROFILE_START( "nearest triangle search", 1 );
    auto index1 = new int[N1];
    auto index2 = new int[N2];
    data->find_tri( (int) N1, x1, 0, index1, false );
    data->find_tri( (int) N2, x2, 0, index2, true );
    bool pass_find_tri = true;
    for ( size_t i = 0; i < N1; i++ ) {
        if ( index1[i] != static_cast<int>( order[i] ) )
            pass_find_tri = false;
    }
    delete[] order;
    delete[] triangles;
    if ( pass_find_tri )
        ut->passes( "Found triangle centroids" );
    else
        ut->failure( "Found triangle centroids" );
    PROFILE_STOP( "nearest triangle search", 1 );

    for ( int p = 0; p < 3; p++ ) {
        if ( N < 10 )
            continue;

        // Evaluate f and the gradient at each point
        std::string problem;
        auto f  = new double[N];         // The function at the verticies
        auto g  = new double[N * ndim];  // The gradient at the verticies
        auto f1 = new double[N1];        // The function at points within the triangles
        auto g1 = new double[N1 * ndim]; // The gradient at points within the triangles
        auto f2 = new double[N2];        // The function at points "near" the verticies
        auto g2 = new double[N2 * ndim]; // The gradient at points "near" the verticies
        double tol_grad, tol_linear, tol_cubic, tol_c_grad, tol_ce1, tol_ce2;
        initialize_problem( p,
                            ndim,
                            N1,
                            x1,
                            xrange,
                            f1,
                            g1,
                            tol_grad,
                            tol_linear,
                            tol_cubic,
                            tol_c_grad,
                            tol_ce1,
                            tol_ce2 );
        initialize_problem( p,
                            ndim,
                            N2,
                            x2,
                            xrange,
                            f2,
                            g2,
                            tol_grad,
                            tol_linear,
                            tol_cubic,
                            tol_c_grad,
                            tol_ce1,
                            tol_ce2 );
        problem = initialize_problem( p,
                                      ndim,
                                      N,
                                      &x[0],
                                      xrange,
                                      f,
                                      g,
                                      tol_grad,
                                      tol_linear,
                                      tol_cubic,
                                      tol_c_grad,
                                      tol_ce1,
                                      tol_ce2 );
        char tmp[256];
        sprintf( tmp,
                 "%s (%i,%i,%s)",
                 problem.c_str(),
                 ndim,
                 static_cast<int>( x.size() ) / ndim,
                 getName<TYPE>() );
        std::string msg( tmp );

        // Test the gradient
        PROFILE_START( "test gradient", 1 );
        auto grad       = new double[ndim * N];
        double error[4] = { 0 };
        for ( int method = 1; method <= 4; method++ ) {
            if ( method == 2 || N < 20 )
                continue;
            data->calc_node_gradient( f, method, grad );
            double err = 0.0;
            for ( size_t i = 0; i < N; i++ ) {
                for ( int d = 0; d < ndim; d++ )
                    err += ( grad[d + i * ndim] - g[d + i * ndim] ) *
                           ( grad[d + i * ndim] - g[d + i * ndim] );
            }
            err               = sqrt( err / N );
            error[method - 1] = err;
            char message[512];
            sprintf( message,
                     "calc_gradient (%i,%s) (%i,%i) (%e)",
                     method,
                     problem.c_str(),
                     ndim,
                     static_cast<int>( x.size() ) / ndim,
                     err );
            if ( err <= tol_grad )
                ut->passes( message );
            else
                ut->failure( message );
        }
        NULL_USE( error );
        PROFILE_STOP( "test gradient", 1 );

        // Test nearest-neighbor interpolation
        PROFILE_START( "test nearest interpolation", 1 );
        auto fi = new double[N];
        data->interp_nearest( f, (int) N, &x[0], nearest, fi );
        bool pass_interp_nearest = true;
        for ( size_t i = 0; i < N; i++ ) {
            if ( !approx_equal( f[i], fi[i] ) )
                pass_interp_nearest = false;
        }
        if ( pass_interp_nearest )
            ut->passes( "nearest-neighbor interpolation: " + msg );
        else
            ut->failure( "nearest-neighbor interpolation: " + msg );
        delete[] fi;
        PROFILE_STOP( "test nearest interpolation", 1 );

        // Allocate variables for interpolation
        auto fi1 = new double[N1];
        auto fi2 = new double[N2];
        auto fi3 = new double[N2];
        auto fi4 = new double[N2];
        auto gi1 = new double[N1 * ndim];
        auto gi2 = new double[N2 * ndim];
        auto gi3 = new double[N2 * ndim];
        auto gi4 = new double[N2 * ndim];
        double f1_norm, f2_norm, f3_norm, f4_norm;
        double g1_norm, g2_norm, g3_norm, g4_norm;
        double f1_err, f2_err, f3_err, f4_err;
        double g1_err, g2_err, g3_err, g4_err;
        auto f_mask         = new bool[N2];
        auto g_mask         = new bool[N2 * ndim];
        bool *f_extrap_mask = nullptr, *g_extrap_mask = nullptr;
        if ( !check_extrap || N < 50 ) {
            f_extrap_mask = f_mask;
            g_extrap_mask = g_mask;
        }

        // Test the linear interpolation
        PROFILE_START( "test linear interpolation", 1 );
        memset( fi1, 0, N1 * sizeof( double ) );
        memset( fi2, 0, N2 * sizeof( double ) );
        memset( fi3, 0, N2 * sizeof( double ) );
        data->interp_linear( f, (int) N1, x1, index1, fi1, nullptr, false );
        data->interp_linear( f, (int) N2, x2, index2, fi2, nullptr, false );
        data->interp_linear( f, (int) N2, x2, index2, fi3, nullptr, true );
        for ( unsigned int i = 0; i < N2; i++ ) {
            f_mask[i] = fi2[i] == fi2[i];
            for ( int d = 0; d < ndim; d++ )
                g_mask[d + i * ndim] = f_mask[i];
        }
        f1_norm = L2norm( N1, f1 );
        f2_norm = L2norm( N2, f2, f_mask );
        f3_norm = L2norm( N2, f2, f_extrap_mask );
        f1_err  = L2errNorm( N1, f1, fi1 );
        f2_err  = L2errNorm( N2, f2, fi2, f_mask );
        f3_err  = L2errNorm( N2, f2, fi3, f_extrap_mask );
        if ( f1_err <= tol_linear * f1_norm && f2_err <= tol_linear * f2_norm &&
             f3_err <= tol_linear * f3_norm )
            ut->passes( "linear interpolation: " + msg );
        else
            ut->failure( "linear interpolation: " + msg );
        PROFILE_STOP( "test linear interpolation", 1 );
        // Test the cubic interpolation (using the exact gradient)
        PROFILE_START( "test cubic interpolation", 1 );
        memset( fi1, 0, N1 * sizeof( double ) );
        memset( fi2, 0, N2 * sizeof( double ) );
        memset( fi3, 0, N2 * sizeof( double ) );
        memset( fi4, 0, N2 * sizeof( double ) );
        memset( gi1, 0, N1 * ndim * sizeof( double ) );
        memset( gi2, 0, N2 * ndim * sizeof( double ) );
        data->interp_cubic( f, g, (int) N1, x1, index1, fi1, gi1, false );
        data->interp_cubic( f, g, (int) N2, x2, index2, fi2, gi2, 0 );
        data->interp_cubic( f, g, (int) N2, x2, index2, fi3, gi3, 1 );
        data->interp_cubic( f, g, (int) N2, x2, index2, fi4, gi4, 2 );
        f1_norm                 = L2norm( N1, f1 );
        f2_norm                 = L2norm( N2, f2, f_mask );
        f3_norm                 = L2norm( N2, f2, f_extrap_mask );
        f4_norm                 = L2norm( N2, f2, f_extrap_mask );
        g1_norm                 = std::max( L2norm( N1 * ndim, g1 ), 1.0 );
        g2_norm                 = std::max( L2norm( N2 * ndim, g2, g_mask ), 1.0 );
        g3_norm                 = std::max( L2norm( N2 * ndim, g2 ), 1.0 );
        g4_norm                 = std::max( L2norm( N2 * ndim, g2 ), 1.0 );
        f1_err                  = L2errNorm( N1, f1, fi1 );
        f2_err                  = L2errNorm( N2, f2, fi2, f_mask );
        f3_err                  = L2errNorm( N2, f2, fi3, f_extrap_mask );
        f4_err                  = L2errNorm( N2, f2, fi4, f_extrap_mask );
        g1_err                  = L2errNorm( N1 * ndim, g1, gi1 );
        g2_err                  = L2errNorm( N2 * ndim, g2, gi2, g_mask );
        g3_err                  = L2errNorm( N2 * ndim, g2, gi3, g_extrap_mask );
        g4_err                  = L2errNorm( N2 * ndim, g2, gi4, g_extrap_mask );
        bool pass_interp_cubic1 = f1_err <= tol_cubic * f1_norm && g1_err <= tol_c_grad * g1_norm;
        bool pass_interp_cubic2 = f2_err <= tol_cubic * f2_norm && g2_err <= tol_c_grad * g2_norm;
        bool pass_interp_cubic3 = f3_err <= tol_linear * f3_norm; // The error in the gradient (and
                                                                  // points) is large for points
                                                                  // outside the domain
        bool pass_interp_cubic4 = f4_err <= tol_cubic * f4_norm && g4_err <= tol_c_grad * g4_norm;
        NULL_USE( g3_norm );
        NULL_USE( g3_err );
        if ( pass_interp_cubic1 && pass_interp_cubic2 && pass_interp_cubic3 &&
             pass_interp_cubic4 ) {
            ut->passes( "cubic interpolation: " + msg );
        } else {
            ut->failure( "cubic interpolation: " + msg );
        }
        PROFILE_STOP( "test cubic interpolation", 1 );

        // Free some memory
        delete[] f;
        delete[] f1;
        delete[] f2;
        delete[] g;
        delete[] g1;
        delete[] g2;
        delete[] grad;
        delete[] fi1;
        delete[] fi2;
        delete[] fi3;
        delete[] fi4;
        delete[] gi1;
        delete[] gi2;
        delete[] gi3;
        delete[] gi4;
        delete[] f_mask;
        delete[] g_mask;
    }
    // Free memory
    delete[] nearest;
    delete[] index1;
    delete[] index2;
    delete[] x1;
    delete[] x2;
}


// This function creates a set of points within box1 that are not contained within box2
// Each box is specified as [xmin,xmax,ymin,ymax,zmin,zmax] and box2 may be NULL
// The grid spacing is given by dx
template<class TYPE>
std::vector<TYPE>
createBoxPoints( int ndim, const int *box1, const int *dx, const int *box2 = nullptr )
{
    int N[5] = { 0, 0, 0, 0, 0 };
    int Nt   = 1;
    for ( int i = 0; i < ndim; i++ ) {
        N[i] = static_cast<int>( floor( ( box1[2 * i + 1] - box1[2 * i + 0] ) / dx[i] ) ) + 1;
        Nt *= N[i];
    }
    std::vector<TYPE> points;
    // Create all points in the box1 and not within box2
    points.reserve( ndim * Nt );
    if ( ndim == 1 ) {
        for ( int i = 0; i < N[0]; i++ ) {
            int x = box1[0] + i * dx[0];
            if ( box2 != nullptr ) {
                if ( x >= box2[0] && x <= box2[1] )
                    continue;
            }
            points.push_back( x );
        }
    } else if ( ndim == 2 ) {
        for ( int i = 0; i < N[0]; i++ ) {
            int x = box1[0] + i * dx[0];
            for ( int j = 0; j < N[1]; j++ ) {
                int y = box1[2] + j * dx[1];
                if ( box2 != nullptr ) {
                    if ( x >= box2[0] && x <= box2[1] && y >= box2[3] && y <= box2[4] )
                        continue;
                }
                points.push_back( x );
                points.push_back( y );
            }
        }
    } else if ( ndim == 3 ) {
        for ( int i = 0; i < N[0]; i++ ) {
            int x = box1[0] + i * dx[0];
            for ( int j = 0; j < N[1]; j++ ) {
                int y = box1[2] + j * dx[1];
                for ( int k = 0; k < N[2]; k++ ) {
                    int z = box1[4] + k * dx[2];
                    if ( box2 != nullptr ) {
                        if ( x >= box2[0] && x <= box2[1] && y >= box2[2] && y <= box2[3] &&
                             z >= box2[4] && z <= box2[5] )
                            continue;
                    }
                    points.push_back( x );
                    points.push_back( y );
                    points.push_back( z );
                }
            }
        }
    }
    return points;
}


template<class TYPE>
std::pair<int, std::vector<TYPE>> createProblem( int problem )
{
    int ndim = 0;
    std::vector<TYPE> points;
    if ( problem == 1 ) {
        // This is a simple problem of a uniform mesh
        ndim       = 3;
        int box[6] = { 0, 4, 0, 4, 0, 4 };
        int dx[3]  = { 1, 1, 1 };
        points     = createBoxPoints<TYPE>( 3, box, dx );
    } else if ( problem == 2 ) {
        // This set was taken from a failed hierarchy that caused an infinite loop
        ndim                   = 3;
        int box1[6]            = { -3, 65, 61, 129, 29, 65 };
        int box2[6]            = { 0, 62, 60, 126, 28, 62 };
        int box3[6]            = { 0, 62, 62, 126, 30, 62 };
        int dx1[3]             = { 4, 4, 4 };
        int dx3[3]             = { 2, 2, 2 };
        std::vector<TYPE> set1 = createBoxPoints<TYPE>( 3, box1, dx1, box2 );
        std::vector<TYPE> set2 = createBoxPoints<TYPE>( 3, box3, dx3 );
        points.reserve( set1.size() + set2.size() );
        for ( auto &elem : set1 )
            points.push_back( elem );
        for ( auto &elem : set2 )
            points.push_back( elem );
    } else if ( problem == 3 ) {
        // This is reduced version of problem 2
        ndim                   = 3;
        int box1[6]            = { -3, 57, -1, 51, -1, 11 };
        int box2[6]            = { 0, 55, -1, 51, -1, 11 };
        int box3[6]            = { 0, 52, 0, 38, 0, 12 };
        int box4[6]            = { 0, 50, 0, 38, 2, 12 };
        int dx1[3]             = { 4, 4, 4 };
        int dx3[3]             = { 2, 2, 4 };
        std::vector<TYPE> set1 = createBoxPoints<TYPE>( 3, box1, dx1, box2 );
        std::vector<TYPE> set2 = createBoxPoints<TYPE>( 3, box3, dx3, box4 );
        points.reserve( set1.size() + set2.size() );
        for ( auto &elem : set1 )
            points.push_back( elem );
        for ( auto &elem : set2 )
            points.push_back( elem );
    } else if ( problem == 4 ) {
        // Invalid triangle neighbor
        ndim                   = 3;
        int box1[6]            = { -3, 15, -1, 43, -1, 11 };
        int box2[6]            = { 0, 15, -1, 43, -1, 11 };
        int box3[6]            = { 0, 12, 36, 44, 0, 12 };
        int box4[6]            = { 0, 10, 36, 44, 2, 12 };
        int dx1[3]             = { 4, 4, 4 };
        int dx3[3]             = { 2, 2, 4 };
        std::vector<TYPE> set1 = createBoxPoints<TYPE>( 3, box1, dx1, box2 );
        std::vector<TYPE> set2 = createBoxPoints<TYPE>( 3, box3, dx3, box4 );
        points.reserve( set1.size() + set2.size() );
        for ( auto &elem : set1 )
            points.push_back( elem );
        for ( auto &elem : set2 )
            points.push_back( elem );
    } else if ( problem == 5 ) {
        // This is a problem that resulted in no valid flips
        ndim                   = 3;
        int box1[6]            = { -3, 55, -1, 55, -1, 15 };
        int box2[6]            = { 0, 52, -2, 52, -2, 12 };
        int box3[6]            = { 0, 52, 0, 52, 0, 12 };
        int dx1[3]             = { 4, 4, 4 };
        int dx3[3]             = { 2, 2, 2 };
        std::vector<TYPE> set1 = createBoxPoints<TYPE>( 3, box1, dx1, box2 );
        std::vector<TYPE> set2 = createBoxPoints<TYPE>( 3, box3, dx3 );
        points.reserve( set1.size() + set2.size() );
        for ( auto &elem : set1 )
            points.push_back( elem );
        for ( auto &elem : set2 )
            points.push_back( elem );
    } else if ( problem == 6 ) {
        // A set of points that caused the nearest-neighbor point search to fail (01/18/13)
        ndim = 3;
        if ( typeid( TYPE ) == typeid( double ) ) {
            double const_points[30] = { 0.228134639925963, 0.107455922361451, 0.865652664047092,
                                        0.184470611307978, 0.838401675869785, 0.668757150119921,
                                        0.369538598042889, 0.123727772818859, 0.640994658126969,
                                        0.499916995230264, 0.448223924105178, 0.071055349040041,
                                        0.234558133130565, 0.459801521022238, 0.439658494159647,
                                        0.053472364370121, 0.256549138105335, 0.295241277565778,
                                        0.375079637861114, 0.623130377070684, 0.385421525726560,
                                        0.965203651821110, 0.545386155470929, 0.355259056526340,
                                        0.034658433895993, 0.917076353750111, 0.761392645215632,
                                        0.350240130632096, 0.309851148106666, 0.654482311175926 };
            points.resize( 30 );
            for ( int i = 0; i < 30; i++ )
                points[i] = const_points[i];
        }
    }
    return std::pair<int, std::vector<TYPE>>( ndim, points );
}


// Test the random number generator
void testRandomNumber( AMP::UnitTest *ut, size_t N )
{
    PROFILE_START( "testRandomNumber" );
    std::vector<double> data( N );
    for ( auto &elem : data )
        elem = rand_double();
    std::sort( data.begin(), data.end() );
    size_t N_dup = 0;
    for ( size_t i = 1; i < data.size(); i++ ) {
        if ( data[i] == data[i - 1] )
            N_dup++;
    }
    if ( data[0] < 0.0 || data[N - 1] > 1.0 )
        ut->failure( "Point outside of [0,1] found" );
    if ( N_dup == 0 ) {
        ut->passes( "No duplicate random numbers" );
    } else {
        char message[100];
        sprintf( message,
                 "Found %i duplicates, expected %e\n",
                 static_cast<int>( N_dup ),
                 N * 2.2204e-16 );
        if ( N_dup < 5 )
            ut->expected_failure( message );
        else
            ut->failure( message );
    }
    PROFILE_STOP( "testRandomNumber" );
}


// Get the convergence
template<class TYPE>
void testConvergence( AMP::UnitTest *ut, int ndim )
{
    // Create two grids with different resolution
    int box[6] = { -1000, 1000, -1000, 1000, -1000, 1000 };
    int dx1[6] = { 100, 100, 100 };
    int dx2[6] = { 50, 50, 50 };
    auto set1  = createBoxPoints<TYPE>( 3, box, dx1 );
    auto set2  = createBoxPoints<TYPE>( 3, box, dx2 );
    // Create the tessellation
    auto data1 = createAndTestDelaunayInterpolation<TYPE>( ut, ndim, set1 );
    auto data2 = createAndTestDelaunayInterpolation<TYPE>( ut, ndim, set2 );
    // Check the convergence for different interpolation method
    for ( int method = 0; method < 3; method++ ) {}
    NULL_USE( data1 );
    NULL_USE( data2 );
}


// Main
int main( int argc, char *argv[] )
{
    AMP::AMPManager::startup( argc, argv );
    AMP::UnitTest ut;
    srand( static_cast<unsigned int>( time( nullptr ) ) );
    PROFILE_ENABLE( 3 ); // 0: code, 1: tests, 3: basic timers, 5: all timers

    // Check that we can create "random" points
    printp( "Running random number test\n" );
    testRandomNumber( &ut, 100000 );

    // Run some basic tests on the point search
    printp( "Running point search tests\n" );
    for ( int d = 1; d <= 5; d++ ) {
        testPointSearch( &ut, d, createRandomPoints<double>( d, 10 ) );
        testPointSearch( &ut, d, createRandomPoints<double>( d, 20 ) );
        testPointSearch( &ut, d, createRandomPoints<double>( d, 50 ) );
        testPointSearch( &ut, d, createRandomPoints<double>( d, 3000 ) );
        testPointSearch( &ut, d, createRandomPoints<int>( d, 100 ) );
    }

    // Test DelaunayInterpolation with random points in 1D, 2D and 3D
    for ( int d = 1; d <= NDIM_MAX; d++ ) {
        printp( "Running interpolation tests with random points %iD (double)\n", d );
        testInterpolation<double>(
            &ut, d, createRandomPoints<double>( d, d + 1 ) ); // minimum # of points
        testInterpolation<double>(
            &ut, d, createRandomPoints<double>( d, 10 ) ); // small # of points
        testInterpolation<double>(
            &ut, d, createRandomPoints<double>( d, 1000 ) ); // medium # of points
        testInterpolation<double>(
            &ut, d, createRandomPoints<double>( d, 10000 ) ); // large # of points
        printp( "Running interpolation tests with random points %iD (int)\n", d );
        testInterpolation<int>(
            &ut, d, createRandomPoints<int>( d, d + 1 ) );                    // minimum # of points
        testInterpolation<int>( &ut, d, createRandomPoints<int>( d, 10 ) );   // small # of points
        testInterpolation<int>( &ut, d, createRandomPoints<int>( d, 1000 ) ); // medium # of points
        // testInterpolation<int>( &ut, d, createRandomPoints<int>(d,10000) );         // large # of
        // points
    }

    // Run some predefined tests
    for ( int i = 1; i <= 6; i++ ) {
        printp( "Running predefined test %i\n", i );
        char tmp[100];
        sprintf( tmp, "Problem %i", i );
        PROFILE_START( tmp );
        std::pair<int, std::vector<int>> data1 = createProblem<int>( i );
        testInterpolation<int>( &ut, data1.first, data1.second );
        std::pair<int, std::vector<double>> data2 = createProblem<double>( i );
        testInterpolation<double>( &ut, data2.first, data2.second );
        PROFILE_STOP( tmp );
    }

    // Run any input file problems
    for ( int i = 1; i < argc; i++ ) {
        printp( "Running input: %s\n", argv[i] );
        std::pair<int, std::vector<double>> data = readPoints( argv[i] );
        int ndim                                 = data.first;
        if ( ndim <= NDIM_MAX )
            testInterpolation<double>( &ut, data.first, data.second, false );
        else
            testPointSearch( &ut, data.first, data.second );
    }

    PROFILE_SAVE( "test_DelaunayInterpolation" );
    ut.report();
    auto N_errors = static_cast<int>( ut.NumFailGlobal() );
    AMP::AMPManager::shutdown();
    return N_errors;
}
