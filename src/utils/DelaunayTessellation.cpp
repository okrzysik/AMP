#include "AMP/utils/DelaunayTessellation.h"
#include "AMP/utils/DelaunayHelpers.h"
#include "AMP/utils/NearestPairSearch.h"
#include "AMP/utils/Utilities.h"

#include "ProfilerApp.h"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <stdexcept>

#define ASSERT AMP_ASSERT
#define ERROR AMP_ERROR

#define DEBUG_CHECK 0 // Flag to enable extra checks (1: some additional cost, 2: very expensive)

#define printp printf
using AMP::Utilities::quicksort;


// Macros to define 3 levels of profilers
#define PROFILE_START_L1( MSG ) PROFILE_START( MSG, 3 )
#define PROFILE_STOP_L1( MSG ) PROFILE_STOP( MSG, 3 )
#define PROFILE_STOP2_L1( MSG ) PROFILE_STOP2( MSG, 3 )
#if 0
#define PROFILE_START_L2( MSG ) PROFILE_START( MSG, 4 )
#define PROFILE_STOP_L2( MSG ) PROFILE_STOP( MSG, 4 )
#define PROFILE_STOP2_L2( MSG ) PROFILE_STOP2( MSG, 4 )
#else
#define PROFILE_START_L2( MSG ) \
    do {                        \
    } while ( 0 )
#define PROFILE_STOP_L2( MSG ) \
    do {                       \
    } while ( 0 )
#define PROFILE_STOP2_L2( MSG ) \
    do {                        \
    } while ( 0 )
#endif
#if 0
#define PROFILE_START_L3( MSG ) PROFILE_START( MSG, 5 )
#define PROFILE_STOP_L3( MSG ) PROFILE_STOP( MSG, 5 )
#define PROFILE_STOP2_L3( MSG ) PROFILE_STOP2( MSG, 5 )
#else
#define PROFILE_START_L3( MSG ) \
    do {                        \
    } while ( 0 )
#define PROFILE_STOP_L3( MSG ) \
    do {                       \
    } while ( 0 )
#define PROFILE_STOP2_L3( MSG ) \
    do {                        \
    } while ( 0 )
#endif


// Choose our numerical format and precision
// Note: for exact integer arithmetic, we need to support at least O(N^D)
//    where N is the index range and D is the dimension.
//    Note: I have not determined the small constant, but is is likely ~4
//    Note: Some routines use higher precisision internally (test_in_circumsphere)
static inline double get_double( const int &x ) { return static_cast<double>( x ); }
static inline double get_double( const int64_t &x ) { return static_cast<double>( x ); }
static inline double get_double( const double &x ) { return static_cast<double>( x ); }
static inline double get_double( const long double &x ) { return static_cast<double>( x ); }
#if 0
#include "samrutils/utilities/extended_int.h"
    typedef extended::int128_t int128_t;
    typedef extended::int256_t int256_t;
    typedef extended::int512_t int512_t;
    static inline double get_double(const int128_t& x) { return x.get_double(); }
    static inline double get_double(const int256_t& x) { return x.get_double(); }
    static inline double get_double(const int512_t& x) { return x.get_double(); }
#elif 0 // USE_BOOST
#include "boost/multiprecision/cpp_int.hpp"
typedef boost::multiprecision::int128_t int128_t;
typedef boost::multiprecision::int256_t int256_t;
typedef boost::multiprecision::int512_t int512_t;
static inline double get_double( const int128_t &x ) { return x.convert_to<double>(); }
static inline double get_double( const int256_t &x ) { return x.convert_to<double>(); }
static inline double get_double( const int512_t &x ) { return x.convert_to<double>(); }
#else
#define DISABLE_EXTENDED
#endif


#ifndef USE_LAPACK
#define USE_LAPACK 0
#endif


#ifdef USE_WINDOWS
#define log2( x ) ( log( x ) / log( 2.0 ) )
#endif

#define MIN( a, b ) ( ( a <= b ) ? a : b )
#define MAX( a, b ) ( ( a >= b ) ? a : b )


static bool conserved_neighbors( const int N1, const int list1[], const int N2, const int list2[] );
bool are_tri_neighbors( const int ndim, const int tri1[], const int tri2[], int *f1, int *f2 );


template<class T>
static inline T *ptr( std::vector<T> x )
{
    return ( x.empty() ) ? nullptr : &x[0];
}


/********************************************************************
 * Compute the dot product of two vectors                            *
 * Note: We are already using an increased precision, and want to    *
 * maintain the maximum degree of accuracy possible.                 *
 ********************************************************************/
template<class TYPE>
inline double dot( int N, const TYPE *x, const TYPE *y );
// Approximate dot product for long double precision
template<>
inline double dot<long double>( int N, const long double *x, const long double *y )
{
    bool sign[8];
    long double z[8];
    for ( int i = 0; i < N; i++ ) {
        z[i]    = x[i] * y[i];
        sign[i] = z[i] >= 0;
    }
    long double ans = z[0];
    z[0]            = 0;
    while ( true ) {
        int index     = -1;
        bool sign_ans = ans >= 0;
        for ( int i = 1; i < N; i++ ) {
            if ( ( sign_ans != sign[i] ) && z[i] != 0 )
                index = i;
        }
        if ( index == -1 )
            break;
        ans += z[index];
        z[index] = 0;
    }
    for ( int i = 1; i < N; i++ ) {
        if ( z[i] != 0 )
            ans += z[i];
    }
    return get_double( ans );
}
// Integer based dot products
template<>
inline double dot<int>( int N, const int *x, const int *y )
{
    int64_t ans( 0 );
    for ( int i = 0; i < N; i++ )
        ans += int64_t( x[i] ) * int64_t( y[i] );
    return get_double( ans );
}
#ifndef DISABLE_EXTENDED
template<>
inline double dot<int64_t>( int N, const int64_t *x, const int64_t *y )
{
    int128_t ans( 0 );
    for ( int i = 0; i < N; i++ )
        ans += int128_t( x[i] ) * int128_t( y[i] );
    return get_double( ans );
}
template<>
inline double dot<int128_t>( int N, const int128_t *x, const int128_t *y )
{
    int256_t ans( 0 );
    for ( int i = 0; i < N; i++ )
        ans += int256_t( x[i] ) * int256_t( y[i] );
    return get_double( ans );
}
template<>
inline double dot<int256_t>( int N, const int256_t *x, const int256_t *y )
{
    int512_t ans( 0 );
    for ( int i = 0; i < N; i++ )
        ans += int512_t( x[i] ) * int512_t( y[i] );
    return get_double( ans );
}
#else
// Exact precision dot product when no extended precision class is availible
inline void split64bit( int64_t x, int64_t &a, int64_t &b )
{
    uint64_t tmp = std::abs( x );
    int64_t s    = ( x < 0 ) ? -1 : 1;
    a            = s * static_cast<int64_t>( tmp >> 31 );
    b            = s * static_cast<int64_t>( tmp & 0x7FFFFFFF );
}
template<>
inline double dot<int64_t>( int N, const int64_t *x, const int64_t *y )
{
    ASSERT( N <= 4 );
    int64_t a, b, c, d, e, f, v62 = 0, v31 = 0, v0 = 0;
    for ( int i = 0; i < N; i++ ) {
        split64bit( x[i], a, b );
        split64bit( y[i], c, d );
        v62 += a * c;
        v31 += a * d + b * c;
        v0 += b * d;
    }
    split64bit( v62, a, b );
    split64bit( v31, c, d );
    split64bit( v0, e, f );
    const long double tmp = 2147483648; // 2^31
    long double ans       = tmp * tmp * tmp * ( a + ( b + c + ( d + e + f / tmp ) / tmp ) / tmp );
    auto ans2             = static_cast<double>( ans );
    return ans2;
}
#endif


/********************************************************************
 * Test if 3 points are co-linear                                    *
 ********************************************************************/
template<int NDIM, class TYPE, class ETYPE>
static inline bool collinear( const TYPE x[3][NDIM], double tol )
{
    double r1[NDIM];
    double r2[NDIM];
    double tmp1 = 0, tmp2 = 0;
    for ( int j = 0; j < NDIM; j++ ) {
        r1[j] = static_cast<double>( x[1][j] - x[0][j] );
        r2[j] = static_cast<double>( x[2][j] - x[0][j] );
        tmp1 += r1[j] * r1[j];
        tmp2 += r2[j] * r2[j];
    }
    tmp1 = 1.0 / sqrt( tmp1 );
    tmp2 = 1.0 / sqrt( tmp2 );
    for ( int j = 0; j < NDIM; j++ ) {
        r1[j] *= tmp1;
        r2[j] *= tmp2;
    }
    double dot = 0.0;
    for ( int j = 0; j < NDIM; j++ ) {
        dot += r1[j] * r2[j];
    }
    bool is_collinear = fabs( 1.0 - fabs( dot ) ) <= tol;
    return is_collinear;
}


/********************************************************************
 * Test if 4 points are co-planar                                    *
 *    |  x1  y1  z1  1 |                                             *
 *    |  x2  y2  z2  1 | = 0                                         *
 *    |  x3  y3  z3  1 |                                             *
 *    |  x4  y4  z4  1 |                                             *
 * Note: for exact math this requires N^D precision                  *
 ********************************************************************/
template<int NDIM, class TYPE, class ETYPE>
static inline bool coplanar( const TYPE x[4][NDIM], TYPE tol )
{
    bool is_coplanar = true;
    if ( NDIM == 3 ) {
        ETYPE det( 0 ), one( 1 ), neg( -1 );
        for ( int d = 0; d < NDIM + 1; d++ ) {
            ETYPE M[9];
            for ( int i1 = 0; i1 < d; i1++ ) {
                for ( int i2 = 0; i2 < NDIM; i2++ ) {
                    M[i1 + i2 * NDIM] = ETYPE( x[i1][i2] );
                }
            }
            for ( int i1 = d + 1; i1 < NDIM + 1; i1++ ) {
                for ( int i2 = 0; i2 < NDIM; i2++ ) {
                    M[i1 - 1 + i2 * NDIM] = ETYPE( x[i1][i2] );
                }
            }
            const ETYPE &sign = ( ( NDIM + d ) % 2 == 0 ) ? one : neg;
            det += sign * DelaunayHelpers<NDIM>::det( M );
        }
        is_coplanar = fabs( get_double( det ) ) <= tol;
    } else {
        ERROR( "Not programmed for dimensions != 3" );
    }
    return is_coplanar;
}


/********************************************************************
 * Increase the storage for tri and tri_nab to hold N_tri triangles  *
 ********************************************************************/
static inline size_t
check_tri_size( size_t N_tri, int NDIM, size_t size_old, int *tri[], int *tri_nab[] )
{
    size_t size_new = size_old;
    while ( ( (size_t) N_tri ) * ( ( size_t )( NDIM + 1 ) ) > size_new )
        size_new *= 2;
    if ( size_new != size_old ) {
        int *tri_old     = *tri;
        int *tri_nab_old = *tri_nab;
        int *tri_new = *tri = new int[size_new];
        int *tri_nab_new = *tri_nab = new int[size_new];
        for ( size_t s = 0; s < size_new; s++ )
            tri_new[s] = -1;
        for ( size_t s = 0; s < size_new; s++ )
            tri_nab_new[s] = -1;
        for ( size_t s = 0; s < size_old; s++ )
            tri_new[s] = tri_old[s];
        for ( size_t s = 0; s < size_old; s++ )
            tri_nab_new[s] = tri_nab_old[s];
        delete[] tri_old;
        delete[] tri_nab_old;
    }
    return size_new;
}


/********************************************************************
 * This is the main function that creates the tessellation           *
 ********************************************************************/
template<int NDIM, class TYPE, class ETYPE>
int DelaunayTessellation::create_tessellation( const int N,
                                               const TYPE x[],
                                               int *tri_out[],
                                               int *tri_nab_out[] )
{
    if ( N < NDIM + 1 )
        throw std::logic_error( "Insufficient number of points" );

    PROFILE_START_L1( "create_tessellation" );

    // Check that no two points match and get the closest pair of points
    std::pair<int, int> index_pair( -1, -1 );
    double r_min = 0;
    {
        // Get the bounding box for the domain
        TYPE xmin[3] = { 0, 0, 0 }, xmax[3] = { 0, 0, 0 };
        for ( int i = 0; i < N; i++ ) {
            for ( int d = 0; d < NDIM; d++ ) {
                xmin[d] = MIN( xmin[d], x[d + i * NDIM] );
                xmax[d] = MAX( xmax[d], x[d + i * NDIM] );
            }
        }
        TYPE domain_size = 0;
        for ( int d = 0; d < NDIM; d++ )
            domain_size += ( xmax[d] - xmin[d] ) * ( xmax[d] - xmin[d] );
        domain_size = static_cast<TYPE>( sqrt( static_cast<double>( domain_size ) ) );

        // First, get the two closest points and check that they are not the same
        PROFILE_START_L1( "create-find_min_dist" );
        index_pair = find_min_dist<NDIM, TYPE>( N, x );
        PROFILE_STOP_L1( "create-find_min_dist" );
        int i1 = index_pair.first;
        int i2 = index_pair.second;
        r_min  = 0;
        for ( int d = 0; d < NDIM; d++ ) {
            auto tmp = static_cast<double>( x[d + NDIM * i1] - x[d + NDIM * i2] );
            r_min += tmp * tmp;
        }
        r_min = sqrt( r_min );
        if ( r_min < 1e-8 * domain_size ) {
            printp( "Duplicate or nearly duplicate points detected (%e)\n", r_min );
            PROFILE_STOP2_L1( "create_tessellation" );
            return -7;
        }
    }
    double TOL_VOL       = 0.0;
    double TOL_COLLINEAR = 1e-12;
    double TOL_COPLANAR  = 0.0;
    if ( !std::numeric_limits<TYPE>::is_integer ) {
        TOL_VOL       = ( NDIM <= 2 ? 1e-6 : 1e-3 ) * pow( r_min, NDIM );
        TOL_VOL       = std::min( TOL_VOL, 0.1 );
        TOL_COLLINEAR = 1e-3;
        TOL_COPLANAR  = 1e-8 * r_min * r_min;
    }

    // Next we need to create a list of the order in which we want to insert the values
    // We will do this by sorting the points by order of their distance from the closest pair
    // Compute the square of the radius
    PROFILE_START_L1( "create-order" );
    auto R2 = new ETYPE[N];
    int i1  = index_pair.first;
    for ( int i = 0; i < N; i++ ) {
        R2[i] = ETYPE( 0 );
        for ( int d = 0; d < NDIM; d++ ) {
            ETYPE tmp( x[d + i * NDIM] - x[d + i1 * NDIM] );
            R2[i] += tmp * tmp;
        }
    }
    // Sort the points, keeping track of the index
    auto I = new int[N];
    for ( int i = 0; i < N; i++ ) {
        I[i] = i;
    }
    quicksort( N, R2, I );
    delete[] R2;
    PROFILE_STOP_L1( "create-order" );

    // Resort the first few points so that the first ndim+1 points are not collinear or coplanar
    int ik    = 2;
    int error = 0;
    for ( int i = 2; i <= NDIM; i++ ) {
        switch ( i ) {
        case 2:
            // Find the first point that is not collinear with the first 2 points in I
            TYPE x2[3][NDIM];
            for ( int j = 0; j < NDIM; j++ ) {
                x2[0][j] = x[j + I[0] * NDIM];
                x2[1][j] = x[j + I[1] * NDIM];
            }
            while ( true ) {
                for ( int j = 0; j < NDIM; j++ )
                    x2[2][j] = x[j + I[ik] * NDIM];
                bool is_collinear = collinear<NDIM, TYPE, ETYPE>( x2, TOL_COLLINEAR );
                if ( !is_collinear ) {
                    break;
                }
                ik++;
                if ( ik >= N ) {
                    // No 3 non-collinear points were found
                    printp( "Error: All points are collinear\n" );
                    error = -2;
                    break;
                }
            }
            break;
        case 3:
            if ( NDIM == 3 ) {
                // Find the first point that is not coplanar with the first 3 points in I
                TYPE x2[4][NDIM];
                for ( int i1 = 0; i1 < 3; i1++ ) {
                    for ( int j = 0; j < NDIM; j++ )
                        x2[i1][j] = x[j + I[i1] * NDIM];
                }
                while ( true ) {
                    for ( int j = 0; j < NDIM; j++ )
                        x2[3][j] = x[j + I[ik] * NDIM];
                    bool is_coplanar = coplanar<NDIM, TYPE, ETYPE>( x2, (TYPE) TOL_COPLANAR );
                    if ( !is_coplanar ) {
                        break;
                    }
                    ik++;
                    if ( ik >= N ) {
                        // No 4 non-coplanar points were found
                        printp( "Error: All points are coplanar\n" );
                        error = -2;
                        break;
                    }
                }
            } else {
                printp( "Error: Co-planar check is not programmed for dimensions other than 3\n" );
                error = -2;
                break;
            }
            break;
        default:
            printp( "Error: Not programmed for this number of dimensions\n" );
            error = -1;
            break;
        }
        if ( error != 0 )
            break;
        // Move the points if necessary
        if ( i != ik ) {
            int tmp = I[ik];
            for ( int j = ik; j > i; j-- )
                I[j] = I[j - 1];
            I[i] = tmp;
        }
    }
    if ( error != 0 ) {
        delete[] I;
        PROFILE_STOP2_L1( "create_tessellation" );
        return error;
    }

    size_t size_tri = ( (size_t) 1 ) << (int) ceil(
                          log2( 2.0 * ( (double) N ) *
                                ( NDIM + 1 ) ) ); // Initial ammount of memory to allocate for tri
    auto tri = new int[size_tri];
    for ( size_t i = 0; i < size_tri; i++ )
        tri[i] = -1;
    size_t N_tri = 1;
    for ( int i = 0; i <= NDIM; i++ )
        tri[i] = I[i];
    TYPE x2[NDIM * ( NDIM + 1 )];
    for ( int i = 0; i <= NDIM; i++ ) {
        for ( int j = 0; j < NDIM; j++ )
            x2[j + i * NDIM] = x[j + NDIM * tri[i]];
    }
    double volume = calc_volume<NDIM, TYPE, ETYPE>( x2 );
    if ( fabs( volume ) <= TOL_VOL ) {
        printp( "Error creating initial triangle\n" );
        delete[] I;
        delete[] tri;
        PROFILE_STOP2_L1( "create_tessellation" );
        return -3;
    } else if ( volume < 0 ) {
        // The volume is negitive, swap the last two indicies
        int tmp       = tri[NDIM - 1];
        tri[NDIM - 1] = tri[NDIM];
        tri[NDIM]     = tmp;
    }

    // Maintain a list of the triangle neighbors
    auto tri_nab = new int[size_tri];
    for ( size_t i = 0; i < size_tri; i++ )
        tri_nab[i] = -1;

    // Maintain a list of the triangle faces on the convex hull
    DelaunayTessellation::FaceList<NDIM, TYPE, ETYPE> face_list( N, x, 0, tri, TOL_VOL );

    // Maintain a list of the unused triangles (those that are all -1, but less than N_tri)
    std::vector<size_t> unused;
    unused.reserve( 512 );

    // Maintain a list of surfaces to check
    std::vector<check_surface_struct> check_surface;
    check_surface.reserve( 256 );

    // Subsequently add each point to the convex hull
    PROFILE_START_L1( "create-add_points" );
    int N_new_max    = 1024;
    auto new_tri     = new int[N_new_max * ( NDIM + 1 )];
    auto new_tri_nab = new int[N_new_max * ( NDIM + 1 )];
    auto neighbor    = new int[N_new_max];
    auto face        = new int[N_new_max];
    auto new_tri_id  = new unsigned int[N_new_max];
    for ( int i = NDIM + 1; i < N; i++ ) {
        PROFILE_START_L2( "create-add_triangles" );
        // Add a point to the convex hull and create the new triangles
        while ( face_list.get_N_face() > N_new_max ) {
            N_new_max *= 2;
            delete[] new_tri;
            new_tri = new int[N_new_max * ( NDIM + 1 )];
            delete[] new_tri_nab;
            new_tri_nab = new int[N_new_max * ( NDIM + 1 )];
            delete[] neighbor;
            neighbor = new int[N_new_max];
            delete[] face;
            face = new int[N_new_max];
            delete[] new_tri_id;
            new_tri_id = new unsigned int[N_new_max];
        }
        int N_tri_new = face_list.add_node(
            I[i], unused, N_tri, new_tri_id, new_tri, new_tri_nab, neighbor, face );
        if ( N_tri_new == 0 ) {
            error = -3;
            printp( "Error:  point is inside or on the convex hull\n" );
        } else if ( N_tri_new == -1 ) {
            error = -3;
            printp( "Error:  Invalid triangle created\n" );
        } else if ( N_tri_new < 0 ) {
            error = -4;
            printp( "Error:  unknown error calling add_node\n" );
        }
        if ( error != 0 ) {
            PROFILE_STOP2_L2( "create-add_triangles" );
            break;
        }
        // Increase the storage for tri and tri_nab if necessary
        size_tri = check_tri_size( N_tri, NDIM, size_tri, &tri, &tri_nab );
        // Add each triangle and update the structures
        for ( int j = 0; j < N_tri_new; j++ ) {
            int index_new = new_tri_id[j];
            for ( int j1 = 0; j1 <= NDIM; j1++ )
                tri[j1 + index_new * ( NDIM + 1 )] = new_tri[j1 + j * ( NDIM + 1 )];
            for ( int j1 = 0; j1 <= NDIM; j1++ )
                tri_nab[j1 + index_new * ( NDIM + 1 )] = new_tri_nab[j1 + j * ( NDIM + 1 )];
            std::swap( tri_nab[face[j] + neighbor[j] * ( NDIM + 1 )], index_new );
            ASSERT( index_new == -1 );
        }
#if DEBUG_CHECK == 2
        bool all_valid = check_current_triangles<NDIM, TYPE, ETYPE>(
            N, x, N_tri, tri, tri_nab, unused, TOL_VOL );
        if ( !all_valid ) {
            error = -4;
            break;
        }
#endif
        // Get a list of the surfaces we need to check for a valid tesselation
        for ( int j = 0; j < N_tri_new; j++ ) {
            int index_new = new_tri_id[j];
            for ( int j1 = 0; j1 <= NDIM; j1++ ) {
                if ( tri_nab[j1 + index_new * ( NDIM + 1 )] != -1 ) {
                    bool found = false;
                    for ( auto &elem : check_surface ) {
                        if ( index_new == elem.t1 && j1 == elem.f1 ) {
                            // The surface is already in the list to check
                            found = true;
                            break;
                        }
                        if ( index_new == elem.t2 && j1 == elem.f2 ) {
                            // The surface is already in the list to check
                            found = true;
                            break;
                        }
                    }
                    if ( !found ) {
                        check_surface_struct tmp;
                        tmp.test = 0x00;
                        tmp.t1   = index_new;
                        tmp.f1   = j1;
                        tmp.t2   = tri_nab[j1 + index_new * ( NDIM + 1 )];
                        int m    = tri_nab[j1 + index_new * ( NDIM + 1 )];
                        for ( int j2 = 0; j2 <= NDIM; j2++ ) {
                            if ( tri_nab[j2 + m * ( NDIM + 1 )] == index_new )
                                tmp.f2 = j2;
                        }
                        check_surface.push_back( tmp );
                    }
                }
            }
        }
        PROFILE_STOP_L2( "create-add_triangles" );
        if ( error != 0 )
            break;

        // Now that we have created a new triangle, perform any edge flips that
        // are necessary to create a valid tesselation
        PROFILE_START_L2( "create-edge_flips" );
        size_t it = 1;
        ASSERT( !check_surface.empty() );
        while ( !check_surface.empty() ) {
            if ( it > std::max<size_t>( 500, N_tri ) ) {
                printp( "Error: infinite loop detected\n" );
                error = -5;
                break;
            }
            // First, lets eliminate all the surfaces that are fine
            PROFILE_START_L3( "create-check_surface" );
            for ( auto &elem : check_surface ) {
                /* Check the surface to see if the triangle pairs need to undergo a flip
                 * The surface is invalid and the triangles need to undergo a flip
                 * if the vertex of either triangle lies within the circumsphere of the other.
                 * Note: if the vertex of one triangle lies within the circumsphere of the
                 * other, then the vertex of the other triangle will lie within the circumsphere
                 * of the current triangle.  Consequently we only need to check one vertex/triangle
                 * pair
                 */
                if ( ( elem.test & 0x01 ) != 0 ) {
                    // We already checked this surface
                    continue;
                }
                for ( int j1 = 0; j1 < NDIM + 1; j1++ ) {
                    int m = tri[j1 + elem.t1 * ( NDIM + 1 )];
                    for ( int j2 = 0; j2 < NDIM; j2++ )
                        x2[j2 + j1 * NDIM] = x[j2 + m * NDIM];
                }
                TYPE v[NDIM];
                int m = tri[elem.f2 + elem.t2 * ( NDIM + 1 )];
                for ( int j = 0; j < NDIM; j++ )
                    v[j] = x[j + m * NDIM];
                int test  = test_in_circumsphere<NDIM, TYPE, ETYPE>( x2, v, TOL_VOL );
                elem.test = ( test != 1 ) ? 0xFF : 0x01;
            }
            // Remove all surfaces that are good
            size_t n = 0;
            for ( size_t k = 0; k < check_surface.size(); k++ ) {
                if ( check_surface[k].test != 0xFF ) {
                    std::swap( check_surface[n], check_surface[k] );
                    n++;
                }
            }
            check_surface.resize( n );
            PROFILE_STOP_L3( "create-check_surface" );
            if ( check_surface.size() == 0 ) {
                // All surfaces are good, we are finished
                break;
            }
            // Find a valid flip
            // The maximum flip currently supported is a 4-4 flip
            int index_old[5], new_tri[( NDIM + 1 ) * 4], new_tri_nab[( NDIM + 1 ) * 4];
            int N_tri_old     = 0;
            int N_tri_new     = 0;
            bool flipped_edge = find_flip<NDIM, TYPE, ETYPE>( x,
                                                              tri,
                                                              tri_nab,
                                                              TOL_VOL,
                                                              check_surface,
                                                              N_tri_old,
                                                              N_tri_new,
                                                              index_old,
                                                              new_tri,
                                                              new_tri_nab );
            if ( !flipped_edge ) {
                // We did not find any valid flips, this is an error
                for ( auto &elem : check_surface )
                    elem.test = 0;
                bool test = find_flip<NDIM, TYPE, ETYPE>( x,
                                                          tri,
                                                          tri_nab,
                                                          TOL_VOL,
                                                          check_surface,
                                                          N_tri_old,
                                                          N_tri_new,
                                                          index_old,
                                                          new_tri,
                                                          new_tri_nab );
                printp( "Error: no valid flips detected\n" );
                if ( test )
                    printp( "   Valid flips were detected if we reset check_surface\n" );
                error = -6;
                break;
            }
            // Check that we conserved the boundary triangles
            PROFILE_START_L3( "create-check_tri_neighbor" );
            int old_tri_nab[( NDIM + 1 ) * 4]; // The maximum flip currently supported is a 4-4 flip
            for ( int j1 = 0; j1 < N_tri_old; j1++ ) {
                for ( int j2 = 0; j2 <= NDIM; j2++ ) {
                    int tmp = tri_nab[j2 + index_old[j1] * ( NDIM + 1 )];
                    for ( int j3 = 0; j3 < N_tri_old; j3++ ) {
                        if ( tmp == index_old[j3] )
                            tmp = -2 - tmp;
                    }
                    old_tri_nab[j2 + j1 * ( NDIM + 1 )] = tmp;
                }
            }
            bool pass = conserved_neighbors(
                N_tri_old * ( NDIM + 1 ), old_tri_nab, N_tri_new * ( NDIM + 1 ), new_tri_nab );
            PROFILE_STOP_L3( "create-check_tri_neighbor" );
            if ( !pass ) {
                printp( "Error: triangle neighbors not conserved\n" );
                error = -4;
                break;
            }
            // Delete the old triangles, add the new ones, and update the structures
            PROFILE_START_L3( "create-update_tri" );
            // First get the indicies where we will store the new triangles
            int index_new[4];
            for ( int j = 0; j < N_tri_new; j++ ) {
                if ( j < N_tri_old ) {
                    // We can just reuse the storage of the old indicies
                    index_new[j] = index_old[j];
                } else if ( !unused.empty() ) {
                    // Use a previously freed place
                    index_new[j] = (int) unused.back();
                    unused.pop_back();
                } else {
                    // We will need to expand N_tri
                    size_tri     = check_tri_size( N_tri + 1, NDIM, size_tri, &tri, &tri_nab );
                    index_new[j] = (int) N_tri;
                    N_tri++;
                }
            }
            // Delete the old triangles
            int N_face_update = 0;
            int old_tri_id[4 * ( NDIM + 1 )], old_face_id[4 * ( NDIM + 1 )],
                new_tri_id[4 * ( NDIM + 1 )], new_face_id[4 * ( NDIM + 1 )];
            for ( int j = N_tri_new; j < N_tri_old; j++ )
                unused.push_back( index_old[j] );
            for ( int j1 = 0; j1 < N_tri_old; j1++ ) {
                for ( int j2 = 0; j2 <= NDIM; j2++ ) {
                    int k  = j2 + index_old[j1] * ( NDIM + 1 );
                    tri[k] = -1;
                    if ( tri_nab[k] == -1 ) {
                        // The given face is on the convex hull
                        old_tri_id[N_face_update]  = index_old[j1];
                        old_face_id[N_face_update] = j2;
                        N_face_update++;
                    } else {
                        // The given face is not on the convex hull
                        /*bool is_old = false;
                        for (int j3=0; j3<N_tri_old; j3++) {
                            if ( tri_nab[k]==index_old[j3] )
                                is_old = true;
                        }*/
                    }
                    tri_nab[k] = -2;
                }
            }
            for ( int j1 = 0; j1 < N_tri_old; j1++ ) {
                size_t j2 = 0;
                while ( j2 < check_surface.size() ) {
                    if ( check_surface[j2].t1 == index_old[j1] ||
                         check_surface[j2].t2 == index_old[j1] ) {
                        std::swap( check_surface[j2], check_surface.back() );
                        check_surface.pop_back();
                    } else {
                        j2++;
                    }
                }
            }
            // Create the new triangles
            for ( int j1 = 0; j1 < N_tri_new; j1++ ) {
                for ( int j2 = 0; j2 <= NDIM; j2++ )
                    tri[j2 + index_new[j1] * ( NDIM + 1 )] = new_tri[j2 + j1 * ( NDIM + 1 )];
            }
            // Update the triangle neighbors and determine the new lists of surfaces to check
            int N_face_update2 = 0;
            for ( int j1 = 0; j1 < N_tri_new; j1++ ) {
                for ( int j2 = 0; j2 <= NDIM; j2++ ) {
                    int k1 = index_new[j1];
                    int k2 = new_tri_nab[j2 + j1 * ( NDIM + 1 )];
                    if ( k2 < -1 ) {
                        // Neighbor triangle is one of the new triangles, we only need to update
                        // tr_nab
                        tri_nab[j2 + k1 * ( NDIM + 1 )] = index_new[-k2 - 2];
                    } else if ( k2 == -1 ) {
                        // Face is on the convex hull, we need to update tri_nab and store the face
                        // for updating face_list
                        tri_nab[j2 + k1 * ( NDIM + 1 )] = -1;
                        new_tri_id[N_face_update2]      = index_new[j1];
                        new_face_id[N_face_update2]     = j2;
                        N_face_update2++;
                    } else {
                        // Neighbor triangle is an existing triangle, we need to update tr_nab for
                        // the new triangle,
                        // the exisiting triangle, and add the face to check_surface
                        tri_nab[j2 + k1 * ( NDIM + 1 )] = k2;
                        int tmp[NDIM];
                        for ( int m = 0; m < j2; m++ )
                            tmp[m] = tri[m + k1 * ( NDIM + 1 )];
                        for ( int m = j2 + 1; m <= NDIM; m++ )
                            tmp[m - 1] = tri[m + k1 * ( NDIM + 1 )];
                        int face = -1;
                        for ( int m1 = 0; m1 <= NDIM; m1++ ) {
                            bool found = false;
                            for ( auto &elem : tmp ) {
                                if ( tri[m1 + k2 * ( NDIM + 1 )] == elem )
                                    found = true;
                            }
                            if ( !found ) {
                                face = m1;
                                break;
                            }
                        }
                        tri_nab[face + k2 * ( NDIM + 1 )] = index_new[j1];
                        check_surface_struct tmp2;
                        tmp2.test = 0x00;
                        tmp2.t1   = k1;
                        tmp2.f1   = j2;
                        tmp2.t2   = k2;
                        tmp2.f2   = face;
                        check_surface.push_back( tmp2 );
                    }
                }
            }
            PROFILE_STOP_L3( "create-update_tri" );
            // Update the faces on the convex hull
            if ( N_face_update != N_face_update2 ) {
                printp( "N_face_update = %i, k = %i\n", N_face_update, N_face_update2 );
                printp( "internal error\n" );
                error = -4;
                break;
            }
            face_list.update_face(
                N_face_update, old_tri_id, old_face_id, new_tri_id, new_face_id, tri );
// Check the current triangles for errors (Only when debug is set, very expensive)
#if DEBUG_CHECK == 2
            all_valid = check_current_triangles<NDIM, TYPE, ETYPE>(
                N, x, N_tri, tri, tri_nab, unused, TOL_VOL );
            if ( !all_valid ) {
                error = -4;
                break;
            }
#endif
            it++;
        }
        PROFILE_STOP_L2( "create-edge_flips" );
        if ( error != 0 )
            break;
// Check the current triangles for errors (Only when debug is set, very expensive)
#if DEBUG_CHECK == 2
        all_valid = check_current_triangles<NDIM, TYPE, ETYPE>(
            N, x, N_tri, tri, tri_nab, unused, TOL_VOL );
        if ( !all_valid ) {
            error = -4;
            break;
        }
#endif
    }
    delete[] new_tri_id;
    delete[] new_tri;
    delete[] new_tri_nab;
    delete[] neighbor;
    delete[] face;
    PROFILE_STOP_L1( "create-add_points" );
    check_surface.clear();
    delete[] I;
    if ( error != 0 ) {
        // Errors occured, delete everything and return
        PROFILE_STOP2_L1( "create_tessellation" );
        delete[] tri_nab;
        delete[] tri;
        return error;
    }

    // Delete any unused triangles
    PROFILE_START_L1( "clean up" );
    while ( !unused.empty() ) {
        // First, see if any unused slots are at the end of the list
        size_t index = unused.size();
        for ( size_t i = 0; i < unused.size(); i++ ) {
            if ( unused[i] == N_tri - 1 ) {
                index = i;
                break;
            }
        }
        if ( index < unused.size() ) {
            unused[index] = unused.back();
            unused.pop_back();
            N_tri--;
            continue;
        }
        // Move the last triangle to fill the last gap
        auto new_tri_id = static_cast<int>( unused.back() ); // The new triangle number
        auto old_tri_id = static_cast<int>( N_tri - 1 );     // The old triangle number
        for ( int j = 0; j <= NDIM; j++ ) {
            size_t k1   = j + new_tri_id * ( NDIM + 1 );
            size_t k2   = j + old_tri_id * ( NDIM + 1 );
            tri[k1]     = tri[k2];
            tri_nab[k1] = tri_nab[k2];
            // Update the neighbors
            if ( tri_nab[k1] != -1 ) {
                bool found = false;
                for ( int j2 = 0; j2 <= NDIM; j2++ ) {
                    if ( tri_nab[j2 + tri_nab[k1] * ( NDIM + 1 )] == old_tri_id ) {
                        tri_nab[j2 + tri_nab[k1] * ( NDIM + 1 )] = new_tri_id;
                        found                                    = true;
                    }
                }
                if ( !found ) {
                    printp( "Error with internal structures (delete any unused triangles)\n" );
                    error = -4;
                }
            }
        }
        // Update the face list
        for ( int j = 0; j <= NDIM; j++ ) {
            size_t k1 = j + new_tri_id * ( NDIM + 1 );
            if ( tri_nab[k1] == -1 )
                face_list.update_face( 1, &old_tri_id, &j, &new_tri_id, &j, tri );
        }
        unused.pop_back();
        N_tri--;
    }
    bool null_tri = false;
    for ( size_t i = 0; i < N_tri * ( NDIM + 1 ); i++ ) {
        if ( tri[i] == -1 )
            null_tri = true;
    }
    if ( null_tri ) {
        // We should have removed all the NULL triangles in the previous step
        printp( "Error with internal structures (NULL tri)\n" );
        error = -4;
    }
    PROFILE_STOP_L1( "clean up" );
    if ( error != 0 ) {
        delete[] tri_nab;
        delete[] tri;
        return error;
    }

    // Remove triangles that are only likely to create problesm
    // clean_triangles<NDIM,TYPE,ETYPE>( N, x, N_tri, tri, tri_nab );

    // Resort the triangle indicies so the smallest index is first (should help with caching)
    // Note: this causes tests to fail (not sure why)
    /*if ( NDIM==2 ) {
        for (size_t i=0; i<N_tri; i++) {
            int *t1 = &tri[i*(NDIM+1)];
            int *n1 = &tri_nab[i*(NDIM+1)];
            int t2[NDIM+1], n2[NDIM+1];
            memcpy(t2,t1,(NDIM+1)*sizeof(int));
            memcpy(n2,n1,(NDIM+1)*sizeof(int));
            int j = 0;
            for (int d=1; d<=NDIM; d++) {
                if ( t1[d] < t1[j] )
                    j = d;
            }
            for (int d=0; d<=NDIM; d++) {
                t1[d] = t2[(j+d)%(NDIM+1)];
                n1[d] = n2[(j+d)%(NDIM+1)];
            }
        }
    }*/

    // Resort the triangles so the are sorted by the first index (improves caching)

    // Check the final set of triangles to make sure they are all valid
    bool all_valid =
        check_current_triangles<NDIM, TYPE, ETYPE>( N, x, N_tri, tri, tri_nab, unused, TOL_VOL );
    if ( !all_valid ) {
        printp( "Final check of triangles failed\n" );
        return -4;
    }

    // Copy the results and delete the structures
    // printf("# of triangles at exit = %i\n",(int)N_tri);
    // printf("# of faces on convex hull at exit = %i\n",face_list.get_N_face());
    *tri_out = new int[( NDIM + 1 ) * N_tri];
    for ( size_t i = 0; i < ( NDIM + 1 ) * N_tri; i++ )
        ( *tri_out )[i] = tri[i];
    delete[] tri;
    if ( tri_nab_out != nullptr ) {
        *tri_nab_out = new int[( NDIM + 1 ) * N_tri];
        for ( size_t i = 0; i < ( NDIM + 1 ) * N_tri; i++ )
            ( *tri_nab_out )[i] = tri_nab[i];
    }
    delete[] tri_nab;

    // Finished
    PROFILE_STOP_L1( "create_tessellation" );
    return (int) N_tri;
}


/********************************************************************
 * Function to remove sliver triangles on the surface                *
 ********************************************************************/
template<int NDIM>
inline double vol_sphere( double r );
template<>
inline double vol_sphere<1>( double r )
{
    return 2.0 * r;
}
template<>
inline double vol_sphere<2>( double r )
{
    return 3.141592653589793 * r * r;
}
template<>
inline double vol_sphere<3>( double r )
{
    return 4.188790204786391 * r * r * r;
}
template<int NDIM, class TYPE, class ETYPE>
void DelaunayTessellation::clean_triangles(
    const int N, const TYPE *x, size_t &N_tri, int *tri, int *tri_nab )
{
    // Get a list of all triangles on the boundary and a figure of merit
    // We will use the ratio of the volume of the circumsphere to the volume of the simplex
    std::vector<std::pair<double, size_t>> index;
    index.reserve( 1000 );
    for ( size_t i = 0; i < N_tri; i++ ) {
        bool on_boundary = false;
        for ( int j = 0; j < NDIM + 1; j++ ) {
            if ( tri_nab[j + i * ( NDIM + 1 )] == -1 )
                on_boundary = true;
        }
        if ( on_boundary ) {
            double x2[NDIM * ( NDIM + 1 )];
            for ( int j = 0; j < NDIM + 1; j++ ) {
                int k = tri[j + i * ( NDIM + 1 )];
                for ( int d = 0; d < NDIM; d++ )
                    x2[d + j * NDIM] = x[d + k * NDIM];
            }
            double vol = calc_volume<NDIM, TYPE, ETYPE>( x2 );
            double R, center[NDIM];
            get_circumsphere<NDIM, TYPE, ETYPE>( x2, R, center );
            double quality = vol_sphere<NDIM, TYPE, ETYPE>( R ) / vol;
            index.emplace_back( quality, i );
        }
    }
    // Remove any triangles that fill more than tol of the sphere
    const double tol = 1e-4;
    for ( size_t i = 0; i < index.size(); i++ ) {
        if ( index[i].first > 1.0 / tol ) {
            std::swap( index[i], index.back() );
            index.pop_back();
        }
    }
    if ( index.empty() )
        return;
    // Sort the triangles based on the quality (worst last)
    std::sort( index.begin(), index.end() );
    // Remove the triangles
    while ( !index.empty() ) {
        size_t i = index.back().second;
        index.pop_back();
        // Check that all nodes are shared by at least 1 other triangle
        // we cannot remove any triangles that wil orphan a node
        bool remove = true;
        for ( int j1 = 0; j1 < NDIM + 1; j1++ ) {
            int n        = tri[j1 + i * ( NDIM + 1 )];
            bool found_n = false;
            for ( int j2 = 0; j2 < NDIM + 1; j2++ ) {
                int k = tri_nab[j2 + i * ( NDIM + 1 )];
                if ( k == -1 )
                    continue;
                for ( int j3 = 0; j3 < NDIM + 1; j3++ ) {
                    if ( tri[j3 + k * ( NDIM + 1 )] == n )
                        found_n = true;
                }
            }
            if ( !found_n )
                remove = false;
        }
        // Remove the triangle
        if ( remove ) {
            size_t i2 = N_tri - 1;
            swap_triangles<NDIM, TYPE, ETYPE>( N_tri, i, i2, tri, tri_nab );
            for ( int j = 0; j < NDIM + 1; j++ ) {
                int k                          = tri_nab[j + i2 * ( NDIM + 1 )];
                tri[j + i2 * ( NDIM + 1 )]     = -1;
                tri_nab[j + i2 * ( NDIM + 1 )] = -1;
                for ( int j2 = 0; j2 < NDIM + 1; j2++ ) {
                    if ( tri_nab[j2 + k * ( NDIM + 1 )] == static_cast<int>( i2 ) )
                        tri_nab[j2 + k * ( NDIM + 1 )] = -1;
                }
            }
            N_tri--;
        }
    }
}


/********************************************************************
 * Function to swap two triangle indicies                            *
 ********************************************************************/
template<int NDIM>
void DelaunayTessellation::swap_triangles( size_t N_tri, int i1, int i2, int *tri, int *tri_nab )
{
    // First swap the triangle data
    for ( int j = 0; j < NDIM + 1; j++ ) {
        std::swap( tri[j + i1 * ( NDIM + 1 )], tri[j + i2 * ( NDIM + 1 )] );
        std::swap( tri_nab[j + i1 * ( NDIM + 1 )], tri_nab[j + i2 * ( NDIM + 1 )] );
    }
    // Get a unique list of all neighbors of either triangle (and the triangles themselves)
    int neighbors[2 * NDIM + 4];
    for ( int j = 0; j < NDIM + 1; j++ ) {
        neighbors[2 * j + 0] = tri_nab[j + i1 * ( NDIM + 1 )];
        neighbors[2 * j + 1] = tri_nab[j + i2 * ( NDIM + 1 )];
    }
    neighbors[2 * NDIM + 2] = i1;
    neighbors[2 * NDIM + 3] = i2;
    int N                   = 2 * NDIM + 4;
    int i                   = 0;
    while ( i < N ) {
        bool keep = neighbors[i] != -1;
        for ( int j = 0; j < i; j++ )
            keep = keep && neighbors[i] != neighbors[j];
        if ( !keep ) {
            std::swap( neighbors[i], neighbors[N - 1] );
            N--;
        } else {
            i++;
        }
    }
    // For each triangle swap any values for i1 and i2
    for ( i = 0; i < N; i++ ) {
        int k = neighbors[i];
        for ( int j = 0; j < NDIM + 1; j++ ) {
            if ( tri_nab[j + k * ( NDIM + 1 )] == i1 ) {
                tri_nab[j + k * ( NDIM + 1 )] = i2;
            } else if ( tri_nab[j + k * ( NDIM + 1 )] == i2 ) {
                tri_nab[j + k * ( NDIM + 1 )] = i1;
            }
        }
    }
}


/********************************************************************
 * This function calculates the volume of a N-dimensional simplex    *
 *         1  |  x1-x4   x2-x4   x3-x4  |                            *
 *    V = --  |  y1-y4   y2-y4   y3-y4  |   (3D)                     *
 *        n!  |  z1-z4   z2-z4   z3-z4  |                            *
 * Note: the sign of the volume depends on the order of the points.  *
 *   It will be positive for points stored in a clockwise mannor     *
 * This version is optimized for performance.                        *
 * Note:  If the volume is zero, then the simplex is invalid         *
 *   Eg. a line in 2D or a plane in 3D.                              *
 * Note: exact math requires N^D precision                           *
 ********************************************************************/
template<int NDIM, class TYPE, class ETYPE>
double DelaunayTessellation::calc_volume( const TYPE x[] )
{
    if ( NDIM == 1 ) {
        return static_cast<double>( x[1] - x[0] );
    }
    ETYPE M[NDIM * NDIM];
    for ( int i = 0; i < NDIM * NDIM; i++ )
        M[i] = ETYPE( x[i] );
    for ( int i = 0; i < NDIM; i++ ) {
        ETYPE tmp( x[i + NDIM * NDIM] );
        for ( int j = 0; j < NDIM; j++ )
            M[i + j * NDIM] -= tmp;
    }
    double Volume = get_double( DelaunayHelpers<NDIM>::det( M ) );
    for ( int i = 2; i <= NDIM; i++ )
        Volume /= i;
    return Volume;
}


/************************************************************************
 * This function tests if a point is inside the circumsphere of an       *
 *    nd-simplex.                                                        *
 * For performance, I assume the points are ordered properly such that   *
 * the volume of the simplex (as calculated by calc_volume) is positive. *
 *                                                                       *
 * The point is inside the circumsphere if the determinant is positive   *
 * for points stored in a clockwise manner.  If the order is not known,  *
 * we can compare to a point we know is inside the cicumsphere.          *
 *    |  x1-xi   y1-yi   z1-zi   (x1-xi)^2+(y1-yi)^2+(z1-yi)^2  |        *
 *    |  x2-xi   y2-yi   z2-zi   (x2-xi)^2+(y2-yi)^2+(z2-yi)^2  |        *
 *    |  x3-xi   y3-yi   z3-zi   (x3-xi)^2+(y3-yi)^2+(z3-yi)^2  |        *
 *    |  x4-xi   y4-yi   z4-zi   (x4-xi)^2+(y4-yi)^2+(z4-yi)^2  |        *
 * det(A) == 0:  We are on the circumsphere                              *
 * det(A) > 0:   We are inside the circumsphere                          *
 * det(A) < 0:   We are outside the circumsphere                         *
 *                                                                       *
 * Note: this implimentation requires N^D precision                      *
 ************************************************************************/
template<class TYPE>
static inline int test_in_circumsphere_1D( const TYPE x[], const TYPE xi[], const double TOL_VOL )
{
    if ( fabs( static_cast<double>( *xi - x[0] ) ) <= TOL_VOL ||
         fabs( static_cast<double>( *xi - x[1] ) ) <= TOL_VOL ) {
        return 0;
    }
    if ( ( *xi < x[0] && *xi < x[1] ) || ( *xi > x[0] && *xi > x[1] ) )
        return -1;
    return 1;
}
template<int NDIM, class TYPE, class ETYPE>
int DelaunayTessellation::test_in_circumsphere( const TYPE x[],
                                                const TYPE xi[],
                                                const double TOL_VOL )
{
    if ( NDIM == 1 ) {
        return test_in_circumsphere_1D( x, xi, TOL_VOL );
    }
    // Solve the sub-determinants (requires N^NDIM precision)
    double R2 = 0.0;
    const ETYPE one( 1 ), neg( -1 );
    ETYPE det2[NDIM + 1], R[NDIM + 1];
    for ( int d = 0; d <= NDIM; d++ ) {
        ETYPE A2[NDIM * NDIM];
        ETYPE sum( 0 );
        for ( int j = 0; j < NDIM; j++ ) {
            ETYPE tmp( x[j + d * NDIM] - xi[j] );
            sum += tmp * tmp;
            for ( int i = 0; i < d; i++ )
                A2[i + j * NDIM] = ETYPE( x[j + i * NDIM] - xi[j] );
            for ( int i = d + 1; i <= NDIM; i++ )
                A2[i - 1 + j * NDIM] = ETYPE( x[j + i * NDIM] - xi[j] );
        }
        R[d] = sum;
        R2 += get_double( R[d] );
        const ETYPE &sign = ( ( NDIM + d ) % 2 == 0 ) ? one : neg;
        det2[d]           = sign * DelaunayHelpers<NDIM>::det( A2 );
    }
    // Compute the determinate (requires N^(NDIM+2) precision, used internally in dot)
    double det_A = dot<ETYPE>( NDIM + 1, det2, R );
    if ( fabs( det_A ) <= 0.1 * R2 * TOL_VOL ) {
        // We are on the circumsphere
        return 0;
    } else if ( det_A > 0 ) {
        // We inside the circumsphere
        return 1;
    } else {
        // We outside the circumsphere
        return -1;
    }
}


/************************************************************************
 * This function computes the circumsphere containing the simplex        *
 *     http://mathworld.wolfram.com/Circumsphere.html                    *
 * We have 4 points that define the circumsphere (x1, x2, x3, x4).       *
 * We will work in a reduced coordinate system with x1 at the center     *
 *    so that we can reduce the size of the determinant and the number   *
 *    of significant digits.                                             *
 *                                                                       *
 *     | x2-x1  y2-y1  z2-z1 |                                           *
 * a = | x3-x1  y3-y1  z3-z1 |                                           *
 *     | x4-x1  y4-y1  z4-z1 |                                           *
 *                                                                       *
 *      | (x2-x1)^2+(y2-y1)^2+(z2-z1)^2  y2  z2 |                        *
 * Dx = | (x3-x1)^2+(y3-y1)^2+(z3-z1)^2  y3  z3 |                        *
 *      | (x4-x1)^2+(y4-y1)^2+(z4-z1)^2  y4  z4 |                        *
 *                                                                       *
 *        | (x2-x1)^2+(y2-y1)^2+(z2-z1)^2  x2  z2 |                      *
 * Dy = - | (x3-x1)^2+(y3-y1)^2+(z3-z1)^2  x3  z3 |                      *
 *        | (x4-x1)^2+(y4-y1)^2+(z4-z1)^2  x4  z4 |                      *
 *                                                                       *
 *      | (x2-x1)^2+(y2-y1)^2+(z2-z1)^2  x2  y2 |                        *
 * Dz = | (x3-x1)^2+(y3-y1)^2+(z3-z1)^2  x3  y3 |                        *
 *      | (x4-x1)^2+(y4-y1)^2+(z4-z1)^2  x4  y4 |                        *
 *                                                                       *
 * c = 0                                                                 *
 *                                                                       *
 * Note: This version does not use exact math                            *
 ************************************************************************/
template<class TYPE>
static inline void get_circumsphere_1D( const TYPE x[], double &R, double *center )
{
    center[0] = 0.5 * static_cast<double>( x[0] + x[1] );
    R         = 0.5 * fabs( static_cast<double>( x[0] - x[1] ) );
}
template<int NDIM, class TYPE>
void DelaunayTessellation::get_circumsphere( const TYPE x0[], double &R, double *center )
{
    if ( NDIM == 1 ) {
        get_circumsphere_1D<TYPE>( x0, R, center );
        return;
    }
    long double x[NDIM * NDIM];
    for ( int i = 0; i < NDIM; i++ ) {
        for ( int j = 0; j < NDIM; j++ )
            x[j + i * NDIM] = get_double( x0[j + ( i + 1 ) * NDIM] - x0[j] );
    }
    long double A[NDIM * NDIM], D[NDIM][NDIM * NDIM];
    for ( int i = 0; i < NDIM; i++ ) {
        long double tmp( 0 );
        for ( int j = 0; j < NDIM; j++ ) {
            long double x2 = get_double( x[j + i * NDIM] );
            tmp += x2 * x2;
            A[i + j * NDIM] = x2;
            for ( int k = j + 1; k < NDIM; k++ )
                D[k][i + ( j + 1 ) * NDIM] = x2;
            for ( int k = 0; k < j; k++ )
                D[k][i + j * NDIM] = x2;
        }
        for ( auto &elem : D )
            elem[i] = tmp;
    }
    long double a = get_double( DelaunayHelpers<NDIM>::det( A ) );
    R             = 0.0;
    for ( int i = 0; i < NDIM; i++ ) {
        long double d =
            ( ( i % 2 == 0 ) ? 1 : -1 ) * get_double( DelaunayHelpers<NDIM>::det( D[i] ) );
        center[i] = static_cast<double>( d / ( 2 * a ) + static_cast<long double>( x0[i] ) );
        R += d * d;
    }
    R = sqrt( R ) / fabs( static_cast<double>( 2 * a ) );
}


/********************************************************************
 * This function tests if the edge flip is valid.                    *
 * To be valid, the line between the two points that are not on the  *
 * surface and the intersection of the face must lie within the face *
 * The point of intersection lies within the face if the Barycentric *
 * coordinates for all but the given face are positive               *
 ********************************************************************/
template<class TYPE>
inline double getFlipTOL();
template<>
inline double getFlipTOL<int>()
{
    return 0;
}
template<>
inline double getFlipTOL<int64_t>()
{
    return 0;
}
template<>
inline double getFlipTOL<double>()
{
    return 1e-12;
}
template<>
inline double getFlipTOL<long double>()
{
    return 1e-12;
}
template<int NDIM, class TYPE, class ETYPE>
bool DelaunayTessellation::test_flip_valid( const TYPE x[], const int i, const TYPE xi[] )
{
    const double TOL = getFlipTOL<TYPE>();
    double L[NDIM + 1];
    compute_Barycentric<NDIM, TYPE, ETYPE>( x, xi, L );
    bool is_valid = true;
    for ( int j = 0; j <= NDIM; j++ )
        is_valid = is_valid && ( j == i || L[j] >= -TOL );
    return is_valid;
}


/****************************************************************
 * Function to compute the Barycentric coordinates               *
 * Note: we use exact math until we perform the normalization    *
 *    The exact math component requires N^(D-1) precision        *
 ****************************************************************/
template<int NDIM, class TYPE, class ETYPE>
void DelaunayTessellation::compute_Barycentric( const TYPE *x, const TYPE *xi, double *L )
{
    // Compute the barycentric coordinates T*L=r-r0
    // http://en.wikipedia.org/wiki/Barycentric_coordinate_system_(mathematics)
    ETYPE T[NDIM * NDIM];
    for ( int i = 0; i < NDIM; i++ ) {
        for ( int j = 0; j < NDIM; j++ ) {
            T[j + i * NDIM] = ETYPE( x[j + i * NDIM] - x[j + NDIM * NDIM] );
        }
    }
    ETYPE r[NDIM];
    for ( int i = 0; i < NDIM; i++ )
        r[i] = ETYPE( xi[i] - x[i + NDIM * NDIM] );
    ETYPE L2[NDIM + 1], det( 0 );
    DelaunayHelpers<NDIM>::solve_system( T, r, L2, det );
    L2[NDIM] = det;
    for ( int i = 0; i < NDIM; i++ )
        L2[NDIM] -= L2[i];
    // Perform the normalization (will require inexact math)
    double scale = 1.0 / get_double( det );
    for ( int i = 0; i < NDIM + 1; i++ )
        L[i] = get_double( L2[i] ) * scale;
}


/************************************************************************
 * Function to find a valid flip                                         *
 * Note: the flip must also calculate the triangle neighbors.            *
 *   While this requires more effort when programming flips, the         *
 *   of the flips allow them to potentially calculate the                *
 *   neighbors much faster than a more generic method.                   *
 *      new_tri_nab >= 0   - Triangle neighbor is an existing triangle   *
 *                           (with the given index)                      *
 *      new_tri_nab -(i+1) - Triangle neighbor is the ith new triangle   *
 *      new_tri_nab ==-1   - Triangle face is on the convex hull         *
 ************************************************************************/
template<class TYPE, class ETYPE>
inline bool DelaunayTessellation::find_flip_2D(
    const TYPE *x,
    const int *tri,
    const int *tri_nab,
    const double TOL_VOL,
    std::vector<DelaunayTessellation::check_surface_struct> &check_surface,
    int &N_tri_old,
    int &N_tri_new,
    int *index_old,
    int *new_tri,
    int *new_tri_nab )
{
    PROFILE_START_L3( "find_flip<2>" );
    // In 2D we only have one type of flip (2-2 flip)
    bool found = false;
    for ( auto &elem : check_surface ) {
        if ( ( elem.test & 0x02 ) != 0 ) {
            // We already tried this flip
            continue;
        }
        int s1 = elem.f1;
        int s2 = elem.f2;
        int t1 = elem.t1;
        int t2 = elem.t2;
#if DEBUG_CHECK >= 1
        TYPE x1[6], x2[6], xi1[2], xi2[2];
        for ( int j = 0; j < 3; j++ ) {
            int m1        = tri[j + t1 * 3];
            int m2        = tri[j + t2 * 3];
            x1[0 + j * 2] = x[0 + m1 * 2];
            x1[1 + j * 2] = x[1 + m1 * 2];
            x2[0 + j * 2] = x[0 + m2 * 2];
            x2[1 + j * 2] = x[1 + m2 * 2];
        }
        xi1[0]    = x2[0 + 2 * s2];
        xi1[1]    = x2[1 + 2 * s2];
        xi2[0]    = x1[0 + 2 * s1];
        xi2[1]    = x1[1 + 2 * s1];
        int test1 = test_in_circumsphere<2, TYPE, ETYPE>( x1, xi1, TOL_VOL );
        int test2 = test_in_circumsphere<2, TYPE, ETYPE>( x2, xi2, TOL_VOL );
        if ( test1 != 1 || test2 != 1 )
            throw std::logic_error( "Internal error" );
#endif
        bool flipped_edge = flip_2D<TYPE, ETYPE>(
            x, tri, tri_nab, t1, s1, t2, s2, index_old, new_tri, new_tri_nab, TOL_VOL );
        if ( flipped_edge ) {
            N_tri_old = 2;
            N_tri_new = 2;
            found     = true;
            break;
        } else {
            // The 2-2 flip is not valid, set the second bit to indicate we checked this
            elem.test |= 0x02;
            printp( "Warning: necessary flips in 2D should always be valid\n" );
        }
    }
    PROFILE_STOP_L3( "find_flip<2>" );
    return found;
}
template<class TYPE, class ETYPE>
inline bool DelaunayTessellation::find_flip_3D(
    const TYPE *x,
    const int *tri,
    const int *tri_nab,
    const double TOL_VOL,
    std::vector<DelaunayTessellation::check_surface_struct> &check_surface,
    int &N_tri_old,
    int &N_tri_new,
    int *index_old,
    int *new_tri,
    int *new_tri_nab )
{
    PROFILE_START_L3( "find_flip<3>" );
    bool found = false;
    // In 3D there are several possible types of flips
    // First let's see if there are any valid 2-2 flips
    if ( !found ) {
        for ( auto &elem : check_surface ) {
            if ( ( elem.test & 0x02 ) != 0 ) {
                // We already tried this flip
                continue;
            }
            int s1            = elem.f1;
            int s2            = elem.f2;
            int t1            = elem.t1;
            int t2            = elem.t2;
            bool flipped_edge = flip_3D_22<TYPE, ETYPE>(
                x, tri, tri_nab, t1, s1, t2, s2, index_old, new_tri, new_tri_nab, TOL_VOL );
            if ( flipped_edge ) {
                N_tri_old = 2;
                N_tri_new = 2;
                found     = true;
                break;
            } else {
                // The 2-2 flip is not valid, set the second bit to indicate we checked this
                elem.test |= 0x02;
            }
        }
    }
    // If we did not find a valid flip, check for 2-3 flips
    if ( !found ) {
        for ( auto &elem : check_surface ) {
            if ( ( elem.test & 0x08 ) != 0 ) {
                // We already tried this flip
                continue;
            }
            int s1            = elem.f1;
            int s2            = elem.f2;
            int t1            = elem.t1;
            int t2            = elem.t2;
            bool flipped_edge = flip_3D_23<TYPE, ETYPE>(
                x, tri, tri_nab, t1, s1, t2, s2, index_old, new_tri, new_tri_nab, TOL_VOL );
            if ( flipped_edge ) {
                N_tri_old = 2;
                N_tri_new = 3;
                found     = true;
                break;
            } else {
                // The 2-3 flip is not valid, set the fourth bit to indicate we checked this
                elem.test |= 0x08;
            }
        }
    }
    // If we did not find a valid flip, check for 3-2 flips
    // Note: we must recheck this if we modified a neighboring triangle
    if ( !found ) {
        for ( auto &elem : check_surface ) {
            int s1            = elem.f1;
            int s2            = elem.f2;
            int t1            = elem.t1;
            int t2            = elem.t2;
            bool flipped_edge = flip_3D_32<TYPE, ETYPE>(
                x, tri, tri_nab, t1, s1, t2, s2, index_old, new_tri, new_tri_nab, TOL_VOL );
            if ( flipped_edge ) {
                N_tri_old = 3;
                N_tri_new = 2;
                found     = true;
                break;
            } else {
                // The 3-2 flip is not valid, set the third bit to indicate we checked this
                elem.test |= 0x04;
            }
        }
    }
    // If we did not find a valid flip, check for 4-4 flips
    // Note: we must recheck this if we modified a neighboring triangle
    if ( !found ) {
        for ( auto &elem : check_surface ) {
            int s1            = elem.f1;
            int s2            = elem.f2;
            int t1            = elem.t1;
            int t2            = elem.t2;
            bool flipped_edge = flip_3D_44<TYPE, ETYPE>(
                x, tri, tri_nab, t1, s1, t2, s2, index_old, new_tri, new_tri_nab, TOL_VOL );
            if ( flipped_edge ) {
                N_tri_old = 4;
                N_tri_new = 4;
                found     = true;
                break;
            } else {
                // The 4-4 flip is not valid, set the fifth bit to indicate we checked this
                elem.test |= 0x10;
            }
        }
    }
    PROFILE_STOP_L3( "find_flip<3>" );
    return found;
}
template<int NDIM, class TYPE, class ETYPE>
bool DelaunayTessellation::find_flip(
    const TYPE *x,
    const int *tri,
    const int *tri_nab,
    const double TOL_VOL,
    std::vector<DelaunayTessellation::check_surface_struct> &check_surface,
    int &N_tri_old,
    int &N_tri_new,
    int *index_old,
    int *new_tri,
    int *new_tri_nab )
{
    bool valid = false;
    if ( NDIM == 2 ) {
        valid = find_flip_2D<TYPE, ETYPE>( x,
                                           tri,
                                           tri_nab,
                                           TOL_VOL,
                                           check_surface,
                                           N_tri_old,
                                           N_tri_new,
                                           index_old,
                                           new_tri,
                                           new_tri_nab );
    } else if ( NDIM == 3 ) {
        valid = find_flip_3D<TYPE, ETYPE>( x,
                                           tri,
                                           tri_nab,
                                           TOL_VOL,
                                           check_surface,
                                           N_tri_old,
                                           N_tri_new,
                                           index_old,
                                           new_tri,
                                           new_tri_nab );
    }
    return valid;
}


/************************************************************************
 * This function performs a flip in 2D                                   *
 * Note: all indicies are hard-coded for ndim=2 and some loops have been *
 * unrolled to simplify the code and improve performance                 *
 ************************************************************************/
template<class TYPE, class ETYPE>
bool DelaunayTessellation::flip_2D( const TYPE x[],
                                    const int tri[],
                                    const int tri_nab[],
                                    const int t1,
                                    const int s1,
                                    const int t2,
                                    const int s2,
                                    int *index_old,
                                    int *new_tri,
                                    int *new_tri_nab,
                                    const double TOL_VOL )
{
    // Check if the flip is valid (it should always be valid in 2D if it is necessary)
    TYPE x2[6], xi[2];
    for ( int i = 0; i < 3; i++ ) {
        int k         = tri[i + t1 * 3];
        x2[0 + i * 2] = x[0 + k * 2];
        x2[1 + i * 2] = x[1 + k * 2];
    }
    int k        = tri[s2 + t2 * 3];
    xi[0]        = x[0 + 2 * k];
    xi[1]        = x[1 + 2 * k];
    bool isvalid = test_flip_valid<2, TYPE, ETYPE>( x2, s1, xi );
    if ( !isvalid ) {
        // We need to do an edge flip, and the edge flip is not valid
        // This should not actually occur
        return false;
    }
    // Create the new triangles that are formed by the edge flips and the neighbor list
    index_old[0] = t1;
    index_old[1] = t2;
    int is[2], tn[4];
    if ( s1 == 0 ) {
        is[0] = tri[1 + 3 * t1];
        is[1] = tri[2 + 3 * t1];
        tn[0] = tri_nab[1 + 3 * t1];
        tn[1] = tri_nab[2 + 3 * t1];
    } else if ( s1 == 1 ) {
        is[0] = tri[0 + 3 * t1];
        is[1] = tri[2 + 3 * t1];
        tn[0] = tri_nab[0 + 3 * t1];
        tn[1] = tri_nab[2 + 3 * t1];
    } else {
        is[0] = tri[0 + 3 * t1];
        is[1] = tri[1 + 3 * t1];
        tn[0] = tri_nab[0 + 3 * t1];
        tn[1] = tri_nab[1 + 3 * t1];
    }
    if ( is[0] == tri[0 + 3 * t2] )
        tn[2] = tri_nab[0 + 3 * t2];
    else if ( is[0] == tri[1 + 3 * t2] )
        tn[2] = tri_nab[1 + 3 * t2];
    else
        tn[2] = tri_nab[2 + 3 * t2];
    if ( is[1] == tri[0 + 3 * t2] )
        tn[3] = tri_nab[0 + 3 * t2];
    else if ( is[1] == tri[1 + 3 * t2] )
        tn[3] = tri_nab[1 + 3 * t2];
    else
        tn[3] = tri_nab[2 + 3 * t2];
    int in1        = tri[s1 + 3 * t1];
    int in2        = tri[s2 + 3 * t2];
    new_tri[0]     = in1;
    new_tri[1]     = in2;
    new_tri[2]     = is[0];
    new_tri[3]     = in1;
    new_tri[4]     = in2;
    new_tri[5]     = is[1];
    new_tri_nab[0] = tn[3];
    new_tri_nab[1] = tn[1];
    new_tri_nab[2] = -3;
    new_tri_nab[3] = tn[2];
    new_tri_nab[4] = tn[0];
    new_tri_nab[5] = -2;
    // Check that the new triagles are valid (and have the proper ordering)
    isvalid = true;
    for ( int it = 0; it < 2; it++ ) {
        for ( int i = 0; i < 3; i++ ) {
            int k         = new_tri[i + it * 3];
            x2[0 + i * 2] = x[0 + k * 2];
            x2[1 + i * 2] = x[1 + k * 2];
        }
        double volume = calc_volume<2, TYPE, ETYPE>( x2 );
        if ( fabs( volume ) <= TOL_VOL ) {
            // The triangle is invalid (collinear)
            isvalid = false;
        } else if ( volume < 0 ) {
            // The ordering of the points is invalid, swap the last two points (we want the volume
            // to be positive)
            int tmp                 = new_tri[1 + it * 3];
            new_tri[1 + it * 3]     = new_tri[2 + it * 3];
            new_tri[2 + it * 3]     = tmp;
            tmp                     = new_tri_nab[1 + it * 3];
            new_tri_nab[1 + it * 3] = new_tri_nab[2 + it * 3];
            new_tri_nab[2 + it * 3] = tmp;
        }
    }
    return isvalid;
}


/************************************************************************
 * This function performs a 2-2 flip in 3D                               *
 * Note: all indicies are hard-coded for ndim=3 and some loops may be    *
 * unrolled to simplify the code or improve performance                  *
 ************************************************************************/
template<class TYPE, class ETYPE>
bool DelaunayTessellation::flip_3D_22( const TYPE x[],
                                       const int tri[],
                                       const int tri_nab[],
                                       const int t1,
                                       const int s1,
                                       const int t2,
                                       const int s2,
                                       int *index_old,
                                       int *new_tri,
                                       int *new_tri_nab,
                                       const double TOL_VOL )
{
    // The 2-2 fip is only valid when 4 of the verticies lie on a plane on the convex hull
    // This is likely for structured grids
    for ( int i1 = 0; i1 < 4; i1++ ) {
        if ( tri_nab[i1 + 4 * t1] != -1 ) {
            // The face does not lie on the convex hull
            continue;
        }
        for ( int i2 = 0; i2 < 4; i2++ ) {
            if ( tri_nab[i2 + 4 * t2] != -1 ) {
                // The face does not lie on the convex hull
                continue;
            }
            if ( tri[i1 + 4 * t1] != tri[i2 + 4 * t2] ) {
                // The 5th vertex does not match for the two triangles
                continue;
            }
            // Get the list of the 4 verticies on the convex hull
            int v1, v2, v3, v4;
            v1 = tri[s1 + 4 * t1];
            v2 = tri[s2 + 4 * t2];
            v3 = -1;
            v4 = -1;
            for ( int j = 0; j < 4; j++ ) {
                if ( j == i1 || j == s1 )
                    continue;
                if ( v3 == -1 )
                    v3 = tri[j + 4 * t1];
                else
                    v4 = tri[j + 4 * t1];
            }
            // Check if the 4 verticies are coplanar (the resulting simplex will have a volume of 0)
            TYPE x2[12];
            x2[0]      = x[0 + 3 * v1];
            x2[1]      = x[1 + 3 * v1];
            x2[2]      = x[2 + 3 * v1];
            x2[3]      = x[0 + 3 * v2];
            x2[4]      = x[1 + 3 * v2];
            x2[5]      = x[2 + 3 * v2];
            x2[6]      = x[0 + 3 * v3];
            x2[7]      = x[1 + 3 * v3];
            x2[8]      = x[2 + 3 * v3];
            x2[9]      = x[0 + 3 * v4];
            x2[10]     = x[1 + 3 * v4];
            x2[11]     = x[2 + 3 * v4];
            double vol = fabs( calc_volume<3, TYPE, ETYPE>( x2 ) );
            if ( vol > TOL_VOL ) {
                // The points are not coplanar
                continue;
            }
            // Create the new triangles that are formed by the 2-2 edge flip
            int v5       = tri[i1 + 4 * t1];
            index_old[0] = t1;
            index_old[1] = t2;
            new_tri[0]   = v1;
            new_tri[1]   = v2;
            new_tri[2]   = v3;
            new_tri[3]   = v5; // The first triangle
            new_tri[4]   = v1;
            new_tri[5]   = v2;
            new_tri[6]   = v4;
            new_tri[7]   = v5; // The second triangle
            // Check that the new triagles are valid (this can occur if they are planar)
            bool isvalid = true;
            for ( int it = 0; it < 2; it++ ) {
                for ( int i = 0; i < 4; i++ ) {
                    int k         = new_tri[i + it * 4];
                    x2[0 + i * 3] = x[0 + k * 3];
                    x2[1 + i * 3] = x[1 + k * 3];
                    x2[2 + i * 3] = x[2 + k * 3];
                }
                double volume = calc_volume<3, TYPE, ETYPE>( x2 );
                if ( fabs( volume ) <= TOL_VOL ) {
                    // The triangle is invalid (collinear)
                    isvalid = false;
                } else if ( volume < 0 ) {
                    // The ordering of the points is invalid, swap the last two points (we want the
                    // volume to be positive)
                    int tmp             = new_tri[2 + it * 4];
                    new_tri[2 + it * 4] = new_tri[3 + it * 4];
                    new_tri[3 + it * 4] = tmp;
                }
            }
            if ( !isvalid )
                continue;
            // Check that the new triangles are Delanuay (we have already checked that they are
            // valid)
            TYPE xv[3];
            for ( int j = 0; j < 4; j++ ) {
                int k         = new_tri[j];
                x2[0 + 3 * j] = x[0 + 3 * k];
                x2[1 + 3 * j] = x[1 + 3 * k];
                x2[2 + 3 * j] = x[2 + 3 * k];
            }
            xv[0]    = x[0 + 3 * v4];
            xv[1]    = x[1 + 3 * v4];
            xv[2]    = x[2 + 3 * v4];
            int test = test_in_circumsphere<3, TYPE, ETYPE>( x2, xv, TOL_VOL );
            if ( test == 1 ) {
                // The flip did not fix the Delaunay condition
                continue;
            }
            // Compute the triangle neighbors (this is a less efficient, but general method)
            const int N_old = 2;
            const int N_new = 2;
            int f1, f2;
            for ( int j = 0; j < N_new * 4; j++ )
                new_tri_nab[j] = -1;
            for ( int j1 = 0; j1 < N_new; j1++ ) {
                for ( int j2 = j1 + 1; j2 < N_new;
                      j2++ ) { // We only need to loop through unchecked triangles
                    bool are_neighbors =
                        are_tri_neighbors( 3, &new_tri[j1 * 4], &new_tri[j2 * 4], &f1, &f2 );
                    if ( are_neighbors ) {
                        new_tri_nab[j1 * 4 + f1] = -j2 - 2;
                        new_tri_nab[j2 * 4 + f2] = -j1 - 2;
                    } else {
                        printp( "Internal error (flip_3D_22)\n" );
                        return false;
                    }
                }
            }
            int neighbors[8];
            for ( int j1 = 0; j1 < N_old; j1++ ) {
                for ( int j2 = 0; j2 < 4; j2++ )
                    neighbors[j2 + j1 * 4] = tri_nab[j2 + index_old[j1] * 4];
            }
            for ( auto &neighbor : neighbors ) {
                for ( int j2 = 0; j2 < N_old; j2++ ) {
                    if ( neighbor == index_old[j2] )
                        neighbor = -1;
                }
            }
            for ( int j1 = 0; j1 < N_new; j1++ ) {
                for ( auto &neighbor : neighbors ) {
                    if ( neighbor == -1 )
                        continue;
                    bool are_neighbors =
                        are_tri_neighbors( 3, &new_tri[j1 * 4], &tri[neighbor * 4], &f1, &f2 );
                    if ( are_neighbors )
                        new_tri_nab[j1 * 4 + f1] = neighbor;
                }
            }
            return true;
        }
    }
    return false;
}


/************************************************************************
 * This function performs a 3-2 flip in 3D                               *
 * Note: all indicies are hard-coded for ndim=3 and some loops may be    *
 * unrolled to simplify the code or improve performance                  *
 ************************************************************************/
template<class TYPE, class ETYPE>
bool DelaunayTessellation::flip_3D_32( const TYPE x[],
                                       const int tri[],
                                       const int tri_nab[],
                                       const int t1,
                                       const int,
                                       const int t2,
                                       const int,
                                       int *index_old,
                                       int *new_tri,
                                       int *new_tri_nab,
                                       const double TOL_VOL )
{
    // Search for triangles that are neighbors to both t1 and t2
    int nab_list[4];
    int N_nab = 0;
    for ( int i1 = 0; i1 < 4; i1++ ) {
        for ( int i2 = 0; i2 < 4; i2++ ) {
            if ( tri_nab[i1 + 4 * t1] == tri_nab[i2 + 4 * t2] && tri_nab[i1 + 4 * t1] != -1 ) {
                nab_list[N_nab] = tri_nab[i1 + 4 * t1];
                N_nab++;
            }
        }
    }
    if ( N_nab == 0 ) {
        // No triangle neighbors both t1 and t2
        return false;
    }
    for ( int i = 0; i < N_nab; i++ ) {
        // The nodes that are common to all 3 triangles are the unique verticies on the 3 new
        // triangles
        int nodes[4];
        for ( int j = 0; j < 4; j++ )
            nodes[j] = tri[j + 4 * nab_list[i]];
        for ( auto &node : nodes ) {
            bool found1 = false;
            bool found2 = false;
            for ( int j2 = 0; j2 < 4; j2++ ) {
                if ( node == tri[j2 + 4 * t1] )
                    found1 = true;
                if ( node == tri[j2 + 4 * t2] )
                    found2 = true;
            }
            if ( !found1 || !found2 )
                node = -1;
        }
        int N_nodes = 0;
        for ( auto &node : nodes ) {
            if ( node != -1 ) {
                nodes[N_nodes] = node;
                N_nodes++;
            }
        }
        if ( N_nodes != 2 ) {
            // This should not occur
            printp( "Unexpected error\n" );
            return false;
        }
        int t[3];
        t[0] = t1;
        t[1] = t2;
        t[2] = nab_list[i];
        int surf[4];
        int k = 0;
        for ( auto &elem : t ) {
            for ( int j2 = 0; j2 < 4; j2++ ) {
                int tmp = tri[j2 + 4 * elem];
                if ( tmp != nodes[0] && tmp != nodes[1] ) {
                    bool in_surf = false;
                    for ( int j3 = 0; j3 < k; j3++ ) {
                        if ( tmp == surf[j3] )
                            in_surf = true;
                    }
                    if ( !in_surf ) {
                        surf[k] = tmp;
                        k++;
                    }
                }
            }
        }
        if ( k != 3 ) {
            // This should not occur
            printp( "Unexpected error\n" );
            return false;
        }
        // Create the new triangles
        index_old[0] = t1;
        index_old[1] = t2;
        index_old[2] = nab_list[i];
        new_tri[0]   = nodes[0];
        new_tri[1]   = surf[0];
        new_tri[2]   = surf[1];
        new_tri[3]   = surf[2];
        new_tri[4]   = nodes[1];
        new_tri[5]   = surf[0];
        new_tri[6]   = surf[1];
        new_tri[7]   = surf[2];
        // Check that the new triangles are valid (this can occur if they are planar)
        TYPE x2[12];
        bool isvalid = true;
        for ( int it = 0; it < 2; it++ ) {
            for ( int i = 0; i < 4; i++ ) {
                int k         = new_tri[i + it * 4];
                x2[0 + i * 3] = x[0 + k * 3];
                x2[1 + i * 3] = x[1 + k * 3];
                x2[2 + i * 3] = x[2 + k * 3];
            }
            double volume = calc_volume<3, TYPE, ETYPE>( x2 );
            if ( fabs( volume ) <= TOL_VOL ) {
                // The triangle is invalid (collinear)
                isvalid = false;
            } else if ( volume < 0 ) {
                // The ordering of the points is invalid, swap the last two points (we want the
                // volume to be positive)
                int tmp             = new_tri[2 + it * 4];
                new_tri[2 + it * 4] = new_tri[3 + it * 4];
                new_tri[3 + it * 4] = tmp;
            }
        }
        if ( !isvalid )
            continue;
        // Check that the new triangles are Delanuay (we have already checked that they are valid)
        TYPE xv[3];
        for ( int j = 0; j < 4; j++ ) {
            int k         = new_tri[j];
            x2[0 + 3 * j] = x[0 + 3 * k];
            x2[1 + 3 * j] = x[1 + 3 * k];
            x2[2 + 3 * j] = x[2 + 3 * k];
        }
        xv[0]    = x[0 + 3 * nodes[1]];
        xv[1]    = x[1 + 3 * nodes[1]];
        xv[2]    = x[2 + 3 * nodes[1]];
        int test = test_in_circumsphere<3, TYPE, ETYPE>( x2, xv, TOL_VOL );
        if ( test == 1 ) {
            // The new triangles did not fixed the surface
            continue;
        }
        // Compute the triangle neighbors (this is a less efficient, but general method)
        const int N_old = 3;
        const int N_new = 2;
        int f1, f2;
        for ( int j = 0; j < N_new * 4; j++ )
            new_tri_nab[j] = -1;
        for ( int j1 = 0; j1 < N_new; j1++ ) {
            for ( int j2 = j1 + 1; j2 < N_new;
                  j2++ ) { // We only need to loop through unchecked triangles
                bool are_neighbors =
                    are_tri_neighbors( 3, &new_tri[j1 * 4], &new_tri[j2 * 4], &f1, &f2 );
                if ( are_neighbors ) {
                    new_tri_nab[j1 * 4 + f1] = -j2 - 2;
                    new_tri_nab[j2 * 4 + f2] = -j1 - 2;
                } else {
                    printp( "Internal error (flip_3D_32)\n" );
                    return false;
                }
            }
        }
        int neighbors[12];
        for ( int j1 = 0; j1 < N_old; j1++ ) {
            for ( int j2 = 0; j2 < 4; j2++ )
                neighbors[j2 + j1 * 4] = tri_nab[j2 + index_old[j1] * 4];
        }
        for ( auto &neighbor : neighbors ) {
            for ( int j2 = 0; j2 < N_old; j2++ ) {
                if ( neighbor == index_old[j2] )
                    neighbor = -1;
            }
        }
        for ( int j1 = 0; j1 < N_new; j1++ ) {
            for ( auto &neighbor : neighbors ) {
                if ( neighbor == -1 )
                    continue;
                bool are_neighbors =
                    are_tri_neighbors( 3, &new_tri[j1 * 4], &tri[neighbor * 4], &f1, &f2 );
                if ( are_neighbors )
                    new_tri_nab[j1 * 4 + f1] = neighbor;
            }
        }
        // Finished
        return true;
    }
    return false;
}


/************************************************************************
 * This function performs a 2-3 flip in 3D                               *
 * Note: all indicies are hard-coded for ndim=3 and some loops may be    *
 * unrolled to simplify the code or improve performance                  *
 ************************************************************************/
template<class TYPE, class ETYPE>
bool DelaunayTessellation::flip_3D_23( const TYPE x[],
                                       const int tri[],
                                       const int tri_nab[],
                                       const int t1,
                                       const int s1,
                                       const int t2,
                                       const int s2,
                                       int *index_old,
                                       int *new_tri,
                                       int *new_tri_nab,
                                       const double TOL_VOL )
{
    // First lets check if the flip is valid
    TYPE x2[12], xi[3];
    for ( int i = 0; i < 4; i++ ) {
        int k         = tri[i + 4 * t1];
        x2[0 + 3 * i] = x[0 + 3 * k];
        x2[1 + 3 * i] = x[1 + 3 * k];
        x2[2 + 3 * i] = x[2 + 3 * k];
    }
    int k        = tri[s2 + 4 * t2];
    xi[0]        = x[0 + 3 * k];
    xi[1]        = x[1 + 3 * k];
    xi[2]        = x[2 + 3 * k];
    bool isvalid = test_flip_valid<3, TYPE, ETYPE>( x2, s1, xi );
    if ( !isvalid ) {
        // We need to do an edge flip, and the edge flip is not valid
        return false;
    }
    // Form the new triangles
    int is[3], in[2] = { 0, 0 };
    k = 0;
    for ( int i = 0; i < 4; i++ ) {
        if ( i == s1 ) {
            in[0] = tri[i + 4 * t1];
        } else {
            is[k] = tri[i + 4 * t1];
            k++;
        }
    }
    in[1]        = tri[s2 + 4 * t2];
    index_old[0] = t1;
    index_old[1] = t2;
    new_tri[0]   = in[0];
    new_tri[1]   = in[1];
    new_tri[2]   = is[0];
    new_tri[3]   = is[1];
    new_tri[4]   = in[0];
    new_tri[5]   = in[1];
    new_tri[6]   = is[0];
    new_tri[7]   = is[2];
    new_tri[8]   = in[0];
    new_tri[9]   = in[1];
    new_tri[10]  = is[1];
    new_tri[11]  = is[2];
    // Check that the new triagles are valid (this can occur if they are planar)
    isvalid = true;
    for ( int it = 0; it < 3; it++ ) {
        for ( int i = 0; i < 4; i++ ) {
            int k         = new_tri[i + it * 4];
            x2[0 + i * 3] = x[0 + k * 3];
            x2[1 + i * 3] = x[1 + k * 3];
            x2[2 + i * 3] = x[2 + k * 3];
        }
        double volume = calc_volume<3, TYPE, ETYPE>( x2 );
        if ( fabs( volume ) <= TOL_VOL ) {
            // The triangle is invalid (collinear)
            isvalid = false;
        } else if ( volume < 0 ) {
            // The ordering of the points is invalid, swap the last two points (we want the volume
            // to be positive)
            int tmp             = new_tri[2 + it * 4];
            new_tri[2 + it * 4] = new_tri[3 + it * 4];
            new_tri[3 + it * 4] = tmp;
        }
    }
    if ( !isvalid )
        return false;
    // Check that the new triangles are Delanuay (we have already checked that they are valid)
    for ( int it = 0; it < 3; it++ ) {
        for ( int i = 0; i < 4; i++ ) {
            int k         = new_tri[i + it * 4];
            x2[0 + i * 3] = x[0 + k * 3];
            x2[1 + i * 3] = x[1 + k * 3];
            x2[2 + i * 3] = x[2 + k * 3];
        }
        int k;
        if ( it == 0 ) {
            k = is[2];
        } else if ( it == 1 ) {
            k = is[1];
        } else {
            k = is[0];
        }
        xi[0]    = x[0 + k * 3];
        xi[1]    = x[1 + k * 3];
        xi[2]    = x[2 + k * 3];
        int test = test_in_circumsphere<3, TYPE, ETYPE>( x2, xi, TOL_VOL );
        if ( test == 1 ) {
            // The new triangles did not fix the surface
            isvalid = false;
        }
    }
    if ( !isvalid )
        return false;
    // Compute the triangle neighbors (this is a less efficient, but general method)
    const int N_old = 2;
    const int N_new = 3;
    int f1, f2;
    for ( int j = 0; j < N_new * 4; j++ )
        new_tri_nab[j] = -1;
    for ( int j1 = 0; j1 < N_new; j1++ ) {
        for ( int j2 = j1 + 1; j2 < N_new;
              j2++ ) { // We only need to loop through unchecked triangles
            bool are_neighbors =
                are_tri_neighbors( 3, &new_tri[j1 * 4], &new_tri[j2 * 4], &f1, &f2 );
            if ( are_neighbors ) {
                new_tri_nab[j1 * 4 + f1] = -j2 - 2;
                new_tri_nab[j2 * 4 + f2] = -j1 - 2;
            } else {
                printp( "Internal error (flip_3D_23)\n" );
                return false;
            }
        }
    }
    int neighbors[8];
    for ( int j1 = 0; j1 < N_old; j1++ ) {
        for ( int j2 = 0; j2 < 4; j2++ )
            neighbors[j2 + j1 * 4] = tri_nab[j2 + index_old[j1] * 4];
    }
    for ( auto &neighbor : neighbors ) {
        for ( int j2 = 0; j2 < N_old; j2++ ) {
            if ( neighbor == index_old[j2] )
                neighbor = -1;
        }
    }
    for ( int j1 = 0; j1 < N_new; j1++ ) {
        for ( auto &neighbor : neighbors ) {
            if ( neighbor == -1 )
                continue;
            bool are_neighbors =
                are_tri_neighbors( 3, &new_tri[j1 * 4], &tri[neighbor * 4], &f1, &f2 );
            if ( are_neighbors )
                new_tri_nab[j1 * 4 + f1] = neighbor;
        }
    }
    return true;
}


/************************************************************************
 * This function performs a 4-4 flip in 3D                               *
 * Note: all indicies are hard-coded for ndim=3 and some loops may be    *
 * unrolled to simplify the code or improve performance                  *
 ************************************************************************/
template<class TYPE, class ETYPE>
bool DelaunayTessellation::flip_3D_44( const TYPE x[],
                                       const int tri[],
                                       const int tri_nab[],
                                       const int t1,
                                       const int s1,
                                       const int t2,
                                       const int s2,
                                       int *index_old,
                                       int *new_tri,
                                       int *new_tri_nab,
                                       const double TOL_VOL )
{
    // Loop through the triangle neighbors for both t1 and t2, to find
    // a pair of triangles that form and octahedron
    for ( int i1 = 0; i1 < 4; i1++ ) {
        int t3 = tri_nab[i1 + 4 * t1];
        for ( int i2 = 0; i2 < 4; i2++ ) {
            int t4 = tri_nab[i2 + 4 * t2];
            if ( t3 == -1 || t3 == t2 || t4 == -1 || t4 == t1 || t3 == t4 ) {
                // We do not have 4 unique triangles
                continue;
            }
            if ( tri_nab[0 + 4 * t3] != t4 && tri_nab[1 + 4 * t3] != t4 &&
                 tri_nab[2 + 4 * t3] != t4 && tri_nab[3 + 4 * t3] != t4 ) {
                // The 4 triangles do not form an octahedron
                continue;
            }
            // Get a list of the surface ids
            int s12, s13, s21, s24, s31 = 0, s34 = 0, s42 = 0, s43 = 0;
            s12 = s1;
            s13 = i1;
            s21 = s2;
            s24 = i2;
            for ( int i = 0; i < 4; i++ ) {
                if ( tri_nab[i + 4 * t3] == t1 )
                    s31 = i;
                if ( tri_nab[i + 4 * t3] == t4 )
                    s34 = i;
                if ( tri_nab[i + 4 * t4] == t2 )
                    s42 = i;
                if ( tri_nab[i + 4 * t4] == t3 )
                    s43 = i;
            }
            // Get a list of the unique nodes
            int nodes[16];
            for ( auto &node : nodes )
                node = -1;
            for ( int i = 0; i < 4; i++ )
                nodes[i] = tri[i + 4 * t1];
            int N_nodes = 4;
            for ( int j1 = 0; j1 < 4; j1++ ) {
                bool found = false;
                for ( int j2 = 0; j2 < N_nodes; j2++ ) {
                    if ( tri[j1 + 4 * t2] == nodes[j2] )
                        found = true;
                }
                if ( !found ) {
                    nodes[N_nodes] = tri[j1 + 4 * t2];
                    N_nodes++;
                }
            }
            for ( int j1 = 0; j1 < 4; j1++ ) {
                bool found = false;
                for ( int j2 = 0; j2 < N_nodes; j2++ ) {
                    if ( tri[j1 + 4 * t3] == nodes[j2] )
                        found = true;
                }
                if ( !found ) {
                    nodes[N_nodes] = tri[j1 + 4 * t3];
                    N_nodes++;
                }
            }
            for ( int j1 = 0; j1 < 4; j1++ ) {
                bool found = false;
                for ( int j2 = 0; j2 < N_nodes; j2++ ) {
                    if ( tri[j1 + 4 * t4] == nodes[j2] )
                        found = true;
                }
                if ( !found ) {
                    nodes[N_nodes] = tri[j1 + 4 * t4];
                    N_nodes++;
                }
            }
            if ( N_nodes != 6 ) {
                // We still do not have a valid octahedron
                printp( "Unexpected number of nodes\n" );
                continue;
            }
            // Check if any 4 of the nodes are coplanar (the resulting simplex will have a volume of
            // 0)
            // Note: while there are 15 combinations of 4 points, in reality
            // because the triangle faces already specify groups of 3 that
            // are coplanar, we only need to check 2 groups of points
            int set1[4], set2[4];
            for ( int i = 0; i < 4; i++ ) {
                set1[i] = tri[i + 4 * t1];
                set2[i] = tri[i + 4 * t1];
            }
            set1[s13] = tri[s21 + t2 * 4];
            set2[s12] = tri[s31 + t3 * 4];
            TYPE x2[12], xi[3];
            double vol1, vol2;
            for ( int j1 = 0; j1 < 4; j1++ ) {
                for ( int j2 = 0; j2 < 3; j2++ )
                    x2[j2 + 3 * j1] = x[j2 + 3 * set1[j1]];
            }
            vol1 = fabs( calc_volume<3, TYPE, ETYPE>( x2 ) );
            for ( int j1 = 0; j1 < 4; j1++ ) {
                for ( int j2 = 0; j2 < 3; j2++ )
                    x2[j2 + 3 * j1] = x[j2 + 3 * set2[j1]];
            }
            vol2 = fabs( calc_volume<3, TYPE, ETYPE>( x2 ) );
            for ( int it = 0; it < 2; it++ ) { // Loop through the sets (2)
                if ( it == 0 ) {
                    // Check set 1
                    if ( vol1 > TOL_VOL ) {
                        // Set 1 is not coplanar
                        continue;
                    }
                    // Set 1 is coplanar, check if the line between the two unused points passes
                    // through the 4 points in the set. If it does, then test_flip_valid will be
                    // true
                    // using either the current triangle and neighbor, or using the other 2
                    // triangles
                    for ( int j1 = 0; j1 < 4; j1++ ) {
                        int k = tri[j1 + 4 * t1];
                        for ( int j2 = 0; j2 < 3; j2++ )
                            x2[j2 + 3 * j1] = x[j2 + 3 * k];
                    }
                    for ( int j = 0; j < 3; j++ )
                        xi[j] = x[j + 3 * tri[s21 + 4 * t2]];
                    bool is_valid1 = test_flip_valid<3, TYPE, ETYPE>( x2, s12, xi );
                    for ( int j1 = 0; j1 < 4; j1++ ) {
                        int k = tri[j1 + 4 * t3];
                        for ( int j2 = 0; j2 < 3; j2++ )
                            x2[j2 + 3 * j1] = x[j2 + 3 * k];
                    }
                    bool is_valid2 = test_flip_valid<3, TYPE, ETYPE>( x2, s34, xi );
                    if ( !is_valid1 && !is_valid2 ) {
                        // The flip is not valid
                        continue;
                    }
                } else if ( it == 1 ) {
                    // Check set 2
                    if ( vol2 > TOL_VOL ) {
                        // Set 1 is not coplanar
                        continue;
                    }
                    // Set 2 is coplanar, check if the line between the two unused points passes
                    // through the 4 points in the set. If it does, then test_flip_valid will be
                    // true
                    // using either the current triangle and neighbor, or using the other 2
                    // triangles
                    for ( int j1 = 0; j1 < 4; j1++ ) {
                        int k = tri[j1 + 4 * t1];
                        for ( int j2 = 0; j2 < 3; j2++ )
                            x2[j2 + 3 * j1] = x[j2 + 3 * k];
                    }
                    for ( int j = 0; j < 3; j++ )
                        xi[j] = x[j + 3 * tri[s31 + 4 * t3]];
                    bool is_valid1 = test_flip_valid<3, TYPE, ETYPE>( x2, s13, xi );
                    for ( int j1 = 0; j1 < 4; j1++ ) {
                        int k = tri[j1 + 4 * t2];
                        for ( int j2 = 0; j2 < 3; j2++ )
                            x2[j2 + 3 * j1] = x[j2 + 3 * k];
                    }
                    bool is_valid2 = test_flip_valid<3, TYPE, ETYPE>( x2, s24, xi );
                    if ( !is_valid1 && !is_valid2 ) {
                        // The flip is not valid
                        continue;
                    }
                } else {
                    printp( "Unexpected error\n" );
                    return false;
                }
                // Get a list of the nodes common to all triangles
                int node_all[4];
                int k = 0;
                for ( int i = 0; i < 4; i++ ) {
                    int tmp     = tri[i + 4 * t1];
                    bool found2 = false;
                    for ( int j = 0; j < 4; j++ ) {
                        if ( tri[j + 4 * t2] == tmp )
                            found2 = true;
                    }
                    bool found3 = false;
                    for ( int j = 0; j < 4; j++ ) {
                        if ( tri[j + 4 * t3] == tmp )
                            found3 = true;
                    }
                    bool found4 = false;
                    for ( int j = 0; j < 4; j++ ) {
                        if ( tri[j + 4 * t3] == tmp )
                            found4 = true;
                    }
                    if ( found2 && found3 && found4 ) {
                        node_all[k] = tmp;
                        k++;
                    }
                }
                if ( k != 2 ) {
                    printp( "Unexpected error\n" );
                    return false;
                }
                // We have a valid flip, create the triangles
                index_old[0] = t1;
                index_old[1] = t2;
                index_old[2] = t3;
                index_old[3] = t4;
                if ( it == 0 ) {
                    // The 1-3 and 2-4 faces form the plane
                    new_tri[0]  = tri[s12 + 4 * t1];
                    new_tri[1]  = tri[s21 + 4 * t2];
                    new_tri[2]  = tri[s13 + 4 * t1];
                    new_tri[3]  = node_all[0];
                    new_tri[4]  = tri[s12 + 4 * t1];
                    new_tri[5]  = tri[s21 + 4 * t2];
                    new_tri[6]  = tri[s13 + 4 * t1];
                    new_tri[7]  = node_all[1];
                    new_tri[8]  = tri[s34 + 4 * t3];
                    new_tri[9]  = tri[s43 + 4 * t4];
                    new_tri[10] = tri[s31 + 4 * t3];
                    new_tri[11] = node_all[0];
                    new_tri[12] = tri[s34 + 4 * t3];
                    new_tri[13] = tri[s43 + 4 * t4];
                    new_tri[14] = tri[s31 + 4 * t3];
                    new_tri[15] = node_all[1];
                } else {
                    // The 1-2 and 3-4 faces form the plane
                    new_tri[0]  = tri[s13 + 4 * t1];
                    new_tri[1]  = tri[s31 + 4 * t3];
                    new_tri[2]  = tri[s12 + 4 * t1];
                    new_tri[3]  = node_all[0];
                    new_tri[4]  = tri[s13 + 4 * t1];
                    new_tri[5]  = tri[s31 + 4 * t3];
                    new_tri[6]  = tri[s12 + 4 * t1];
                    new_tri[7]  = node_all[1];
                    new_tri[8]  = tri[s24 + 4 * t2];
                    new_tri[9]  = tri[s42 + 4 * t4];
                    new_tri[10] = tri[s21 + 4 * t2];
                    new_tri[11] = node_all[0];
                    new_tri[12] = tri[s24 + 4 * t2];
                    new_tri[13] = tri[s42 + 4 * t4];
                    new_tri[14] = tri[s21 + 4 * t2];
                    new_tri[15] = node_all[1];
                }
                // Check that the new triagles are valid (this can occur if they are planar)
                bool isvalid = true;
                for ( int it2 = 0; it2 < 4; it2++ ) {
                    for ( int i = 0; i < 4; i++ ) {
                        int k         = new_tri[i + it2 * 4];
                        x2[0 + i * 3] = x[0 + k * 3];
                        x2[1 + i * 3] = x[1 + k * 3];
                        x2[2 + i * 3] = x[2 + k * 3];
                    }
                    double volume = calc_volume<3, TYPE, ETYPE>( x2 );
                    if ( fabs( volume ) <= TOL_VOL ) {
                        // The triangle is invalid (coplanar)
                        isvalid = false;
                    } else if ( volume < 0 ) {
                        // The ordering of the points is invalid, swap the last two points (we want
                        // the volume to be positive)
                        int tmp              = new_tri[2 + it2 * 4];
                        new_tri[2 + it2 * 4] = new_tri[3 + it2 * 4];
                        new_tri[3 + it2 * 4] = tmp;
                    }
                }
                if ( !isvalid )
                    continue;
                // Check that the new triangles are Delanuay (we have already checked that they are
                // valid)
                for ( int it2 = 0; it2 < 4; it2++ ) { // Loop through the test cases
                    for ( int i = 0; i < 4; i++ ) {
                        int k         = new_tri[i + it2 * 4];
                        x2[0 + i * 3] = x[0 + k * 3];
                        x2[1 + i * 3] = x[1 + k * 3];
                        x2[2 + i * 3] = x[2 + k * 3];
                    }
                    int k = -1;
                    if ( it2 == 0 || it2 == 2 )
                        k = new_tri[7];
                    else if ( it2 == 1 || it2 == 3 )
                        k = new_tri[3];
                    xi[0]    = x[0 + 3 * k];
                    xi[1]    = x[1 + 3 * k];
                    xi[2]    = x[2 + 3 * k];
                    int test = test_in_circumsphere<3, TYPE, ETYPE>( x2, xi, TOL_VOL );
                    if ( test == 1 )
                        isvalid = false;
                }
                if ( !isvalid )
                    continue;
                // Compute the triangle neighbors (this is a less efficient, but general method)
                const int N_old = 4;
                const int N_new = 4;
                int f1, f2;
                for ( int j = 0; j < N_new * 4; j++ )
                    new_tri_nab[j] = -1;
                for ( int j1 = 0; j1 < N_new; j1++ ) {
                    for ( int j2 = j1 + 1; j2 < N_new;
                          j2++ ) { // We only need to loop through unchecked triangles
                        bool are_neighbors =
                            are_tri_neighbors( 3, &new_tri[j1 * 4], &new_tri[j2 * 4], &f1, &f2 );
                        if ( are_neighbors ) {
                            new_tri_nab[j1 * 4 + f1] = -j2 - 2;
                            new_tri_nab[j2 * 4 + f2] = -j1 - 2;
                        }
                    }
                }
                int neighbors[16];
                for ( int j1 = 0; j1 < N_old; j1++ ) {
                    for ( int j2 = 0; j2 < 4; j2++ )
                        neighbors[j2 + j1 * 4] = tri_nab[j2 + index_old[j1] * 4];
                }
                for ( auto &neighbor : neighbors ) {
                    for ( int j2 = 0; j2 < N_old; j2++ ) {
                        if ( neighbor == index_old[j2] )
                            neighbor = -1;
                    }
                }
                for ( int j1 = 0; j1 < N_new; j1++ ) {
                    for ( auto &neighbor : neighbors ) {
                        if ( neighbor == -1 )
                            continue;
                        bool are_neighbors =
                            are_tri_neighbors( 3, &new_tri[j1 * 4], &tri[neighbor * 4], &f1, &f2 );
                        if ( are_neighbors )
                            new_tri_nab[j1 * 4 + f1] = neighbor;
                    }
                }
                // Finished
                return true;
            }
        }
    }
    return false;
}


/********************************************************************
 * This is the constructor to inialize the FaceList class            *
 ********************************************************************/
static inline unsigned int get_hash_key( size_t index )
{
    auto key = static_cast<unsigned int>( index );
    key *= 0x9E3779B9; // 2^32*0.5*(sqrt(5)-1)
    key >>= 22;        // mod(key,1024)
    return key;
}
template<int NDIM, class TYPE, class ETYPE>
DelaunayTessellation::FaceList<NDIM, TYPE, ETYPE>::FaceList(
    const int Nx_in, const TYPE *x_in, const int tri0_id, const int tri0[], const TYPE TOL_VOL )
    : Nx( Nx_in ), x0( x_in ), TOL_vol( TOL_VOL )
{
    ASSERT( sizeof( unsigned int ) ==
            4 ); // We need an unsigned int to be 4 bytes (32 bits) for the hash function
    N_face = 0;
    for ( auto &elem : hash_table )
        elem = -1; // Initialize the hash table to -1
    size = 1024;   // We will maintain internal storage using blocks of 2^k
    data = new face_data_struct[size / sizeof( face_data_struct )];
    for ( size_t i = 0; i < size / sizeof( face_data_struct ); i++ )
        data[i].reset();
    for ( int i = 0; i <= NDIM; i++ ) {
        unsigned int key = get_hash_key( get_face_index( i, 0 ) );
        int k            = hash_table[key];
        hash_table[key]  = i;
        data[i].prev     = -static_cast<int>( key + 1 );
        data[i].next     = k;
        if ( k != -1 )
            data[k].prev = i;
        data[i].tri_id  = tri0_id;
        data[i].face_id = i;
        int j1          = 0;
        for ( int j2 = 0; j2 <= NDIM; j2++ ) {
            if ( j2 == i )
                continue;
            int k             = tri0[j2];
            data[i].index[j1] = k;
            for ( int d = 0; d < NDIM; d++ )
                data[i].x[j1][d] = x0[d + k * NDIM];
            j1++;
        }
        N_face++;
    }
    // Initialize the "center" of the volume
    for ( auto &elem : xc )
        elem = 0.0;
    for ( int i = 0; i <= NDIM; i++ ) {
        for ( int j = 0; j < NDIM; j++ ) {
            for ( int d = 0; d < NDIM; d++ )
                xc[d] += static_cast<double>( data[i].x[j][d] );
        }
    }
    for ( auto &elem : xc )
        elem /= ( NDIM ) * ( NDIM + 1 );
    check_data();
}


/********************************************************************
 * This is the de-constructor for the NodeTriList class              *
 ********************************************************************/
template<int NDIM, class TYPE, class ETYPE>
DelaunayTessellation::FaceList<NDIM, TYPE, ETYPE>::~FaceList()
{
    delete[] data;
    size   = 0;
    N_face = 0;
}


/********************************************************************
 * Delete a face                                                     *
 ********************************************************************/
template<int NDIM, class TYPE, class ETYPE>
void DelaunayTessellation::FaceList<NDIM, TYPE, ETYPE>::delete_faces( int N_delete, int *ids )
{
    // Sort the ids to delete
    quicksort( N_delete, ids );
    // We delete in reverse order so that we can safely swap entries
    for ( int i = N_delete - 1; i >= 0; i-- ) {
        // Delete the entry and remove it from the hash table
        int k = ids[i];
        if ( data[k].prev >= 0 ) {
            data[data[k].prev].next = data[k].next;
            if ( data[k].next != -1 )
                data[data[k].next].prev = data[k].prev;
        } else {
            unsigned int key = -data[k].prev - 1;
            int j            = data[k].next;
            hash_table[key]  = j;
            if ( j != -1 )
                data[j].prev = -static_cast<int>( key + 1 );
        }
        if ( k < N_face - 1 ) {
            // Move the last entry to the current entry and update the hash table
            data[k] = data[N_face - 1];
            if ( data[k].prev >= 0 ) {
                data[data[k].prev].next = k;
            } else {
                unsigned int key = -data[k].prev - 1;
                hash_table[key]  = k;
            }
            if ( data[k].next != -1 )
                data[data[k].next].prev = k;
        }
        data[N_face - 1].reset();
        N_face--;
    }
}


/********************************************************************
 * This function adds a node to the convex hull                      *
 * Note: currently this algorithum scales as O(N)                    *
 ********************************************************************/
template<int NDIM, class TYPE, class ETYPE>
int DelaunayTessellation::FaceList<NDIM, TYPE, ETYPE>::add_node( const int node_id,
                                                                 std::vector<size_t> &unused,
                                                                 size_t &N_tri,
                                                                 unsigned int *new_tri_id,
                                                                 int *new_tri,
                                                                 int *new_tri_nab,
                                                                 int *neighbor,
                                                                 int *face_id )
{
    PROFILE_START_L2( "FaceList::add_node" );
#if DEBUG_CHECK == 2
    check_data();
#endif
    // Create the new triangles
    auto ids      = new int[N_face];
    int N_tri_new = 0;
    TYPE x2[NDIM * ( NDIM + 1 )];
    for ( int i = 0; i < N_face; i++ ) {
        // Get the triangle number and face id
        int tri_num  = data[i].tri_id;
        int face_num = data[i].face_id;
// We will need to for a new triangle with each face where the distance
// from the face to the point is > 0 and the volume is >= TOL_vol
#if 0
            bool pass = calc_surface_distance(data[i].x,&x0[node_id*NDIM]) > 0;
#else
        bool pass = outside_triangle( data[i].x, &x0[node_id * NDIM] );
#endif
        if ( !pass )
            continue;
        for ( int j1 = 0; j1 < NDIM; j1++ ) {
            for ( int j2 = 0; j2 < NDIM; j2++ )
                x2[j2 + j1 * NDIM] = data[i].x[j1][j2];
        }
        for ( int j2 = 0; j2 < NDIM; j2++ )
            x2[j2 + NDIM * NDIM] = x0[j2 + node_id * NDIM];
        double volume = fabs( calc_volume<NDIM, TYPE, ETYPE>( x2 ) );
        if ( volume <= TOL_vol ) {
            // printf("%e %e\n",dist,volume);
            continue;
        }
        // Create a new triangle consisting of the current point and the current face
        for ( int j = 0; j < NDIM; j++ )
            new_tri[j + N_tri_new * ( NDIM + 1 )] = data[i].index[j];
        new_tri[NDIM + N_tri_new * ( NDIM + 1 )] = node_id;
        for ( int j = 0; j < NDIM; j++ )
            new_tri_nab[j + N_tri_new * ( NDIM + 1 )] = -1;
        new_tri_nab[NDIM + N_tri_new * ( NDIM + 1 )] = tri_num;
        neighbor[N_tri_new]                          = tri_num;
        face_id[N_tri_new]                           = face_num;
        if ( !unused.empty() ) {
            new_tri_id[N_tri_new] = unused.back();
            unused.pop_back();
        } else {
            new_tri_id[N_tri_new] = N_tri;
            N_tri++;
        }
        ids[N_tri_new] = i;
        N_tri_new++;
    }
    if ( N_tri_new == 0 ) {
        // No triangles were created, the point must be inside (or on) the convex hull
        delete[] ids;
        PROFILE_STOP2_L2( "FaceList::add_node" );
        return 0;
    }
    // Check that the new triangles are valid and have the proper point ordering
    for ( int i = 0; i < N_tri_new; i++ ) {
        for ( int j1 = 0; j1 <= NDIM; j1++ ) {
            int k = new_tri[j1 + i * ( NDIM + 1 )];
            for ( int j2 = 0; j2 < NDIM; j2++ )
                x2[j2 + j1 * NDIM] = x0[j2 + k * NDIM];
        }
        double volume = calc_volume<NDIM, TYPE, ETYPE>( x2 );
        if ( fabs( volume ) <= TOL_vol ) {
            // The triangle is invalid
            delete[] ids;
            PROFILE_STOP2_L2( "FaceList::add_node" );
            return -1;
        } else if ( volume < 0 ) {
            // The volume is negitive, swap the order of the last two points
            std::swap( new_tri[NDIM - 1 + i * ( NDIM + 1 )], new_tri[NDIM + i * ( NDIM + 1 )] );
            std::swap( new_tri_nab[NDIM - 1 + i * ( NDIM + 1 )],
                       new_tri_nab[NDIM + i * ( NDIM + 1 )] );
        }
        for ( int j1 = 0; j1 <= NDIM; j1++ ) {
            int k = new_tri[j1 + i * ( NDIM + 1 )];
            for ( int j2 = 0; j2 < NDIM; j2++ )
                x2[j2 + j1 * NDIM] = x0[j2 + k * NDIM];
        }
        volume = calc_volume<NDIM, TYPE, ETYPE>( x2 );
        if ( volume <= 0.0 )
            printp( "Warning: volume is still negitive\n" );
    }
    // Check if the current triangle shares a face with one of the other new triangles
    for ( int i1 = 0; i1 < N_tri_new; i1++ ) {
        for ( int i2 = i1 + 1; i2 < N_tri_new; i2++ ) {
            int f1, f2;
            bool is_neighbor = are_tri_neighbors(
                NDIM, &new_tri[i1 * ( NDIM + 1 )], &new_tri[i2 * ( NDIM + 1 )], &f1, &f2 );
            if ( is_neighbor ) {
#if DEBUG_CHECK > 0
                ASSERT( new_tri_nab[f1 + i1 * ( NDIM + 1 )] == -1 );
                ASSERT( new_tri_nab[f2 + i2 * ( NDIM + 1 )] == -1 );
#endif
                new_tri_nab[f1 + i1 * ( NDIM + 1 )] = new_tri_id[i2];
                new_tri_nab[f2 + i2 * ( NDIM + 1 )] = new_tri_id[i1];
            }
        }
    }
    // Remove the old faces that were used to create the triangles
    delete_faces( N_tri_new, ids );
    // Allocate more memory if necessary (we will add a max of N_tri_new*NDIM faces)
    if ( ( N_face + N_tri_new * NDIM ) > (int) ( size / sizeof( face_data_struct ) ) ) {
        while ( ( N_face + N_tri_new * NDIM ) > (int) ( size / sizeof( face_data_struct ) ) )
            size *= 2;
        face_data_struct *tmp_data = data;
        data                       = new face_data_struct[size / sizeof( face_data_struct )];
        for ( size_t i = 0; i < size / sizeof( face_data_struct ); i++ )
            data[i].reset();
        for ( int i = 0; i < N_face; i++ )
            data[i] = tmp_data[i];
        delete[] tmp_data;
    }
    // Add the new faces for each triangle
    int N_add = 0;
    for ( int i = 0; i < N_tri_new * ( NDIM + 1 ); i++ ) {
        if ( new_tri_nab[i] == -1 )
            N_add++;
    }
    if ( NDIM == 2 ) {
        ASSERT( N_add == 2 );
    } else if ( NDIM == 3 ) {
        ASSERT( N_add >= 3 );
        ASSERT( N_add <= 2 * N_tri_new + 1 );
    }
    for ( int i = 0; i < N_tri_new; i++ ) {
        for ( int j = 0; j <= NDIM; j++ ) {
            if ( new_tri_nab[j + i * ( NDIM + 1 )] == -1 ) {
                // Add the face
                data[N_face].tri_id  = new_tri_id[i];
                data[N_face].face_id = j;
                int j1               = 0;
                for ( int j2 = 0; j2 <= NDIM; j2++ ) {
                    if ( j2 == j )
                        continue;
                    int k                  = new_tri[j2 + i * ( NDIM + 1 )];
                    data[N_face].index[j1] = k;
                    for ( int d = 0; d < NDIM; d++ )
                        data[N_face].x[j1][d] = x0[d + k * NDIM];
                    j1++;
                }
                // Update the hash table
                unsigned int key  = get_hash_key( get_face_index( j, new_tri_id[i] ) );
                int k             = hash_table[key];
                hash_table[key]   = N_face;
                data[N_face].prev = -static_cast<int>( key + 1 );
                data[N_face].next = k;
                if ( k != -1 )
                    data[k].prev = N_face;
                N_face++;
            }
        }
    }
#if DEBUG_CHECK == 2
    check_data();
#endif
    delete[] ids;
// Check the triangle neighbors
#if DEBUG_CHECK == 2
    int *tri_nab_tmp = new int[N_tri_new * ( NDIM + 1 )];
    memcpy( tri_nab_tmp, new_tri_nab, N_tri_new * ( NDIM + 1 ) * sizeof( int ) );
    for ( int i = 0; i < N_tri_new * ( NDIM + 1 ); i++ ) {
        int k = -1;
        for ( int j = 0; j < N_tri_new; j++ )
            if ( tri_nab_tmp[i] == (int) new_tri_id[j] )
                k = j;
        tri_nab_tmp[i] = k;
    }
    bool all_valid = check_current_triangles<NDIM, TYPE, ETYPE>(
        Nx, x0, N_tri_new, new_tri, tri_nab_tmp, std::vector<size_t>(), 0 );
    if ( !all_valid )
        printp( "Warning: new triangle neighbors are inconsistent\n" );
    delete[] tri_nab_tmp;
#endif
    PROFILE_STOP_L2( "FaceList::add_node" );
    return N_tri_new;
}


/********************************************************************
 * Function to check the face data                                   *
 ********************************************************************/
template<int NDIM>
static inline bool compare_index( const int *i1, const int *i2 );
template<>
inline bool compare_index<2>( const int *i1, const int *i2 )
{
    return ( i1[0] == i2[0] || i1[0] == i2[1] ) && ( i1[1] == i2[0] || i1[1] == i2[1] );
}
template<>
inline bool compare_index<3>( const int *i1, const int *i2 )
{
    return ( i1[0] == i2[0] || i1[0] == i2[1] || i1[0] == i2[2] ) &&
           ( i1[1] == i2[0] || i1[1] == i2[1] || i1[1] == i2[2] ) &&
           ( i1[2] == i2[0] || i1[2] == i2[1] || i1[2] == i2[2] );
}
template<int NDIM, class TYPE, class ETYPE>
void DelaunayTessellation::FaceList<NDIM, TYPE, ETYPE>::check_data()
{
    for ( int i = 0; i < N_face; i++ ) {
        int Nf = 0;
        for ( int j = 0; j < N_face; j++ ) {
            bool c1 = data[i].tri_id == data[j].tri_id && data[i].face_id == data[j].face_id;
            bool c2 = compare_index<NDIM>( data[i].index, data[j].index );
            if ( c1 && c2 )
                Nf++;
            else if ( c1 || c2 )
                ERROR( "Internal error" );
        }
        ASSERT( Nf == 1 );
    }
    if ( NDIM == 2 ) {
        // No node should be on the face list more than twice
        auto tmp = new unsigned char[Nx];
        memset( tmp, 0, Nx );
        for ( int i = 0; i < N_face; i++ ) {
            tmp[data[i].index[0]]++;
            tmp[data[i].index[1]]++;
        }
        bool pass = true;
        for ( int i = 0; i < Nx; i++ )
            pass = pass && tmp[i] <= 2;
        delete[] tmp;
        if ( !pass )
            ERROR( "Internal error" );
    }
}


/********************************************************************
 * This function updates the faces on the convex hull                *
 * Note:  Special attention must be paid to triangle indicies that   *
 * are in both the old and new lists.                                *
 ********************************************************************/
template<int NDIM, class TYPE, class ETYPE>
void DelaunayTessellation::FaceList<NDIM, TYPE, ETYPE>::update_face( const int N,
                                                                     const int old_tid[],
                                                                     const int old_fid[],
                                                                     const int new_tid[],
                                                                     const int new_fid[],
                                                                     const int tri_in[] )
{
    PROFILE_START_L3( "update_face" );
    // Check the inputs
    bool valid_inputs = N <= 16;
    for ( int i = 0; i < N; i++ ) {
        for ( int k = 0; k <= NDIM; k++ ) {
            if ( tri_in[k + new_tid[i] * ( NDIM + 1 )] == -1 )
                valid_inputs = false;
        }
    }
    if ( !valid_inputs )
        throw std::logic_error( "Internal error: DelaunayTessellation::FaceList::update_face" );
    // First we want to check for entries that appear in both the old and new lists
    size_t index_old[16];    // List of the unique indicies of the faces in the old list
    size_t index_new[16];    // List of the unique indicies of the faces in the new list
    size_t index_update[16]; // List of the unique indicies of the faces in both lists
    for ( int i = 0; i < N; i++ ) {
        index_old[i] = get_face_index( old_fid[i], old_tid[i] );
        index_new[i] = get_face_index( new_fid[i], new_tid[i] );
    }
    int N2       = N;
    int N_update = 0;
    for ( int i = N2 - 1; i >= 0; i-- ) {
        for ( int j = 0; j < N2; j++ ) {
            if ( index_old[i] == index_new[j] ) {
                // This face occurs in both lists, remove it
                index_update[N_update] = index_old[i];
                index_old[i]           = index_old[N2 - 1];
                index_new[j]           = index_new[N2 - 1];
                N_update++;
                N2--;
                break;
            }
        }
    }
    // Update the faces whose triangle number and face number did not change (the verticies might
    // have changed)
    for ( int i = 0; i < N_update; i++ ) {
        // Find the face
        unsigned int key = get_hash_key( index_update[i] );
        auto tri_id      = static_cast<int>( index_update[i] / ( NDIM + 1 ) );
        auto face_id     = static_cast<int>( index_update[i] % ( NDIM + 1 ) );
        int j            = hash_table[key];
        while ( j != -1 ) {
            if ( data[j].tri_id == tri_id && data[j].face_id == face_id )
                break;
            j = data[j].next;
        }
        ASSERT( j != -1 ); // Check that the face was found
        // Update the current face (note: we do not need to update the hash table data, since the
        // key did not change)
        int j1 = 0;
        for ( int j2 = 0; j2 <= NDIM; j2++ ) {
            if ( j2 == face_id )
                continue;
            int k             = tri_in[j2 + tri_id * ( NDIM + 1 )];
            data[j].index[j1] = k;
            for ( int d = 0; d < NDIM; d++ )
                data[j].x[j1][d] = x0[d + k * NDIM];
            j1++;
        }
    }
    // Replace each face that needs updating
    for ( int i = 0; i < N2; i++ ) {
        // Find the old face
        unsigned int key = get_hash_key( index_old[i] );
        auto old_tri     = static_cast<int>( index_old[i] / ( NDIM + 1 ) );
        auto old_face    = static_cast<int>( index_old[i] % ( NDIM + 1 ) );
        int j            = hash_table[key];
        while ( j != -1 ) {
            if ( data[j].tri_id == old_tri && data[j].face_id == old_face )
                break;
            j = data[j].next;
        }
        ASSERT( j != -1 ); // Check that the face was found
        // Replace the current face with the new face
        auto new_tri    = static_cast<int>( index_new[i] / ( NDIM + 1 ) );
        auto new_face   = static_cast<int>( index_new[i] % ( NDIM + 1 ) );
        data[j].tri_id  = new_tri;
        data[j].face_id = new_face;
        int j1          = 0;
        for ( int j2 = 0; j2 <= NDIM; j2++ ) {
            if ( j2 == new_face )
                continue;
            int k             = tri_in[j2 + new_tri * ( NDIM + 1 )];
            data[j].index[j1] = k;
            for ( int d = 0; d < NDIM; d++ )
                data[j].x[j1][d] = x0[d + k * NDIM];
            j1++;
        }
        // Update the hash table
        if ( data[j].prev >= 0 ) {
            int k = data[j].prev;
            ASSERT( k >= 0 );
            data[k].next = data[j].next;
            if ( data[j].next != -1 )
                data[data[j].next].prev = k;
        } else {
            int k           = data[j].next;
            hash_table[key] = k;
            if ( k != -1 )
                data[k].prev = -static_cast<int>( key + 1 );
        }
        key             = get_hash_key( index_new[i] );
        int k           = hash_table[key];
        hash_table[key] = j;
        data[j].prev    = -static_cast<int>( key + 1 );
        data[j].next    = k;
        if ( k != -1 )
            data[k].prev = j;
    }
#if DEBUG_CHECK == 2
    check_data();
#endif
    PROFILE_STOP_L3( "update_face" );
}


/********************************************************************
 * This function calculates the normal for the face of a simplex     *
 * Note: the normal is not normalized and may not be in the outward  *
 *       direction, but may be computed using exact math             *
 *                                                                   *
 * Compute the normal using the n-dimensional cross product:         *
 *    |     e1          e2      ...      en     |                    *
 *    | (A1-A2).e1  (A1-A2).e2  ...  (A1-A2).en |                    *
 *    | (A1-A3).e1  (A1-A3).e2  ...  (A1-A3).en |                    *
 *    |     ...         ...     ...      ...    |                    *
 *    | (A1-An).e1  (A1-An).e2  ...  (A1-An).en |                    *
 *                                                                   *
 ********************************************************************/
template<class ETYPE>
inline void calc_surface_normal( const ETYPE x[1][1], ETYPE norm[1] )
{
    norm[0] = ETYPE( 1 );
}
template<class ETYPE>
inline void calc_surface_normal( const ETYPE x[2][2], ETYPE norm[2] )
{
    norm[0] = x[0][1] - x[1][1];
    norm[1] = x[1][0] - x[0][0];
}
template<class ETYPE>
inline void calc_surface_normal( const ETYPE x[3][3], ETYPE norm[3] )
{
    norm[0] = ( x[0][1] - x[1][1] ) * ( x[0][2] - x[2][2] ) -
              ( x[0][2] - x[1][2] ) * ( x[0][1] - x[2][1] );
    norm[1] = ( x[0][2] - x[1][2] ) * ( x[0][0] - x[2][0] ) -
              ( x[0][0] - x[1][0] ) * ( x[0][2] - x[2][2] );
    norm[2] = ( x[0][0] - x[1][0] ) * ( x[0][1] - x[2][1] ) -
              ( x[0][1] - x[1][1] ) * ( x[0][0] - x[2][0] );
}


/********************************************************************
 * This function calculates distance between a hyperplane and a      *
 *   point.  The distance can be determined from  D_i=n·(x_0-x_i),   *
 *   while the sign can be determined by comparing to a known point  *
 *   inside the convex hull (xc).                                    *
 * Note: this uses a mixture of exact and inexact math, but should   *
 *   be accurate.  The exact math portion requires N^D precision.    *
 ********************************************************************/
template<int NDIM, class TYPE, class ETYPE>
double
DelaunayTessellation::FaceList<NDIM, TYPE, ETYPE>::calc_surface_distance( const TYPE x[NDIM][NDIM],
                                                                          const TYPE xi[] ) const
{
    // First compute the normal
    // Note: the normal is not normalized and may point inward
    ETYPE norm[NDIM], x2[NDIM][NDIM];
    for ( int i = 0; i < NDIM; i++ ) {
        for ( int j = 0; j < NDIM; j++ ) {
            x2[i][j] = ETYPE( x[i][j] );
        }
    }
    calc_surface_normal( x2, norm );
    // Compute the distance from the surface
    ETYPE dist( 0 );
    long double dot = 0.0;
    long double tmp = 0.0;
    for ( int i = 0; i < NDIM; i++ ) {
        long double norm2 = get_double( norm[i] );
        dot += norm2 * ( x[0][i] - xc[i] );
        dist += norm[i] * ETYPE( xi[i] - x[0][i] );
        tmp += norm2 * norm2;
    }
    double sign = ( dot < 0 ) ? -1.0 : 1.0;
    return sign * get_double( dist ) / sqrt( tmp );
}


/********************************************************************
 * This function determines if a point is outside a triangle on the  *
 *    convex hull.  This is done by using the calculation for the    *
 *    distance to the plane and the point inside the hull in the     *
 *    same mannor as calc_surface_distance, but should be more       *
 *    efficient.                                                     *
 * Note: this uses a mixture of exact and inexact math, but should   *
 *   be accurate.  The exact math portion requires N^D precision.    *
 ********************************************************************/
template<>
bool DelaunayTessellation::FaceList<2, int, int>::outside_triangle( const int x[2][2],
                                                                    const int xi[] ) const
{
    int nx     = x[0][1] - x[1][1];
    int ny     = x[1][0] - x[0][0];
    int dist   = nx * ( xi[0] - x[0][0] ) + ny * ( xi[1] - x[0][1] );
    double dot = nx * ( x[0][0] - xc[0] ) + ny * ( x[0][1] - xc[1] );
    return dot * dist > 0;
}
template<>
bool DelaunayTessellation::FaceList<3, int, int>::outside_triangle( const int x[3][3],
                                                                    const int xi[] ) const
{
    int nx = ( x[0][1] - x[1][1] ) * ( x[0][2] - x[2][2] ) -
             ( x[0][2] - x[1][2] ) * ( x[0][1] - x[2][1] );
    int ny = ( x[0][2] - x[1][2] ) * ( x[0][0] - x[2][0] ) -
             ( x[0][0] - x[1][0] ) * ( x[0][2] - x[2][2] );
    int nz = ( x[0][0] - x[1][0] ) * ( x[0][1] - x[2][1] ) -
             ( x[0][1] - x[1][1] ) * ( x[0][0] - x[2][0] );
    int dist   = nx * ( xi[0] - x[0][0] ) + ny * ( xi[1] - x[0][1] ) + nz * ( xi[2] - x[0][2] );
    double dot = nx * ( x[0][0] - xc[0] ) + ny * ( x[0][1] - xc[1] ) + nz * ( x[0][2] - xc[2] );
    return dot * dist > 0;
}
template<int NDIM, class TYPE, class ETYPE>
bool DelaunayTessellation::FaceList<NDIM, TYPE, ETYPE>::outside_triangle( const TYPE x[NDIM][NDIM],
                                                                          const TYPE xi[] ) const
{
    // First compute the normal
    // Note: the normal is not normalized and may point inward
    ETYPE norm[NDIM], x2[NDIM][NDIM];
    for ( int i = 0; i < NDIM; i++ ) {
        for ( int j = 0; j < NDIM; j++ ) {
            x2[i][j] = ETYPE( x[i][j] );
        }
    }
    calc_surface_normal( x2, norm );
    // Compute the distance from the surface
    ETYPE dist( 0 );
    long double dot = 0.0;
    for ( int i = 0; i < NDIM; i++ ) {
        long double norm2 = get_double( norm[i] );
        dot += norm2 * ( x[0][i] - xc[i] );
        dist += norm[i] * ETYPE( xi[i] - x[0][i] );
    }
    return dot * get_double( dist ) > 0;
}


/********************************************************************
 * FaceList::face_data_struct functions                              *
 ********************************************************************/
template<int NDIM, class TYPE, class ETYPE>
void DelaunayTessellation::FaceList<NDIM, TYPE, ETYPE>::face_data_struct::reset()
{
    prev    = -1;
    next    = -1;
    tri_id  = -1;
    face_id = -1;
    for ( int i = 0; i < NDIM; i++ ) {
        index[i] = -1;
        for ( int j = 0; j < NDIM; j++ )
            x[i][j] = 0;
    }
}


/********************************************************************
 * Function to check the current triangles to see if any are invalid *
 ********************************************************************/
template<int NDIM, class TYPE, class ETYPE>
bool DelaunayTessellation::check_current_triangles( int N,
                                                    const TYPE x[],
                                                    size_t N_tri,
                                                    const int tri[],
                                                    const int tri_nab[],
                                                    const std::vector<size_t> &unused,
                                                    double TOL_VOL )
{
    if ( N_tri <= 0 )
        return false;
    NULL_USE( TOL_VOL );
    PROFILE_START_L2( "check_current_triangles" );
    bool pass = true;
    for ( size_t i = 0; i < N_tri; i++ ) {
        // First check if we are dealing with an unused triangle
        bool used = true;
        for ( auto &elem : unused ) {
            if ( i == elem )
                used = false;
        }
        if ( !used )
            break;
        // Check if the triangle is valid
        TYPE x2[NDIM * ( NDIM + 1 )];
        for ( int d = 0; d <= NDIM; d++ ) {
            int k = tri[d + i * ( NDIM + 1 )];
            if ( k < 0 || k >= N ) {
                printp( "Error with internal structures (invalid triangle index)\n" );
                pass = false;
                break;
            }
            for ( int d2 = 0; d2 < NDIM; d2++ )
                x2[d2 + d * NDIM] = x[d2 + k * NDIM];
        }
        if ( !pass ) {
            break;
        }
        double vol = fabs( calc_volume<NDIM, TYPE, ETYPE>( x2 ) );
        if ( vol < 0.0 ) {
            // The volume of each triangle must be strickly > 0
            printp( "Invalid triangle volume detected\n" );
            pass = false;
        }
        if ( !pass ) {
            break;
        }
        // Check if the triangle neighbors are valid
        for ( int d = 0; d <= NDIM; d++ ) {
            int k = tri_nab[d + i * ( NDIM + 1 )];
            if ( k == -1 )
                continue;
            if ( k < 0 || k >= (int) N_tri )
                pass = false;
            for ( auto &elem : unused ) {
                if ( k == static_cast<int>( elem ) )
                    pass = false;
            }
            for ( int d2 = 0; d2 <= NDIM; d2++ ) {
                if ( tri_nab[d2 + k * ( NDIM + 1 )] == k )
                    pass = false;
            }
            bool shared = false;
            for ( int d2 = 0; d2 <= NDIM; d2++ ) {
                if ( tri_nab[d2 + k * ( NDIM + 1 )] == (int) i )
                    shared = true;
            }
            if ( !shared )
                pass = false;
        }
        if ( !pass ) {
            printp( "Error with internal structures (invalid triangle neighbor)\n" );
            break;
        }
    }
    PROFILE_STOP_L2( "check_current_triangles" );
    return pass;
}


// Subroutine to check if we conserve the neighbor triangles for two sets
static bool conserved_neighbors( const int N1, const int list1[], const int N2, const int list2[] )
{
    if ( N1 > 16 || N2 > 16 ) {
        printp( "Need to change internal data structures\n" );
        return false;
    }
    int tmp[16];
    for ( int i = 0; i < N1; i++ )
        tmp[i] = list1[i];
    for ( int j = 0; j < N2; j++ ) {
        if ( list2[j] < 0 )
            continue;
        bool found = false;
        for ( int i = 0; i < N1; i++ ) {
            if ( tmp[i] < 0 )
                continue;
            if ( tmp[i] == list2[j] ) {
                found  = true;
                tmp[i] = -1;
                break;
            }
        }
        if ( !found )
            return false;
    }
    for ( int i = 0; i < N1; i++ ) {
        if ( tmp[i] >= 0 )
            return false;
    }
    return true;
}


// Subroutine to check if two given triangles share a face
bool are_tri_neighbors( const int ndim, const int tri1[], const int tri2[], int *f1, int *f2 )
{
    // Lets check each node in tri1 to see if it is in tri2
    *f1 = -1;
    for ( int i = 0; i <= ndim; i++ ) {
        int tmp    = tri1[i];
        bool found = false;
        for ( int j = 0; j <= ndim; j++ ) {
            if ( tmp == tri2[j] ) {
                found = true;
                break;
            }
        }
        if ( !found ) {
            if ( *f1 != -1 ) {
                // There are 2 or more nodes not in common
                return false;
            } else {
                *f1 = i;
            }
        }
    }
    if ( *f1 == -1 ) {
        // All nodes are in common
        return false;
    }
    // Now we need to find the missing node in tri2
    *f2 = -1;
    for ( int i = 0; i <= ndim; i++ ) {
        int tmp    = tri2[i];
        bool found = false;
        for ( int j = 0; j <= ndim; j++ ) {
            if ( tmp == tri1[j] ) {
                found = true;
                break;
            }
        }
        if ( !found ) {
            *f2 = i;
            break;
        }
    }
    return true;
}


/********************************************************************
 * Primary interfaces                                                *
 ********************************************************************/
template<class TYPE>
uint64_t maxabs( const int N, const TYPE *x )
{
    uint64_t x_max = 0;
    for ( int i = 0; i < N; i++ ) {
        TYPE x2 = x[i] >= 0 ? x[i] : -x[i];
        x_max   = std::max<uint64_t>( static_cast<uint64_t>( x2 ), x_max );
    }
    return x_max;
}
int DelaunayTessellation::create_tessellation(
    const int ndim, const int N, const double x[], int *tri[], int *tri_nab[] )
{
    PROFILE_START( "create_tessellation (double)", 2 );
    int N_tri = -1;
    if ( ndim == 2 )
        N_tri = create_tessellation<2, double, long double>( N, x, tri, tri_nab );
    else if ( ndim == 3 )
        N_tri = create_tessellation<3, double, long double>( N, x, tri, tri_nab );
    PROFILE_STOP( "create_tessellation (double)", 2 );
    return N_tri;
}
int DelaunayTessellation::create_tessellation(
    const int ndim, const int N, const short int x[], int *tri[], int *tri_nab[] )
{
    PROFILE_START( "create_tessellation (short int)", 2 );
    auto y = new int[N * ndim];
    for ( int i = 0; i < ndim * N; i++ )
        y[i] = x[i];
    int N_tri = -1;
    if ( ndim == 2 ) {
        N_tri = create_tessellation<2, int, int>( N, y, tri, tri_nab );
    } else if ( ndim == 3 ) {
        N_tri = create_tessellation<3, int, int64_t>( N, y, tri, tri_nab );
    }
    delete[] y;
    PROFILE_STOP( "create_tessellation (short int)", 2 );
    return N_tri;
}
int DelaunayTessellation::create_tessellation(
    const int ndim, const int N, const int x[], int *tri[], int *tri_nab[] )
{
    PROFILE_START( "create_tessellation (int)", 2 );
    // Check that all values of x are < 2^30
    const uint64_t x_max = maxabs( ndim * N, x );
    const uint64_t XMAX  = ( (uint64_t) 1 ) << 30;
    if ( x_max >= XMAX ) {
        printp( "To be stable, all vaues of x must < 2^30 for type int" );
        PROFILE_STOP2( "create_tessellation (int)", 2 );
        return -8;
    }
    int N_tri = -1;
    if ( ndim == 2 ) {
        if ( x_max < 32768 ) {
            N_tri = create_tessellation<2, int, int>( N, x, tri, tri_nab );
        } else {
            N_tri = create_tessellation<2, int, int64_t>( N, x, tri, tri_nab );
        }
    } else if ( ndim == 3 ) {
        if ( x_max < 768 ) {
            N_tri = create_tessellation<3, int, int>( N, x, tri, tri_nab );
        } else if ( x_max < 1048576 ) {
            N_tri = create_tessellation<3, int, int64_t>( N, x, tri, tri_nab );
        } else {
#ifdef DISABLE_EXTENDED
            ERROR( "create_tessellation (int) is disabled for large ints" );
            return -8;
#else
            N_tri = create_tessellation<3, int, int128_t>( N, x, tri, tri_nab );
#endif
        }
    }
    PROFILE_STOP( "create_tessellation (int)", 2 );
    return N_tri;
}
int DelaunayTessellation::create_tessellation(
    const int ndim, const int N, const int64_t x[], int *tri[], int *tri_nab[] )
{
    PROFILE_START( "create_tessellation (int64_t)", 2 );
    // Check that all values of x are < 2^62
    const uint64_t x_max = maxabs( ndim * N, x );
    const uint64_t XMAX  = ( (uint64_t) 1 ) << 62;
    if ( x_max >= XMAX ) {
        printp( "To be stable, all vaues of x must < 2^62 for type int\n" );
        PROFILE_STOP2( "create_tessellation (int64_t)", 2 );
        return -8;
    }
    int N_tri = -1;
#ifdef DISABLE_EXTENDED
    ERROR( "create_tessellation (int64_t) is disabled" );
    NULL_USE( tri );
    NULL_USE( tri_nab );
    return -8;
#else
    if ( ndim == 2 ) {
        N_tri = create_tessellation<2, int64_t, int128_t>( N, x, tri, tri_nab );
    } else if ( ndim == 3 ) {
        N_tri = create_tessellation<3, int64_t, int256_t>( N, x, tri, tri_nab );
    }
#endif
    PROFILE_STOP( "create_tessellation (int64_t)", 2 );
    return N_tri;
}
double DelaunayTessellation::calc_volume( int ndim, const double x[] )
{
    double vol = 0.0;
    if ( ndim == 1 ) {
        vol = x[1] - x[0];
    } else if ( ndim == 2 ) {
        vol = calc_volume<2, double, long double>( x );
    } else if ( ndim == 3 ) {
        vol = calc_volume<3, double, long double>( x );
    } else if ( ndim == 4 ) {
        vol = calc_volume<4, double, long double>( x );
    } else {
        throw std::logic_error( "Unsupported dimension" );
    }
    return vol;
}
int DelaunayTessellation::test_in_circumsphere( int ndim,
                                                const double x[],
                                                const double xi[],
                                                const double TOL_VOL )
{
    int test = -2;
    if ( ndim == 1 ) {
        test = test_in_circumsphere<1, double, long double>( x, xi, TOL_VOL );
    } else if ( ndim == 2 ) {
        test = test_in_circumsphere<2, double, long double>( x, xi, TOL_VOL );
    } else if ( ndim == 3 ) {
        test = test_in_circumsphere<3, double, long double>( x, xi, TOL_VOL );
    } else {
        throw std::logic_error( "Unsupported dimension" );
    }
    return test;
}
int DelaunayTessellation::test_in_circumsphere( int ndim,
                                                const int x[],
                                                const int xi[],
                                                const double TOL_VOL )
{
    int test = -2;
    if ( ndim == 1 ) {
        test = test_in_circumsphere<1, int, int>( x, xi, TOL_VOL );
    } else if ( ndim == 2 ) {
        test = test_in_circumsphere<2, int, int64_t>( x, xi, TOL_VOL );
    } else if ( ndim == 3 ) {
        const uint64_t x_max = std::max( maxabs( ndim * ( ndim + 1 ), x ), maxabs( ndim, xi ) );
        if ( x_max < 768 ) {
            test = test_in_circumsphere<3, int, int>( x, xi, TOL_VOL );
        } else if ( x_max < 1048576 ) {
            test = test_in_circumsphere<3, int, int64_t>( x, xi, TOL_VOL );
        } else {
#ifdef DISABLE_EXTENDED
            ERROR( "create_tessellation (int) is disabled for large ints" );
#else
            test = test_in_circumsphere<3, int, int128_t>( x, xi, TOL_VOL );
#endif
        }
    } else {
        throw std::logic_error( "Unsupported dimension" );
    }
    return test;
}
void DelaunayTessellation::get_circumsphere( const int ndim,
                                             const double x[],
                                             double &R,
                                             double *center )
{
    if ( ndim == 1 ) {
        get_circumsphere<1, double>( x, R, center );
    } else if ( ndim == 2 ) {
        get_circumsphere<2, double>( x, R, center );
    } else if ( ndim == 3 ) {
        get_circumsphere<3, double>( x, R, center );
    } else {
        throw std::logic_error( "Unsupported dimension" );
    }
}
void DelaunayTessellation::get_circumsphere( const int ndim,
                                             const int x[],
                                             double &R,
                                             double *center )
{
    if ( ndim == 1 ) {
        get_circumsphere<1, int>( x, R, center );
    } else if ( ndim == 2 ) {
        get_circumsphere<2, int>( x, R, center );
    } else if ( ndim == 3 ) {
        get_circumsphere<3, int>( x, R, center );
    } else {
        throw std::logic_error( "Unsupported dimension" );
    }
}
