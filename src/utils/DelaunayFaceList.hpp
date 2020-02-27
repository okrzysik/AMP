#ifndef included_AMP_DelaunayFaceList_hpp
#define included_AMP_DelaunayFaceList_hpp

#include "AMP/utils/DelaunayFaceList.h"
#include "AMP/utils/DelaunayHelpers.h"
#include "AMP/utils/Utilities.h"

#include <stdint.h>
#include <stdlib.h>
#include <vector>


// Macros to define 3 levels of profilers
//#if ( defined( DEBUG ) || defined( _DEBUG ) ) && !defined( NDEBUG )
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


#define DEBUG_CHECK 0 // Flag to enable extra checks (1: some additional cost, 2: very expensive)


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


namespace AMP::DelaunayTessellation {


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
static constexpr double inv_factorial( int N )
{
    double x = 1;
    for ( int i = 2; i <= N; i++ )
        x *= i;
    return 1.0 / x;
}
template<int NDIM, class TYPE, class ETYPE>
double calc_volume( const std::array<TYPE, NDIM> x[] )
{
    if constexpr ( NDIM == 1 )
        return static_cast<double>( x[1][0] - x[0][0] );
    ETYPE M[NDIM * NDIM];
    for ( int d = 0; d < NDIM; d++ ) {
        ETYPE tmp( x[NDIM][d] );
        for ( int j = 0; j < NDIM; j++ )
            M[d + j * NDIM] = ETYPE( x[j][d] ) - tmp;
    }
    constexpr double C = inv_factorial( NDIM );
    return C * get_double( DelaunayHelpers<NDIM>::det( M ) );
}
/*template<int NDIM, class TYPE, class ETYPE>
double calc_volume( const TYPE x[] )
{
    // This will be removed
    std::array<TYPE, NDIM> x2[NDIM + 1];
    for ( int i = 0; i <= NDIM; i++ ) {
        for ( int d = 0; d < NDIM; d++ )
            x2[i][d] = x[d + i * NDIM];
    }
    return calc_volume<NDIM, TYPE, ETYPE>( x2 );
}*/


/********************************************************************
 * Check if we conserve the neighbor triangles for two sets          *
 ********************************************************************/
static bool conserved_neighbors( const int N1, const int list1[], const int N2, const int list2[] )
{
    if ( N1 > 16 || N2 > 16 ) {
        printf( "Need to change internal data structures\n" );
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


/********************************************************************
 * Check if two given triangles share a face                         *
 ********************************************************************/
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
 * This function performs a set of test taht checks that each        *
 *    triangle is valid and that the triangle neighbors are valid.   *
 * Note: This can be an expensive test and should only be used when  *
 *    debugging.                                                     *
 * Inputs:                                                           *
 *    N         The number of verticies                              *
 *    x         The coordinates of the verticies (NDIMxN)            *
 *    N_tri     The number of triangles                              *
 *    tri       The current list of triangles (NDIM+1xN_tri)         *
 *    tri_nab   The current list of triangle neighbors (NDIM+1xN_tri)*
 *    N_unused  The number of unused coordinates                     *
 *    TOL_VOL   The tolerance to useto check if a simplex is valid   *
 ********************************************************************/
template<int NDIM, class TYPE, class ETYPE>
bool check_current_triangles( int N,
                              const std::array<TYPE, NDIM> x[],
                              size_t N_tri,
                              const std::array<int, NDIM + 1> tri[],
                              const std::array<int, NDIM + 1> tri_nab[],
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
        std::array<TYPE, NDIM> x2[NDIM + 1];
        for ( int d = 0; d <= NDIM; d++ ) {
            int k = tri[i][d];
            if ( k < 0 || k >= N ) {
                printf( "Error with internal structures (invalid triangle index)\n" );
                pass = false;
                break;
            }
            x2[d] = x[k];
        }
        if ( !pass ) {
            break;
        }
        double vol = fabs( calc_volume<NDIM, TYPE, ETYPE>( x2 ) );
        if ( vol < 0.0 ) {
            // The volume of each triangle must be strickly > 0
            printf( "Invalid triangle volume detected\n" );
            pass = false;
        }
        if ( !pass ) {
            break;
        }
        // Check if the triangle neighbors are valid
        for ( int d = 0; d <= NDIM; d++ ) {
            int k = tri_nab[i][d];
            if ( k == -1 )
                continue;
            if ( k < 0 || k >= (int) N_tri )
                pass = false;
            for ( auto &elem : unused ) {
                if ( k == static_cast<int>( elem ) )
                    pass = false;
            }
            for ( int d2 = 0; d2 <= NDIM; d2++ ) {
                if ( tri_nab[k][d2] == k )
                    pass = false;
            }
            bool shared = false;
            for ( int d2 = 0; d2 <= NDIM; d2++ ) {
                if ( tri_nab[k][d2] == (int) i )
                    shared = true;
            }
            if ( !shared )
                pass = false;
        }
        if ( !pass ) {
            printf( "Error with internal structures (invalid triangle neighbor)\n" );
            break;
        }
    }
    PROFILE_STOP_L2( "check_current_triangles" );
    return pass;
}


/********************************************************************
 * FaceList constructor                                              *
 ********************************************************************/
static inline unsigned int get_hash_key( size_t index )
{
    auto key = static_cast<unsigned int>( index );
    key *= 0x9E3779B9; // 2^32*0.5*(sqrt(5)-1)
    key >>= 22;        // mod(key,1024)
    return key;
}
template<int NDIM, class TYPE, class ETYPE>
DelaunayTessellation::FaceList<NDIM, TYPE, ETYPE>::FaceList( const int Nx_in,
                                                             const Point *x_in,
                                                             const int tri0_id,
                                                             const std::array<int, NDIM + 1> &tri0,
                                                             const TYPE TOL_VOL )
    : Nx( Nx_in ), x0( x_in ), TOL_vol( TOL_VOL )
{
    AMP_ASSERT( sizeof( unsigned int ) ==
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
            data[i].x[j1]     = x0[k];
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
    AMP::Utilities::quicksort( N_delete, ids );
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
int DelaunayTessellation::FaceList<NDIM, TYPE, ETYPE>::add_node(
    const int node_id,
    std::vector<size_t> &unused,
    size_t &N_tri,
    unsigned int *new_tri_id,
    std::array<int, NDIM + 1> *new_tri,
    std::array<int, NDIM + 1> *new_tri_nab,
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
    Point x2[NDIM + 1];
    for ( int i = 0; i < N_face; i++ ) {
        // Get the triangle number and face id
        int tri_num  = data[i].tri_id;
        int face_num = data[i].face_id;
// We will need to for a new triangle with each face where the distance
// from the face to the point is > 0 and the volume is >= TOL_vol
#if 0
            bool pass = calc_surface_distance(data[i].x,&x0[node_id*NDIM]) > 0;
#else
        bool pass = outside_triangle( data[i].x, x0[node_id] );
#endif
        if ( !pass )
            continue;
        for ( int j = 0; j < NDIM; j++ )
            x2[j] = data[i].x[j];
        x2[NDIM]      = x0[node_id];
        double volume = fabs( calc_volume<NDIM, TYPE, ETYPE>( x2 ) );
        if ( volume <= TOL_vol ) {
            continue;
        }
        // Create a new triangle consisting of the current point and the current face
        Triangle tri, nab;
        for ( int j = 0; j < NDIM; j++ )
            tri[j] = data[i].index[j];
        tri[NDIM] = node_id;
        for ( int j = 0; j < NDIM; j++ )
            nab[j] = -1;
        nab[NDIM] = tri_num;
        int id;
        if ( !unused.empty() ) {
            id = unused.back();
            unused.pop_back();
        } else {
            id = N_tri;
            N_tri++;
        }
        new_tri_id[N_tri_new]  = id;
        new_tri[N_tri_new]     = tri;
        new_tri_nab[N_tri_new] = nab;
        neighbor[N_tri_new]    = tri_num;
        face_id[N_tri_new]     = face_num;
        ids[N_tri_new]         = i;
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
            int k  = new_tri[i][j1];
            x2[j1] = x0[k];
        }
        double volume = calc_volume<NDIM, TYPE, ETYPE>( x2 );
        if ( fabs( volume ) <= TOL_vol ) {
            // The triangle is invalid
            delete[] ids;
            PROFILE_STOP2_L2( "FaceList::add_node" );
            return -1;
        } else if ( volume < 0 ) {
            // The volume is negitive, swap the order of the last two points
            std::swap( new_tri[i][NDIM - 1], new_tri[i][NDIM] );
            std::swap( new_tri_nab[i][NDIM - 1], new_tri_nab[i][NDIM] );
        }
        for ( int j1 = 0; j1 <= NDIM; j1++ ) {
            int k  = new_tri[i][j1];
            x2[j1] = x0[k];
        }
        volume = calc_volume<NDIM, TYPE, ETYPE>( x2 );
        if ( volume <= 0.0 )
            printf( "Warning: volume is still negitive\n" );
    }
    // Check if the current triangle shares a face with one of the other new triangles
    for ( int i1 = 0; i1 < N_tri_new; i1++ ) {
        for ( int i2 = i1 + 1; i2 < N_tri_new; i2++ ) {
            int f1, f2;
            bool is_neighbor =
                are_tri_neighbors( NDIM, new_tri[i1].data(), new_tri[i2].data(), &f1, &f2 );
            if ( is_neighbor ) {
#if DEBUG_CHECK > 0
                AMP_ASSERT( new_tri_nab[i1][f1] == -1 );
                AMP_ASSERT( new_tri_nab[i2][f2] == -1 );
#endif
                new_tri_nab[i1][f1] = new_tri_id[i2];
                new_tri_nab[i2][f2] = new_tri_id[i1];
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
    for ( int i = 0; i < N_tri_new; i++ ) {
        for ( int d = 0; d <= NDIM; d++ ) {
            if ( new_tri_nab[i][d] == -1 )
                N_add++;
        }
    }
    if constexpr ( NDIM == 2 ) {
        AMP_ASSERT( N_add == 2 );
    } else if constexpr ( NDIM == 3 ) {
        AMP_ASSERT( N_add >= 3 );
        AMP_ASSERT( N_add <= 2 * N_tri_new + 1 );
    }
    for ( int i = 0; i < N_tri_new; i++ ) {
        for ( int j = 0; j <= NDIM; j++ ) {
            if ( new_tri_nab[i][j] == -1 ) {
                // Add the face
                data[N_face].tri_id  = new_tri_id[i];
                data[N_face].face_id = j;
                int j1               = 0;
                for ( int j2 = 0; j2 <= NDIM; j2++ ) {
                    if ( j2 == j )
                        continue;
                    int k                  = new_tri[i][j2];
                    data[N_face].index[j1] = k;
                    data[N_face].x[j1]     = x0[k];
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
    auto tri_nab_tmp = new Triangle[N_tri_new];
    memcpy( tri_nab_tmp, new_tri_nab, N_tri_new * sizeof( Triangle ) );
    for ( int i = 0; i < N_tri_new; i++ ) {
        for ( int d = 0; d <= NDIM; d++ ) {
            int k = -1;
            for ( int j = 0; j < N_tri_new; j++ )
                if ( tri_nab_tmp[i][d] == (int) new_tri_id[j] )
                    k = j;
            tri_nab_tmp[i][d] = k;
        }
    }
    bool all_valid = check_current_triangles<NDIM, TYPE, ETYPE>(
        Nx, x0, N_tri_new, new_tri, tri_nab_tmp, std::vector<size_t>(), 0 );
    if ( !all_valid )
        printf( "Warning: new triangle neighbors are inconsistent\n" );
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
                AMP_ERROR( "Internal error" );
        }
        AMP_ASSERT( Nf == 1 );
    }
    if constexpr ( NDIM == 2 ) {
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
            AMP_ERROR( "Internal error" );
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
                                                                     const Triangle *tri_in )
{
    PROFILE_START_L3( "update_face" );
    // Check the inputs
    bool valid_inputs = N <= 16;
    for ( int i = 0; i < N; i++ ) {
        for ( int k = 0; k <= NDIM; k++ ) {
            if ( tri_in[new_tid[i]][k] == -1 )
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
        AMP_ASSERT( j != -1 ); // Check that the face was found
        // Update the current face (note: we do not need to update the hash table data, since the
        // key did not change)
        int j1 = 0;
        for ( int j2 = 0; j2 <= NDIM; j2++ ) {
            if ( j2 == face_id )
                continue;
            int k             = tri_in[tri_id][j2];
            data[j].index[j1] = k;
            data[j].x[j1]     = x0[k];
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
        AMP_ASSERT( j != -1 ); // Check that the face was found
        // Replace the current face with the new face
        auto new_tri    = static_cast<int>( index_new[i] / ( NDIM + 1 ) );
        auto new_face   = static_cast<int>( index_new[i] % ( NDIM + 1 ) );
        data[j].tri_id  = new_tri;
        data[j].face_id = new_face;
        int j1          = 0;
        for ( int j2 = 0; j2 <= NDIM; j2++ ) {
            if ( j2 == new_face )
                continue;
            int k             = tri_in[new_tri][j2];
            data[j].index[j1] = k;
            data[j].x[j1]     = x0[k];
            j1++;
        }
        // Update the hash table
        if ( data[j].prev >= 0 ) {
            int k = data[j].prev;
            AMP_ASSERT( k >= 0 );
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
double DelaunayTessellation::FaceList<NDIM, TYPE, ETYPE>::calc_surface_distance(
    const std::array<TYPE, NDIM> x[NDIM], const std::array<TYPE, NDIM> &xi ) const
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
bool DelaunayTessellation::FaceList<2, int, int>::outside_triangle(
    const std::array<int, 2> x[2], const std::array<int, 2> &xi ) const
{
    int nx     = x[0][1] - x[1][1];
    int ny     = x[1][0] - x[0][0];
    int dist   = nx * ( xi[0] - x[0][0] ) + ny * ( xi[1] - x[0][1] );
    double dot = nx * ( x[0][0] - xc[0] ) + ny * ( x[0][1] - xc[1] );
    return dot * dist > 0;
}
template<>
bool DelaunayTessellation::FaceList<3, int, int>::outside_triangle(
    const std::array<int, 3> x[3], const std::array<int, 3> &xi ) const
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
bool DelaunayTessellation::FaceList<NDIM, TYPE, ETYPE>::outside_triangle(
    const std::array<TYPE, NDIM> x[NDIM], const std::array<TYPE, NDIM> &xi ) const
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


} // namespace AMP::DelaunayTessellation

#endif
