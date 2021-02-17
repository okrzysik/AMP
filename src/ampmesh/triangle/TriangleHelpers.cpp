#include "AMP/ampmesh/triangle/TriangleHelpers.h"
#include "AMP/ampmesh/MeshGeometry.h"
#include "AMP/ampmesh/MeshParameters.h"
#include "AMP/ampmesh/MeshPoint.h"
#include "AMP/ampmesh/MeshUtilities.h"
#include "AMP/ampmesh/MultiGeometry.h"
#include "AMP/ampmesh/MultiMesh.h"
#include "AMP/ampmesh/shapes/Circle.h"
#include "AMP/ampmesh/shapes/GeometryHelpers.h"
#include "AMP/ampmesh/shapes/Sphere.h"
#include "AMP/ampmesh/triangle/TriangleMesh.h"
#include "AMP/utils/Database.h"
#include "AMP/utils/DelaunayHelpers.h"
#include "AMP/utils/DelaunayTessellation.h"
#include "AMP/utils/NearestPairSearch.h"
#include "AMP/utils/Utilities.h"
#include "AMP/utils/kdtree.h"

#include "ProfilerApp.h"

#include <algorithm>
#include <map>
#include <random>


namespace AMP {
namespace Mesh {
namespace TriangleHelpers {


// Helper function to create constexpr std::array with a single value
template<class T, std::size_t N, std::size_t... I>
static constexpr std::array<std::remove_cv_t<T>, N> to_array_impl( const T *a,
                                                                   std::index_sequence<I...> )
{
    return { { a[I]... } };
}
template<class TYPE, std::size_t N>
static constexpr std::array<TYPE, N> make_array( const TYPE &x )
{
    TYPE tmp[N] = { x };
    return to_array_impl<TYPE, N>( tmp, std::make_index_sequence<N>{} );
}


// Helper function to wrap fread
static inline void fread2( void *ptr, size_t size, size_t count, FILE *stream )
{
    size_t N = fread( ptr, size, count, stream );
    AMP_ASSERT( N == count );
}


// Helper functions to see if two points are ~ the same
template<size_t N>
static inline bool approx_equal( const std::array<double, N> &x,
                                 const std::array<double, N> &y,
                                 const std::array<double, N> &tol )
{
    if constexpr ( N == 1 )
        return fabs( x[0] - y[0] ) <= tol[0];
    else if constexpr ( N == 2 )
        return fabs( x[0] - y[0] ) <= tol[0] && fabs( x[1] - y[1] ) <= tol[1];
    else if constexpr ( N == 3 )
        return fabs( x[0] - y[0] ) <= tol[0] && fabs( x[1] - y[1] ) <= tol[1] &&
               fabs( x[2] - y[2] ) <= tol[2];
}


/****************************************************************
 * Find the first n intersections in multiple lists              *
 * This function assumes the lists are in sorted order           *
 ****************************************************************/
static int intersect_sorted(
    const int N_lists, const int size[], int64_t *list[], const int N_max, int64_t *intersection )
{
    if ( N_max <= 0 )
        return ~( (unsigned int) 0 );
    int N_int = 0;
    std::vector<int> index( N_lists );
    for ( int i = 0; i < N_lists; i++ )
        index[i] = 0;
    unsigned int current_val = list[0][0];
    bool finished            = false;
    while ( true ) {
        unsigned int min_val = 2147483647;
        bool in_intersection = true;
        for ( int i = 0; i < N_lists; i++ ) {
            if ( index[i] >= size[i] ) {
                finished = true;
                break;
            }
            while ( list[i][index[i]] < current_val ) {
                index[i]++;
                if ( index[i] >= size[i] ) {
                    finished = true;
                    break;
                }
            }
            if ( list[i][index[i]] == current_val ) {
                index[i]++;
            } else {
                in_intersection = false;
            }
            if ( index[i] < size[i] ) {
                if ( list[i][index[i]] < min_val )
                    min_val = list[i][index[i]];
            }
        }
        if ( finished )
            break;
        if ( in_intersection ) {
            intersection[N_int] = current_val;
            N_int++;
            if ( N_int >= N_max )
                break;
        }
        current_val = min_val;
    }
    return N_int;
}


/****************************************************************
 * Count the number of unique triangles                          *
 ****************************************************************/
template<size_t NDIM, bool ordered>
static uint64_t hash( const std::array<int64_t, NDIM> &x )
{
    uint64_t hash = 0;
    for ( size_t i = 0; i < NDIM; i++ ) {
        // Use hashing function: 2^64*0.5*(sqrt(5)-1)
        uint64_t z = static_cast<uint64_t>( x[i] ) * 0x9E3779B97F4A7C15;
        if constexpr ( ordered ) {
            hash = ( ( hash << 5 ) + hash ) ^ z;
        } else {
            hash = hash ^ z;
        }
    }
    return hash;
}
template<size_t NG>
size_t count( const std::vector<std::array<int64_t, NG + 1>> &tri )
{
    std::vector<uint64_t> x( tri.size(), 0 );
    for ( size_t i = 0; i < tri.size(); i++ )
        x[i] = hash<NG + 1, false>( tri[i] );
    std::sort( x.begin(), x.end() );
    x.erase( std::unique( x.begin(), x.end() ), x.end() );
    return x.size();
}


/****************************************************************
 * Read stl file                                                 *
 ****************************************************************/
size_t readSTLHeader( const std::string &filename )
{
    char header[80];
    uint32_t N;
    auto fid = fopen( filename.c_str(), "rb" );
    AMP_INSIST( fid, "Unable to open " + filename );
    fread2( header, sizeof( header ), 1, fid );
    fread2( &N, sizeof( N ), 1, fid );
    fclose( fid );
    return N;
}

std::vector<std::array<std::array<double, 3>, 3>> readSTL( const std::string &filename,
                                                           double scale )
{
    char header[80];
    uint32_t N;
    // Read the file
    auto fid = fopen( filename.c_str(), "rb" );
    AMP_INSIST( fid, "Unable to open " + filename );
    fread2( header, sizeof( header ), 1, fid );
    fread2( &N, sizeof( N ), 1, fid );
    auto tmp = new char[N * 50];
    fread2( tmp, 50, N, fid );
    fclose( fid );
    // Get a list of the local triangles based on their coordinates
    std::vector<std::array<std::array<double, 3>, 3>> tri_coord( N );
    for ( size_t i = 0; i < N; i++ ) {
        uint16_t attrib    = 0;
        float normal[3]    = { 0, 0, 0 };
        float vertex[3][3] = { { 0 } };
        memcpy( normal, &tmp[50 * i], sizeof( normal ) );
        memcpy( vertex, &tmp[50 * i + 12], sizeof( vertex ) );
        memcpy( &attrib, &tmp[50 * i + 48], sizeof( attrib ) );
        tri_coord[i][0][0] = scale * vertex[0][0];
        tri_coord[i][0][1] = scale * vertex[0][1];
        tri_coord[i][0][2] = scale * vertex[0][2];
        tri_coord[i][1][0] = scale * vertex[1][0];
        tri_coord[i][1][1] = scale * vertex[1][1];
        tri_coord[i][1][2] = scale * vertex[1][2];
        tri_coord[i][2][0] = scale * vertex[2][0];
        tri_coord[i][2][1] = scale * vertex[2][1];
        tri_coord[i][2][2] = scale * vertex[2][2];
        NULL_USE( attrib );
        NULL_USE( normal );
    }
    delete[] tmp;
    return tri_coord;
}


/****************************************************************
 * Create triangles/verticies from a set of triangles specified  *
 * by their coordinates                                          *
 ****************************************************************/
template<size_t NG, size_t NP>
void createTriangles( const std::vector<std::array<std::array<double, NP>, NG + 1>> &tri_list,
                      std::vector<std::array<double, NP>> &verticies,
                      std::vector<std::array<int64_t, NG + 1>> &triangles,
                      double tol )
{
    // Get the range of points and tolerance to use
    std::array<double, 2 * NP> range;
    for ( size_t d = 0; d < NP; d++ ) {
        range[2 * d + 0] = tri_list[0][0][d];
        range[2 * d + 1] = tri_list[0][0][d];
    }
    for ( const auto &tri : tri_list ) {
        for ( const auto &point : tri ) {
            for ( size_t d = 0; d < NP; d++ ) {
                range[2 * d + 0] = std::min( range[2 * d + 0], point[d] );
                range[2 * d + 1] = std::max( range[2 * d + 1], point[d] );
            }
        }
    }
    std::array<double, NP> tol2;
    for ( size_t d = 0; d < NP; d++ )
        tol2[d] = tol * ( range[2 * d + 1] - range[2 * d + 0] );
    // Get the unique verticies and create triangle indicies
    verticies.clear();
    triangles.clear();
    constexpr auto null_tri = make_array<int64_t, NG + 1>( -1 );
    triangles.resize( tri_list.size(), null_tri );
    for ( size_t i = 0; i < tri_list.size(); i++ ) {
        for ( size_t j = 0; j < NG + 1; j++ ) {
            auto &point   = tri_list[i][j];
            int64_t index = -1;
            for ( size_t k = 0; k < verticies.size() && index == -1; k++ ) {
                if ( approx_equal( point, verticies[k], tol2 ) )
                    index = k;
            }
            if ( index == -1 ) {
                index = verticies.size();
                verticies.push_back( point );
            }
            triangles[i][j] = index;
        }
    }
}


/****************************************************************
 * Create triangles neighbors from the triangles                 *
 ****************************************************************/
template<size_t NG>
std::vector<std::array<int64_t, NG + 1>>
create_tri_neighbors( const std::vector<std::array<int64_t, NG + 1>> &tri )
{
    // Allocate memory
    constexpr auto null_tri = make_array<int64_t, NG + 1>( -1 );
    std::vector<std::array<int64_t, NG + 1>> tri_nab( tri.size(), null_tri );
    if ( tri.size() == 1 )
        return tri_nab;
    // 1D is a special easy case
    if constexpr ( NG == 1 ) {
        for ( size_t i = 0; i < tri.size(); i++ ) {
            tri_nab[i][0] = i + 1;
            tri_nab[i][1] = i - 1;
        }
        tri_nab[0][1]              = -1;
        tri_nab[tri.size() - 1][0] = -1;
        return tri_nab;
    }
    PROFILE_START( "create_tri_neighbors", 1 );
    // Get the number of verticies
    size_t N_vertex = 0;
    for ( const auto &t : tri ) {
        for ( size_t i = 0; i < NG + 1; i++ )
            N_vertex = std::max<size_t>( N_vertex, t[i] + 1 );
    }
    // Count the number of triangles connected to each vertex
    std::vector<int64_t> N_tri_nab( N_vertex, 0 );
    for ( size_t i = 0; i < tri.size(); i++ ) {
        for ( size_t d = 0; d < NG + 1; d++ )
            N_tri_nab[tri[i][d]]++;
    }
    // For each node, get a list of the triangles that connect to that node
    auto tri_list = new int64_t *[N_vertex]; // List of triangles connected each node (N)
    tri_list[0]   = new int64_t[( NG + 1 ) * tri.size()];
    for ( size_t i = 1; i < N_vertex; i++ )
        tri_list[i] = &tri_list[i - 1][N_tri_nab[i - 1]];
    for ( size_t i = 0; i < ( NG + 1 ) * tri.size(); i++ )
        tri_list[0][i] = -1;
    // Create a sorted list of all triangles that have each node as a vertex
    for ( size_t i = 0; i < N_vertex; i++ )
        N_tri_nab[i] = 0;
    for ( size_t i = 0; i < tri.size(); i++ ) {
        for ( size_t j = 0; j <= NG; j++ ) {
            int64_t k                 = tri[i][j];
            tri_list[k][N_tri_nab[k]] = i;
            N_tri_nab[k]++;
        }
    }
    for ( size_t i = 0; i < N_vertex; i++ )
        AMP::Utilities::quicksort( N_tri_nab[i], tri_list[i] );
    int64_t N_tri_max = 0;
    for ( size_t i = 0; i < N_vertex; i++ ) {
        if ( N_tri_nab[i] > N_tri_max )
            N_tri_max = N_tri_nab[i];
    }
    // Note, if a triangle is a neighbor, it will share all but the current node
    int size[NG];
    for ( int64_t i = 0; i < (int64_t) tri.size(); i++ ) {
        // Loop through the different faces of the triangle
        for ( size_t j = 0; j <= NG; j++ ) {
            int64_t *list[NG] = { nullptr };
            int64_t k1        = 0;
            for ( size_t k2 = 0; k2 <= NG; k2++ ) {
                if ( k2 == j )
                    continue;
                int64_t k = tri[i][k2];
                list[k1]  = tri_list[k];
                size[k1]  = N_tri_nab[k];
                k1++;
            }
            // Find the intersection of all triangle lists except the current node
            int64_t intersection[5] = { -1, -1, -1, -1, -1 };
            int64_t N_int           = intersect_sorted( NG, size, list, 5, intersection );
            int64_t m               = 0;
            if ( N_int == 0 || N_int > 2 ) {
                // We cannot have less than 1 triangle or more than 2 triangles sharing NDIM nodes
                AMP_ERROR( "Error in create_tri_neighbors detected" );
            } else if ( intersection[0] == i ) {
                m = intersection[1];
            } else if ( intersection[1] == i ) {
                m = intersection[0];
            } else {
                // One of the triangles must be the current triangle
                AMP_ERROR( "Error in create_tri_neighbors detected" );
            }
            tri_nab[i][j] = m;
        }
    }
    // Check tri_nab
    for ( int64_t i = 0; i < (int64_t) tri.size(); i++ ) {
        for ( size_t d = 0; d <= NG; d++ ) {
            if ( tri_nab[i][d] < -1 || tri_nab[i][d] >= (int64_t) tri.size() || tri_nab[i][d] == i )
                AMP_ERROR( "Internal error" );
        }
    }
    delete[] tri_list[0];
    delete[] tri_list;
    PROFILE_STOP( "create_tri_neighbors", 1 );
    return tri_nab;
}


/****************************************************************
 * Create triangles/verticies from a set of triangles specified  *
 * by their coordinates                                          *
 ****************************************************************/
static inline std::array<double, 3> calcNorm( const std::vector<std::array<double, 3>> &x,
                                              const std::array<int64_t, 3> &tri )
{
    return AMP::Geometry::GeometryHelpers::normal( x[tri[0]], x[tri[1]], x[tri[2]] );
}
static inline double dot( const std::array<double, 3> &x, const std::array<double, 3> &y )
{
    return x[0] * y[0] + x[1] * y[1] + x[2] * y[2];
}
template<size_t NG, size_t NP>
static std::vector<int> createBlockIDs( const std::vector<std::array<double, NP>> &verticies,
                                        const std::vector<std::array<int64_t, NG + 1>> &tri,
                                        const std::vector<std::array<int64_t, NG + 1>> &tri_nab )
{
    if ( tri.empty() )
        return std::vector<int>();
    // Calculate the normal for each triangle face
    typedef std::array<double, NP> Point;
    std::vector<Point> norm( tri.size() );
    for ( size_t i = 0; i < tri.size(); i++ )
        norm[i] = calcNorm( verticies, tri[i] );
    // Identify different blocks by the change in the normal
    int nextBlockID = 0;
    std::vector<int> blockID( tri.size(), -1 );
    std::vector<bool> finished( tri.size(), false );
    std::set<size_t> queued;
    double tol = 0.1;
    for ( size_t i = 0; i < tri.size(); i++ ) {
        if ( finished[i] )
            continue;
        blockID[i] = nextBlockID++;
        queued.insert( i );
        while ( !queued.empty() ) {
            auto it  = queued.begin();
            size_t j = *it;
            queued.erase( it );
            finished[j] = true;
            for ( auto k : tri_nab[j] ) {
                if ( k == -1 )
                    continue; // There is no neighbor
                if ( finished[k] )
                    continue; // We already examined this neighbor
                auto theta = acos( dot( norm[j], norm[k] ) );
                if ( theta <= tol ) {
                    // The norm is within tol, set the block id
                    blockID[k] = blockID[j];
                    queued.insert( k );
                }
            }
        }
    }
    return blockID;
}


/****************************************************************
 * Try to split the mesh into seperate independent domains       *
 ****************************************************************/
static inline std::array<int64_t, 2> getFace( const std::array<int64_t, 3> &tri, size_t i )
{
    return { tri[( i + 1 ) % 3], tri[( i + 2 ) % 3] };
}
static inline std::array<int64_t, 3> getFace( const std::array<int64_t, 4> &tri, size_t i )
{
    return { tri[( i + 1 ) % 4], tri[( i + 2 ) % 4], tri[( i + 3 ) % 4] };
}
template<size_t NG>
static void addFaces( const std::array<int64_t, NG + 1> &tri,
                      int64_t index,
                      std::vector<std::pair<uint64_t, int64_t>> &faces )
{
    for ( size_t i = 0; i <= NG; i++ ) {
        // Get each face of the triangle
        auto face = getFace( tri, i );
        auto id1  = hash<NG, true>( face );
        // Reverse the order
        std::reverse( face.begin(), face.end() );
        auto id2   = hash<NG, true>( face );
        bool found = false;
        size_t k   = 0;
        for ( size_t j = 0; j < faces.size(); j++ ) {
            if ( faces[j].first == id1 )
                found = true;
            else
                faces[k++] = faces[j];
        }
        faces.resize( k );
        if ( !found ) {
            // Face does not exist, add it
            int64_t tmp = ( index << 4 ) + i;
            faces.push_back( std::make_pair( id2, tmp ) );
        }
    }
}
template<class TYPE>
static inline void erase( TYPE &faceMap, int64_t i )
{
    for ( auto it = faceMap.begin(); it != faceMap.end(); ) {
        if ( it->second >> 4 == i )
            it = faceMap.erase( it );
        else
            ++it;
    }
}
template<size_t NG>
static std::vector<std::array<int64_t, NG + 1>>
    removeSubDomain( std::vector<std::array<int64_t, NG + 1>> &tri )
{
    // For each triangle get a hash id for each face
    std::multimap<uint64_t, int64_t> faceMap;
    for ( size_t i = 0, k = 0; i < tri.size(); i++ ) {
        for ( size_t j = 0; j <= NG; j++, k++ ) {
            auto face   = getFace( tri[i], j );
            uint64_t id = hash<NG, true>( face );
            int64_t tmp = ( i << 4 ) + j;
            faceMap.insert( std::make_pair( id, tmp ) );
        }
    }
    // Choose an initial triangle
    size_t i0 = 0;
    int count = 100;
    for ( size_t i = 0; i < tri.size(); i++ ) {
        int Nf_max = 0;
        for ( size_t j = 0; j <= NG; j++ ) {
            // Get each face of the triangle
            auto face = getFace( tri[i], j );
            // auto id   = hash<NG, true>( face );
            // Reverse the order
            std::reverse( face.begin(), face.end() );
            auto id2 = hash<NG, true>( face );
            // Get the number of matching faces
            int Nf = faceMap.count( id2 );
            Nf_max = std::max( Nf_max, Nf );
        }
        if ( Nf_max < count ) {
            count = Nf_max;
            i0    = i;
        }
    }
    // Add the initial triangle store the edges
    std::vector<bool> used( tri.size(), false );
    std::vector<std::array<int64_t, NG + 1>> tri2;
    std::vector<std::pair<uint64_t, int64_t>> faces;
    used[i0] = true;
    tri2.push_back( tri[i0] );
    addFaces<NG>( tri[i0], i0, faces );
    erase( faceMap, i0 );
    // Add triangles until all faces have been filled
    while ( !faces.empty() ) {
        bool found = false;
        for ( size_t i = 0; i < faces.size(); i++ ) {
            int Nf = faceMap.count( faces[i].first );
            AMP_ASSERT( Nf > 0 );
            if ( Nf == 1 ) {
                // We are dealing with a unique match, add the triangle
                auto it = faceMap.find( faces[i].first );
                int j   = it->second >> 4;
                used[j] = true;
                tri2.push_back( tri[j] );
                addFaces<NG>( tri[j], j, faces );
                erase( faceMap, j );
                found = true;
                break;
            }
        }
        if ( found )
            continue;
        // We have multiple faces to choose from, try to remove a subdomain from the remaining faces
        try {
            std::vector<std::array<int64_t, NG + 1>> tri3;
            for ( size_t j = 0; j < tri.size(); j++ ) {
                if ( !used[j] )
                    tri3.push_back( tri[j] );
            }
            auto tri4 = removeSubDomain<NG>( tri3 );
            for ( const auto t : tri2 )
                tri3.push_back( t );
            std::swap( tri, tri3 );
            return tri4;
        } catch ( ... ) {
        }
        // Still no luck
        AMP_ERROR( "Unable to resolve multiple faces" );
    }
    // Remove triangles that were used
    size_t k = 0;
    for ( size_t j = 0; j < tri.size(); j++ ) {
        if ( !used[j] )
            tri[k++] = tri[j];
    }
    tri.resize( k );
    return tri2;
}
template<size_t NG>
std::vector<std::vector<std::array<int64_t, NG + 1>>>
    splitDomains( std::vector<std::array<int64_t, NG + 1>> tri )
{
    std::vector<std::vector<std::array<int64_t, NG + 1>>> tri_sets;
    while ( !tri.empty() )
        tri_sets.emplace_back( removeSubDomain<NG>( tri ) );
    return tri_sets;
}
template<>
std::vector<std::vector<std::array<int64_t, 2>>>
    splitDomains<1>( std::vector<std::array<int64_t, 2>> )
{
    AMP_ERROR( "1D splitting of domains is not supported" );
    return std::vector<std::vector<std::array<int64_t, 2>>>();
}


/********************************************************
 *  Generate mesh for STL file                           *
 ********************************************************/
std::shared_ptr<AMP::Mesh::Mesh> generateSTL( std::shared_ptr<MeshParameters> params )
{
    auto db       = params->getDatabase();
    auto filename = db->getWithDefault<std::string>( "FileName", "" );
    auto name     = db->getWithDefault<std::string>( "MeshName", "NULL" );
    auto comm     = params->getComm();
    // Read the STL file
    typedef std::vector<std::array<double, 3>> pointset;
    typedef std::vector<std::array<int64_t, 3>> triset;
    pointset vert;
    std::vector<triset> tri( 1 ), tri_nab;
    if ( comm.getRank() == 0 ) {
        auto scale     = db->getWithDefault<double>( "scale", 1.0 );
        auto triangles = TriangleHelpers::readSTL( filename, scale );
        // Create triangles from the points
        double tol = 1e-6;
        TriangleHelpers::createTriangles<2, 3>( triangles, vert, tri[0], tol );
        // Find the number of unique triangles (duplicates may indicate multiple objects)
        bool multidomain = TriangleHelpers::count<2>( tri[0] ) > 1;
        if ( multidomain && db->getWithDefault<bool>( "split", true ) ) {
            // Try to split the domains
            tri         = TriangleHelpers::splitDomains<2>( tri[0] );
            multidomain = false;
        }
        // Create the triangle neighbors
        tri_nab.resize( tri.size() );
        if ( !multidomain ) {
            for ( size_t i = 0; i < tri.size(); i++ ) {
                // Get the triangle neighbors
                tri_nab[i] = TriangleHelpers::create_tri_neighbors<2>( tri[i] );
                // Check if the geometry is closed
                bool closed = true;
                for ( const auto &t : tri_nab[i] ) {
                    for ( const auto &p : t )
                        closed = closed && p >= 0;
                }
                if ( !closed )
                    AMP_WARNING( "Geometry is not closed" );
            }
        } else {
            AMP_WARNING( "Not splitting multi-domain, no neighbor info will be created" );
            for ( size_t i = 0; i < tri.size(); i++ ) {
                tri_nab[i].resize( tri[i].size() );
                for ( auto &t : tri_nab[i] )
                    t.fill( -1 );
            }
        }
    }
    size_t N_domains = comm.bcast( tri.size(), 0 );
    tri.resize( N_domains );
    tri_nab.resize( N_domains );
    // Create the block ids
    std::vector<std::vector<int>> blocks( N_domains );
    for ( size_t i = 0; i < N_domains; i++ )
        blocks[i] = createBlockIDs<2, 3>( vert, tri[i], tri_nab[i] );
    // Create the mesh
    std::shared_ptr<AMP::Mesh::Mesh> mesh;
    if ( N_domains == 1 ) {
        mesh = TriangleMesh<2, 3>::generate( vert, tri[0], tri_nab[0], comm, nullptr, blocks[0] );
    } else {
        // For now have all meshes on the same communicator
        std::vector<std::shared_ptr<AMP::Mesh::Mesh>> submeshes( N_domains );
        for ( size_t i = 0; i < N_domains; i++ ) {
            submeshes[i] =
                TriangleMesh<2, 3>::generate( vert, tri[i], tri_nab[i], comm, nullptr, blocks[i] );
            submeshes[i]->setName( name + "_" + std::to_string( i + 1 ) );
        }
        mesh.reset( new MultiMesh( name, comm, submeshes ) );
    }
    // Displace the mesh
    std::vector<double> disp( 3, 0.0 );
    if ( db->keyExists( "x_offset" ) )
        disp[0] = db->getScalar<double>( "x_offset" );
    if ( db->keyExists( "y_offset" ) )
        disp[1] = db->getScalar<double>( "y_offset" );
    if ( db->keyExists( "z_offset" ) )
        disp[2] = db->getScalar<double>( "z_offset" );
    if ( disp[0] != 0.0 && disp[1] != 0.0 && disp[2] != 0.0 )
        mesh->displaceMesh( disp );
    // Set the mesh name
    mesh->setName( name );
    return mesh;
}


/********************************************************
 *  Generate mesh for geometry                           *
 ********************************************************/
static inline void check_nearest( const std::vector<Point> &x )
{
    if ( x.empty() )
        return;
    auto index = find_min_dist( x );
    auto dx    = x[index.first] - x[index.second];
    double d   = dx.abs();
    AMP_ASSERT( d > 1e-8 );
}
static inline std::vector<Point> getVolumePoints( const AMP::Geometry::Geometry &geom,
                                                  double resolution )
{
    // Create interior points from an arbitrary geometry
    // Note: we can adjust the points so that they are not aligned
    //    on xyz grids which may help with the tessellation
    std::vector<Point> points;
    auto [x0, x1] = geom.box();
    // auto xc = 0.7 * ( x0 + x1 );
    // double tol = 0.7 * resolution / ( x1 - x0 ).abs();
    for ( double x = x0.x(); x <= x1.x(); x += resolution ) {
        for ( double y = x0.y(); y <= x1.y(); y += resolution ) {
            for ( double z = x0.z(); z <= x1.z(); z += resolution ) {
                Point p( x0.size(), { x, y, z } );
                // auto p2 = p - xc;
                // p -= ( tol * p2.abs() ) * p2;
                if ( geom.inside( p ) )
                    points.push_back( p );
            }
        }
    }
    return points;
}
static inline std::vector<Point> getSurfacePoints( const AMP::Geometry::Geometry &geom, int N )
{
    // Create surface points for an arbitrary geometry
    std::vector<Point> points;
    const int ndim      = geom.getDim();
    const auto [x1, x2] = geom.box();
    const auto x0       = 0.5 * ( x1 + x2 );
    const auto dx       = x2 - x1;
    if ( ndim == 3 ) {
        double r = sqrt( dx.x() * dx.x() + dx.y() * dx.y() + dx.z() * dx.z() );
        int n    = ceil( sqrt( N ) );
        for ( int i = 0; i < n; i++ ) {
            double x = ( 0.5 + i ) / (double) n;
            for ( int j = 0; j < n; j++ ) {
                double y = ( 0.5 + j ) / (double) n;
                Point s  = AMP::Geometry::GeometryHelpers::map_logical_sphere_surface( r, x, y );
                auto dir = normalize( s );
                double d = geom.distance( x0, dir );
                AMP_ASSERT( d < 0 );
                points.push_back( x0 - d * dir );
            }
        }
    } else {
        AMP_ERROR( "Not finished" );
    }
    check_nearest( points );
    return points;
}
static inline double getPos( int i, int N, bool isPeriodic )
{
    if ( N <= 1 )
        return 0;
    if ( isPeriodic )
        return ( i + 0.5 ) / N;
    return i / ( N - 1.0 );
}
static inline std::vector<Point> getLogicalPoints( const AMP::Geometry::LogicalGeometry &geom,
                                                   double resolution )
{
    // Create surface/interior points for a logical geometry
    std::vector<Point> points;
    int ndim = geom.getDim();
    auto N   = geom.getLogicalGridSize( std::vector<double>( ndim, resolution ) );
    N.resize( 3, 1 );
    auto periodic = geom.getPeriodicDim();
    periodic.resize( 3, false );
    for ( int i = 0; i < N[0]; i++ ) {
        double x = getPos( i, N[0], periodic[0] );
        for ( int j = 0; j < N[1]; j++ ) {
            double y = getPos( j, N[1], periodic[1] );
            for ( int k = 0; k < N[2]; k++ ) {
                double z = getPos( k, N[2], periodic[2] );
                Point p  = geom.physical( { x, y, z } );
                points.push_back( p );
            }
        }
    }
    return points;
}
static inline std::vector<Point> combineSurfaceVolumePoints( const std::vector<Point> &volume,
                                                             const std::vector<Point> &surface,
                                                             const AMP::Geometry::Geometry &geom,
                                                             double resolution )
{
    // Add the points in the volume
    std::vector<Point> points = volume;
    // Remove volume points that are close to the surface
    size_t k    = 0;
    double tol  = 0.6 * resolution;
    double tol2 = tol * tol;
    for ( size_t i = 0; i < points.size(); i++ ) {
        auto ps   = geom.nearest( points[i] );
        double d2 = ( points[i] - ps ).norm();
        if ( d2 >= tol2 )
            points[k++] = points[i];
    }
    // Add the surface points
    points.resize( k );
    for ( size_t i = 0; i < surface.size(); i++ )
        points.push_back( surface[i] );
    // Check the distance
    check_nearest( points );
    return points;
}
template<uint8_t NDIM>
static void removeTriangles( std::vector<std::array<int, NDIM + 1>> &tri,
                             std::vector<std::array<int, NDIM + 1>> &tri_nab,
                             const std::vector<bool> &remove )
{
    std::vector<int64_t> map( tri.size(), -1 );
    size_t N = 0;
    for ( size_t i = 0; i < tri.size(); i++ ) {
        if ( !remove[i] )
            map[i] = N++;
    }
    if ( N == tri.size() )
        return;
    for ( size_t i = 0; i < tri.size(); i++ ) {
        if ( !remove[i] ) {
            tri[map[i]] = tri[i];
            for ( int d = 0; d <= NDIM; d++ ) {
                if ( tri_nab[i][d] != -1 ) {
                    if ( remove[tri_nab[i][d]] )
                        tri_nab[i][d] = -1;
                    else
                        tri_nab[i][d] = map[tri_nab[i][d]];
                }
            }
            tri_nab[map[i]] = tri_nab[i];
        }
    }
    tri.resize( N );
    tri_nab.resize( N );
}
template<uint8_t NDIM>
static std::shared_ptr<AMP::Mesh::Mesh> generate( std::shared_ptr<AMP::Geometry::Geometry> geom,
                                                  const std::vector<Point> &points,
                                                  const AMP_MPI &comm )
{
    // Tessellate
    auto [lb, ub] = geom->box();
    auto dx       = ub - lb;
    int N         = points.size();
    std::vector<std::array<double, NDIM>> x1( points.size() );
    std::vector<std::array<int, NDIM>> x2( points.size() );
    for ( int i = 0; i < N; i++ ) {
        for ( int d = 0; d < NDIM; d++ ) {
            double x   = points[i][d];
            double tmp = 2.0 * ( x - lb[d] ) / dx[d] - 1.0;
            int xi     = round( 1000000 * tmp );
            x1[i][d]   = x;
            x2[i][d]   = xi;
        }
    }
    NULL_USE( x2 );
    std::vector<std::array<int, NDIM + 1>> tri, tri_nab;
    try {
        // Try to tessellate with the acutal points
        std::tie( tri, tri_nab ) = DelaunayTessellation::create_tessellation<NDIM>( x1 );
    } catch ( ... ) {
        try {
            // Try to tessellate with the integer points
            std::tie( tri, tri_nab ) = DelaunayTessellation::create_tessellation<NDIM>( x2 );
        } catch ( ... ) {
            // Failed to tessellate
            auto fid = fopen( "failed_points.csv", "wb" );
            for ( const auto &p : x1 ) {
                for ( const auto &v : p )
                    fprintf( fid, "%0.12f ", v );
                fprintf( fid, "\n" );
            }
            fclose( fid );
            AMP_ERROR( "Failed to tessellate" );
        }
    }
    AMP_ASSERT( !tri.empty() );
    // Delete triangles that have duplicate neighbors
    {
        // Identify the triangles that need to be removed
        std::vector<bool> remove( tri.size(), false );
        for ( size_t i = 0; i < tri.size(); i++ ) {
            for ( int d = 0; d <= NDIM; d++ ) {
                if ( tri_nab[i][d] == -1 )
                    continue;
                int count = 0;
                for ( int d2 = 0; d2 <= NDIM; d2++ ) {
                    if ( tri_nab[i][d] == tri_nab[i][d2] )
                        count++;
                }
                if ( count != 1 )
                    remove[i] = true;
            }
        }
        // Remove the triangles
        removeTriangles<NDIM>( tri, tri_nab, remove );
    }
    // Delete surface triangles that have zero volume
    if constexpr ( NDIM == 3 ) {
        // Identify the triangles that need to be removed
        std::vector<bool> remove( tri.size(), false );
        for ( size_t i = 0; i < tri.size(); i++ ) {
            if ( tri_nab[i][0] >= 0 && tri_nab[i][1] >= 0 && tri_nab[i][2] >= 0 &&
                 tri_nab[i][3] >= 0 )
                continue;
            Point x[4];
            for ( int j = 0; j <= NDIM; j++ )
                x[j] = points[tri[i][j]];
            double M[9];
            for ( size_t k = 0; k < 3; k++ ) {
                for ( size_t d = 0; d < 3; d++ )
                    M[d + k * 3] = x[k][d] - x[3][d];
            }
            constexpr double C = 1.0 / 6.0;
            double V           = std::abs( C * DelaunayHelpers<NDIM>::det( M ) );
            remove[i]          = V < 1e-6;
        }
        // Remove the triangles
        removeTriangles<NDIM>( tri, tri_nab, remove );
    }
    // Try to remove triangles outside the domain
    bool isConvex = geom->isConvex();
    if ( !isConvex ) {
        AMP_WARNING( "non-convex domains are not fully supported yet" );
        // Identify the triangles that need to be removed
        std::vector<bool> remove( tri.size(), false );
        const double tmp = 1.0 / ( NDIM + 1.0 );
        for ( size_t i = 0; i < tri.size(); i++ ) {
            Point center( NDIM, { 0, 0, 0 } );
            for ( int j = 0; j <= NDIM; j++ )
                center += points[tri[i][j]];
            center *= tmp;
            remove[i] = !geom->inside( center );
        }
        // Remove the triangles
        removeTriangles<NDIM>( tri, tri_nab, remove );
    }
    // Generate the mesh
    std::vector<std::array<int64_t, NDIM + 1>> tri2( tri.size() ), tri_nab2( tri.size() );
    for ( size_t i = 0; i < tri.size(); i++ ) {
        for ( int d = 0; d <= NDIM; d++ ) {
            tri2[i][d]     = tri[i][d];
            tri_nab2[i][d] = tri_nab[i][d];
        }
    }
    return TriangleMesh<NDIM, NDIM>::generate(
        std::move( x1 ), std::move( tri2 ), std::move( tri_nab2 ), comm, std::move( geom ) );
}
std::shared_ptr<AMP::Mesh::Mesh>
generate( std::shared_ptr<AMP::Geometry::Geometry> geom, const AMP_MPI &comm, double resolution )
{
    auto multigeom = std::dynamic_pointer_cast<AMP::Geometry::MultiGeometry>( geom );
    if ( multigeom ) {
        std::vector<std::shared_ptr<AMP::Mesh::Mesh>> submeshes;
        for ( auto &geom2 : multigeom->getGeometries() )
            submeshes.push_back( generate( geom2, comm, resolution ) );
        return std::make_shared<MultiMesh>( "name", comm, submeshes );
    }
    // Perform some basic checks
    int ndim         = geom->getDim();
    auto meshGeom    = std::dynamic_pointer_cast<AMP::Geometry::MeshGeometry>( geom );
    auto logicalGeom = std::dynamic_pointer_cast<AMP::Geometry::LogicalGeometry>( geom );
    // Create the grid verticies
    std::vector<Point> points;
    if ( logicalGeom ) {
        // We are dealing with a logical geometry
        points = getLogicalPoints( *logicalGeom, resolution );
    } else if ( meshGeom ) {
        // Get the volume points
        auto interior = getVolumePoints( *geom, resolution );
        // Get the surface points
        auto &mesh   = meshGeom->getMesh();
        auto data    = sample( mesh, resolution );
        auto surface = std::get<0>( data );
        // Combine
        points = combineSurfaceVolumePoints( interior, surface, *geom, resolution );
    } else {
        // Get the volume points
        auto interior = getVolumePoints( *geom, resolution );
        // Get the surface points
        auto surface = getSurfacePoints( *geom, 0.1 * interior.size() );
        // Combine
        points = combineSurfaceVolumePoints( interior, surface, *geom, resolution );
    }
    // Smooth the points to try and make the distance between all points ~ equal

    // Tessellate and generate the mesh
    std::shared_ptr<AMP::Mesh::Mesh> mesh;
    if ( ndim == 2 ) {
        mesh = generate<2>( geom, points, comm );
    } else if ( ndim == 3 ) {
        mesh = generate<3>( geom, points, comm );
    } else {
        AMP_ERROR( "Not supported yet" );
    }
    if ( meshGeom )
        mesh->setName( meshGeom->getMesh().getName() );
    return mesh;
}


/********************************************************
 *  Explicit instantiations                              *
 ********************************************************/
// clang-format off
typedef std::array<double,1> point1D;
typedef std::array<double,2> point2D;
typedef std::array<double,3> point3D;
typedef std::vector<point1D> pointset1D;
typedef std::vector<point2D> pointset2D;
typedef std::vector<point3D> pointset3D;
typedef std::vector<std::array<int64_t,2>> triset1D;
typedef std::vector<std::array<int64_t,3>> triset2D;
typedef std::vector<std::array<int64_t,4>> triset3D;
template size_t count<1>( const triset1D & );
template size_t count<2>( const triset2D & );
template size_t count<3>( const triset3D & );
template void createTriangles<1,1>( const std::vector<std::array<point1D,2>>&, pointset1D&, triset1D&, double );
template void createTriangles<1,2>( const std::vector<std::array<point2D,2>>&, pointset2D&, triset1D&, double );
template void createTriangles<1,3>( const std::vector<std::array<point3D,2>>&, pointset3D&, triset1D&, double );
template void createTriangles<2,2>( const std::vector<std::array<point2D,3>>&, pointset2D&, triset2D&, double );
template void createTriangles<2,3>( const std::vector<std::array<point3D,3>>&, pointset3D&, triset2D&, double );
template void createTriangles<3,3>( const std::vector<std::array<point3D,4>>&, pointset3D&, triset3D&, double );
template triset1D create_tri_neighbors<1>( const triset1D& );
template triset2D create_tri_neighbors<2>( const triset2D& );
template triset3D create_tri_neighbors<3>( const triset3D& );
template std::vector<triset2D> splitDomains<2>( triset2D );
template std::vector<triset3D> splitDomains<3>( triset3D );
// clang-format on

} // namespace TriangleHelpers
} // namespace Mesh
} // namespace AMP
