#include "AMP/ampmesh/triangle/TriangleMesh.h"
#include "AMP/ampmesh/triangle/TriangleMeshIterator.h"

#include "AMP/ampmesh/MultiIterator.h"


#include "AMP/utils/Utilities.h"
#ifdef USE_AMP_VECTORS
#include "AMP/vectors/Variable.h"
#include "AMP/vectors/Vector.h"
#include "AMP/vectors/VectorBuilder.h"
#endif

#include "ProfilerApp.h"

#include <algorithm>
#include <cstring>
#include <iostream>


namespace AMP {
namespace Mesh {


// Helper function to create constexpr std::array with a single value
template <class T, std::size_t N, std::size_t... I>
static constexpr std::array<std::remove_cv_t<T>, N>
    to_array_impl( const T* a, std::index_sequence<I...>)
{
    return { {a[I]...} };
}
template<class TYPE, std::size_t N>
static constexpr std::array<TYPE,N> make_array( const TYPE& x )
{
    TYPE tmp[10] = { x, x, x, x, x, x, x, x, x, x };
    return to_array_impl<TYPE,N>(tmp, std::make_index_sequence<N>{});
}


// Helper function to copy arrays if they are the same type
template<class T1, class T2>
static inline typename std::enable_if<std::is_same<T1,T2>::value>::type 
copyIfSame( const std::vector<T1>& src, std::vector<T2>& dst ) { dst = src; }
template<class T1, class T2>
static inline typename std::enable_if<!std::is_same<T1,T2>::value>::type 
copyIfSame( const std::vector<T1>&, std::vector<T2>& ) { }


// Helper functions to see if two points are ~ the same
static inline bool approx_equal( const std::array<double,1>& x, const std::array<double,1>& y, const std::array<double,1>& tol )
{
    return fabs( x[0] - y[0] ) <= tol[0];
}
static inline bool approx_equal( const std::array<double,2>& x, const std::array<double,2>& y, const std::array<double,2>& tol )
{
    return fabs( x[0] - y[0] ) <= tol[0] && fabs( x[1] - y[1] ) <= tol[1];
}
static inline bool approx_equal( const std::array<double,3>& x, const std::array<double,3>& y, const std::array<double,3>& tol )
{
    return fabs( x[0] - y[0] ) <= tol[0] && fabs( x[1] - y[1] ) <= tol[1] && fabs( x[2] - y[2] ) <= tol[2];
}


// Subroutine to find the first n intersections in multiple lists
// This function assumes the lists are in sorted order
static int intersect_sorted( const int N_lists, const int size[], uint64_t *list[],
    const int N_max, uint64_t *intersection )
{
    if ( N_max <= 0 )
        return ~( (unsigned int) 0 );
    int N_int  = 0;
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
 * Read stl file                                                 *
 ****************************************************************/
static std::vector<std::array<std::array<double,3>,3>> readSTL( const std::string& filename, double scale )
{
    char header[80];
    uint32_t N;
    // Read the file
    auto fid = fopen( filename.c_str(), "rb" );
    AMP_INSIST( fid, "Unable to open " + filename );
    fread( header, sizeof( header ), 1, fid );
    fread( &N, sizeof( N ), 1, fid );
    auto tmp = new char[N*50];
    fread( tmp, 50, N, fid );
    fclose( fid );
    // Get a list of the local triangles based on their coordinates
    std::vector<std::array<std::array<double,3>,3>> triangles( N );
    for ( size_t i=0; i<N; i++) {
        float tmp2[9];
        memcpy( tmp2, &tmp[50*i+12], sizeof(tmp2) );
        triangles[i][0][0] = scale * tmp2[0];
        triangles[i][0][1] = scale * tmp2[1];
        triangles[i][0][2] = scale * tmp2[2];
        triangles[i][1][0] = scale * tmp2[3];
        triangles[i][1][1] = scale * tmp2[4];
        triangles[i][1][2] = scale * tmp2[5];
        triangles[i][2][0] = scale * tmp2[6];
        triangles[i][2][1] = scale * tmp2[7];
        triangles[i][2][2] = scale * tmp2[8];
    }
    delete [] tmp;
    return triangles;
}


/****************************************************************
 * Create triangles/verticies from a set of triangles specified  *
 * by their coordinates                                          *
 ****************************************************************/
template<size_t NG, size_t NP>
static void createTriangles( const std::vector<std::array<std::array<double,NP>,NG+1>>& tri_list,
    std::vector<std::array<double,NP>>& verticies,
    std::vector<std::array<uint64_t,NG+1>>& triangles,
    double tol )
{
    // Get the range of points and tolerance to use
    std::array<double,2*NP> range;
    for ( size_t d=0; d<NP; d++) {
        range[2*d+0] = tri_list[0][0][d];
        range[2*d+1] = tri_list[0][0][d];
    }
    for ( const auto& tri : tri_list ) {
        for ( const auto& point : tri ) {
            for ( size_t d=0; d<NP; d++) {
                range[2*d+0] = std::min( range[2*d+0], point[d] );
                range[2*d+1] = std::max( range[2*d+1], point[d] );
            }
        }
    }
    std::array<double,NP> tol2;
    for ( size_t d=0; d<NP; d++)
        tol2[d] = tol * ( range[2*d+1] - range[2*d+0] );
    // Get the unique verticies and create triangle indicies
    verticies.clear();
    triangles.clear();
    constexpr auto null_tri = make_array<uint64_t,NG+1>( static_cast<int64_t>( -1 ) );
    triangles.resize( tri_list.size(), null_tri );
    for ( size_t i=0; i<tri_list.size(); i++) {
        for ( size_t j=0; j<NG+1; j++) {
            auto& point = tri_list[i][j];
            int64_t index = -1;
            for ( size_t j=0; j<verticies.size() && index == -1; j++) {
                if ( approx_equal( point, verticies[j], tol2 ) )
                    index = j;
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
static std::vector<std::array<uint64_t,NG+1>> create_tri_neighbors( const std::vector<std::array<uint64_t,NG+1>>& tri )
{
    // Allocate memory
    uint64_t neg_one = static_cast<int64_t>( -1 );
    std::vector<std::array<uint64_t,NG+1>> tri_nab( tri.size() );
    for ( auto &t : tri_nab )
        t.fill( neg_one );
    if ( tri.size() == 1 )
        return tri_nab;
    // 1D is a special easy case
    if ( NG == 1 ) {
        for ( size_t i = 0; i < tri.size(); i++ ) {
            tri_nab[i][0] = i + 1;
            tri_nab[i][1] = i - 1;
        }
        tri_nab[0][1]              = neg_one;
        tri_nab[tri.size() - 1][0] = neg_one;
        return tri_nab;
    }
    PROFILE_START( "create_tri_neighbors", 1 );
    // Get the number of verticies
    size_t N_vertex = 0;
    for ( const auto &t : tri ) {
        for ( size_t i=0; i<NG+1; i++)
            N_vertex = std::max<size_t>( N_vertex, t[i]+1 );
    }
    // Count the number of triangles connected to each vertex
    std::vector<uint64_t> N_tri_nab( N_vertex, 0 );
    for ( size_t i = 0; i < tri.size(); i++ ) {
        for ( size_t d = 0; d < NG + 1; d++ )
            N_tri_nab[tri[i][d]]++;
    }
    // For each node, get a list of the triangles that connect to that node
    auto tri_list = new uint64_t *[N_vertex]; // List of triangles connected each node (N)
    tri_list[0]   = new uint64_t[( NG + 1 ) * tri.size()];
    for ( size_t i = 1; i < N_vertex; i++ )
        tri_list[i] = &tri_list[i - 1][N_tri_nab[i - 1]];
    for ( size_t i = 0; i < ( NG + 1 ) * tri.size(); i++ )
        tri_list[0][i] = neg_one;
    // Create a sorted list of all triangles that have each node as a vertex
    for ( size_t i = 0; i < N_vertex; i++ )
        N_tri_nab[i] = 0;
    for ( size_t i = 0; i < tri.size(); i++ ) {
        for ( size_t j = 0; j <= NG; j++ ) {
            uint64_t k                 = tri[i][j];
            tri_list[k][N_tri_nab[k]] = i;
            N_tri_nab[k]++;
        }
    }
    for ( size_t i = 0; i < N_vertex; i++ )
        AMP::Utilities::quicksort( N_tri_nab[i], tri_list[i] );
    uint64_t N_tri_max = 0;
    for ( size_t i = 0; i < N_vertex; i++ ) {
        if ( N_tri_nab[i] > N_tri_max )
            N_tri_max = N_tri_nab[i];
    }
    // Note, if a triangle is a neighbor, it will share all but the current node
    int size[NG];
    int error = 0;
    for ( uint64_t i = 0; i < tri.size(); i++ ) {
        // Loop through the different faces of the triangle
        for ( size_t j = 0; j <= NG; j++ ) {
            uint64_t *list[NG] = { nullptr };
            uint64_t k1        = 0;
            for ( size_t k2 = 0; k2 <= NG; k2++ ) {
                if ( k2 == j )
                    continue;
                uint64_t k = tri[i][k2];
                list[k1] = tri_list[k];
                size[k1] = N_tri_nab[k];
                k1++;
            }
            // Find the intersection of all triangle lists except the current node
            uint64_t intersection[5] = { neg_one, neg_one, neg_one, neg_one, neg_one };
            uint64_t N_int           = intersect_sorted( NG, size, list, 5, intersection );
            uint64_t m               = 0;
            if ( N_int == 0 || N_int > 2 ) {
                // We cannot have less than 1 triangle or more than 2 triangles sharing NDIM nodes
                error = 1;
                break;
            } else if ( intersection[0] == i ) {
                m = intersection[1];
            } else if ( intersection[1] == i ) {
                m = intersection[0];
            } else {
                // One of the triangles must be the current triangle
                error = 1;
                break;
            }
            tri_nab[i][j] = m;
        }
        if ( error != 0 )
            break;
    }
    // Check tri_nab
    for ( uint64_t i = 0; i < tri.size(); i++ ) {
        for ( size_t d = 0; d <= NG; d++ ) {
            if ( tri_nab[i][d] < neg_one || tri_nab[i][d] >= tri.size() ||
                 tri_nab[i][d] == i )
                error = 2;
        }
    }
    delete[] tri_list[0];
    delete[] tri_list;
    if ( error == 1 ) {
        throw std::logic_error( "Error in create_tri_neighbors detected" );
    } else if ( error == 2 ) {
        throw std::logic_error( "Internal error" );
    }
    PROFILE_STOP( "create_tri_neighbors", 1 );
    return tri_nab;
}


/****************************************************************
 * Generator                                                     *
 ****************************************************************/
template<size_t NG, size_t NP>
static AMP::shared_ptr<TriangleMesh<NG,NP>> generateSTL( const std::string& filename, const AMP_MPI& comm, double scale );
template<>
AMP::shared_ptr<TriangleMesh<2,3>> generateSTL( const std::string& filename, const AMP_MPI& comm, double scale )
{
    // Read the STL file
    std::vector<std::array<std::array<double,3>,3>> triangles;
    if ( comm.getRank() == 0 )
        triangles = readSTL( filename, scale );
    // Create the mesh from the local triangles
    auto mesh = TriangleMesh<2,3>::generate( triangles, comm, 1e-6 );
    return mesh;
}
template<size_t NG, size_t NP>
static AMP::shared_ptr<TriangleMesh<NG,NP>> generateSTL( const std::string&, const AMP_MPI&, double )
{
    AMP_ERROR( "STL meshes are only supported for NG=2, NP=3" );
    return AMP::shared_ptr<TriangleMesh<NG,NP>>();
}
template<size_t NG, size_t NP>
AMP::shared_ptr<TriangleMesh<NG,NP>> TriangleMesh<NG,NP>::generate( MeshParameters::shared_ptr params )
{
    auto db = params->getDatabase();
    // Create the mesh
    auto filename = db->getStringWithDefault( "FileName", "" );
    AMP::shared_ptr<TriangleMesh<NG,NP>> mesh;
    if ( filename.substr(std::max<int>(filename.length(),4)-4) == ".stl" ) {
        auto scale = db->getDoubleWithDefault( "scale", 1.0 );
        mesh = generateSTL<NG,NP>( filename, params->getComm(), scale );
        mesh->setName( db->getStringWithDefault( "MeshName", "NULL" ) );
    } else {
        AMP_ERROR("Not finished");
    }
    // Displace the mesh
    std::vector<double> displacement( NP, 0.0 );
    if ( db->keyExists( "x_offset" ) && NP >= 1 )
        displacement[0] = db->getDouble( "x_offset" );
    if ( db->keyExists( "y_offset" ) && NP >= 2 )
        displacement[1] = db->getDouble( "y_offset" );
    if ( db->keyExists( "z_offset" ) && NP >= 3 )
        displacement[2] = db->getDouble( "z_offset" );
    bool test = false;
    for ( auto &elem : displacement ) {
        if ( elem != 0.0 )
            test = true;
    }
    if ( test )
        mesh->displaceMesh( displacement );
    return mesh;
}
template<size_t NG, size_t NP>
AMP::shared_ptr<TriangleMesh<NG,NP>> TriangleMesh<NG,NP>::generate(
    const std::vector<std::array<std::array<double,NP>,NG+1>>& tri_list, const AMP_MPI& comm, double tol )
{
    // Get the global list of tri_list
    auto global_list = comm.allGather( tri_list );
    // Create triangles from the points
    std::vector<std::array<double,NP>> verticies;
    std::vector<std::array<uint64_t,NG+1>> triangles;
    createTriangles<NG,NP>( global_list, verticies, triangles, tol );
    // Get the triangle neighbors
    //auto neighbors = create_tri_neighbors<NG>( triangles );
    std::vector<std::array<uint64_t,NG+1>> neighbors( triangles.size(), make_array<uint64_t,NG+1>( -1 ) );
    // Create the mesh
    AMP::shared_ptr<TriangleMesh<NG,NP>> mesh( new TriangleMesh<NG,NP>( verticies, triangles, neighbors, comm ) );
    return mesh;
}
template<size_t N>
static inline void relabel( std::vector<std::array<uint64_t,N>>& data, GeomType type )
{
    for ( auto& shape : data ) {
        for ( auto& i : shape )
            i = MeshElementID( true, type, i, 0, 0 ).elemID();
    }
}
template<size_t NG, size_t NP>
TriangleMesh<NG,NP>::TriangleMesh( const std::vector<std::array<double,NP>>& verticies,
    const std::vector<std::array<uint64_t,NG+1>>& tri, const std::vector<std::array<uint64_t,NG+1>>& tri_nab, const AMP_MPI& comm )
{
    AMP_ASSERT( tri.size() == 0 || comm.getRank() == 0 );
    AMP_ASSERT( tri_nab.size() == tri.size() );
    // Set basic mesh info
    d_db        = nullptr;
    d_params    = nullptr;
    d_geometry  = nullptr;
    GeomDim     = static_cast<GeomType>( NG );
    PhysicalDim = NP;
    d_max_gcw   = 0;
    d_comm      = comm;
    d_name      = "NULL";
    setMeshID();
    // Store the results so far
    d_vert = verticies;
    d_neighbors = tri_nab;
    copyIfSame( tri, d_edge );
    copyIfSame( tri, d_tri );
    copyIfSame( tri, d_tet );
    d_N_global[0] = d_vert.size();
    d_N_global[1] = d_edge.size();
    d_N_global[2] = d_tri.size();
    d_N_global[3] = d_tet.size();
    // Relabel the ids with the type and rank
    relabel( d_edge, GeomType::Vertex );
    relabel( d_tri, GeomType::Vertex );
    relabel( d_tet, GeomType::Vertex );
    relabel( d_neighbors, static_cast<GeomType>( NG ) );
    // Perform load balancing
    if ( comm.getSize() == 1 ) {
        d_max_gcw = 10;
    } else {
        comm.bcast( d_N_global.data(), d_N_global.size(), 0 );
        // Get the triangle centers (used for load balance)
        std::vector<std::array<double,NP>> center( tri.size(), make_array<double,NP>( 0 ) );
        for ( size_t i=0; i<tri.size(); i++) {
            for (size_t j=0; j<NG+1; j++) {
                const auto& point = verticies[tri[i][j]];
                for (size_t d=0; d<NP; d++)
                    center[i][d] += point[d] / NP;
            }
        }
        // Load balance
        d_max_gcw = 3;
        AMP_ERROR("Not finished");
    }
    // Initialize the iterators and some common data
    initialize();
    // Create the surface, block, and boundary iterators
    d_block_iterators.resize( 1 );
    d_block_iterators[0] = d_iterators;
}


/****************************************************************
 * Initialize mesh data                                          *
 ****************************************************************/
static AMP::shared_ptr<std::vector<uint64_t>> createLocalList( size_t N, GeomType type, int rank )
{
    auto list = AMP::make_shared<std::vector<uint64_t>>( N );
    for ( size_t i=0; i<N; i++)
        (*list)[i] = MeshElementID( true, type, i, rank, 0 ).elemID();
    return list;
}
template<size_t NG, size_t NP>
void TriangleMesh<NG,NP>::initialize( )
{
    // Get the bounding boxes
    d_box.resize( 2 * NP );
    d_box_local.resize( 2 * NP );
    for ( const auto& p : d_vert ) {
        for (size_t d=0; d<NP; d++) {
            d_box_local[2*d+0] = std::min( d_box_local[2*d+0], p[d] );
            d_box_local[2*d+1] = std::max( d_box_local[2*d+1], p[d] );
        }
    }
    for (size_t d=0; d<NP; d++) {
        d_box[2*d+0] = d_comm.minReduce( d_box_local[2*d+0] );
        d_box[2*d+1] = d_comm.maxReduce( d_box_local[2*d+1] );
    }
    // Initialize the iterators
    int max_gcw = d_comm.getSize() == 1 ? 0:d_max_gcw;
    d_iterators.resize( ( max_gcw + 1 ) * ( NG + 1 ) );
    d_iterators[0] = TriangleMeshIterator<NG,NP>( this, createLocalList( d_vert.size(), GeomType::Vertex, d_comm.getRank() ) );
    if ( NG >= 1 )
        d_iterators[1] = TriangleMeshIterator<NG,NP>( this, createLocalList( d_edge.size(), GeomType::Edge, d_comm.getRank() ) );
    if ( NG >= 2 )
        d_iterators[2] = TriangleMeshIterator<NG,NP>( this, createLocalList( d_tri.size(), GeomType::Face, d_comm.getRank() ) );
    if ( NG >= 3 )
        d_iterators[3] = TriangleMeshIterator<NG,NP>( this, createLocalList( d_tet.size(), GeomType::Volume, d_comm.getRank() ) );
    for ( int gcw = 1; gcw <= max_gcw; gcw++ ) {
        AMP_ERROR( "Not finished" );
    }
}


/****************************************************************
 * Estimate the mesh size                                        *
 ****************************************************************/
template<size_t NG, size_t NP>
size_t TriangleMesh<NG,NP>::estimateMeshSize( const MeshParameters::shared_ptr &params )
{
    size_t N = 0;
    auto db  = params->getDatabase();
    auto filename = db->getStringWithDefault( "FileName", "" );
    if ( filename.substr(std::max<int>(filename.length(),4)-4) == ".stl" ) {
        // We are reading an stl file
        char header[80];
        uint32_t N2;
        auto fid = fopen( filename.c_str(), "rb" );
        AMP_INSIST( fid, "Unable to open " + filename );
        fread( header, sizeof( header ), 1, fid );
        fread( &N2, sizeof( N2 ), 1, fid );
        fclose( fid );
        N = N2;
    } else {
        AMP_ERROR("Not finished");
    }
    return N;
}


/****************************************************************
 * Constructor                                                   *
 ****************************************************************/
template<size_t NG, size_t NP>
TriangleMesh<NG,NP>::TriangleMesh( MeshParameters::shared_ptr params_in ) : Mesh( params_in )
{
    // Check for valid inputs
    AMP_INSIST( d_params != nullptr, "Params must not be null" );
    AMP_INSIST( !d_comm.isNull(), "Communicator must be set" );
    AMP_INSIST( d_db.get(), "Database must exist" );
    AMP_ERROR("Not finished");
}
template<size_t NG, size_t NP>
TriangleMesh<NG,NP>::TriangleMesh( const TriangleMesh &rhs ) : Mesh( rhs.d_params )
{
    d_N_global         = rhs.d_N_global;
    d_vert             = rhs.d_vert;
    d_edge             = rhs.d_edge;
    d_tri              = rhs.d_tri;
    d_tet              = rhs.d_tet;
    d_neighbors        = rhs.d_neighbors;
    d_remote_vert      = rhs.d_remote_vert;
    d_remote_edge      = rhs.d_remote_edge;
    d_remote_tri       = rhs.d_remote_tri;
    d_remote_tet       = rhs.d_remote_tet;
    d_remote_neighbors = rhs.d_remote_neighbors;
}
template<size_t NG, size_t NP>
AMP::shared_ptr<Mesh> TriangleMesh<NG,NP>::clone() const
{
    return AMP::shared_ptr<TriangleMesh<NG,NP>>( new TriangleMesh<NG,NP>( *this ) );
}


/****************************************************************
 * De-constructor                                                *
 ****************************************************************/
template<size_t NG, size_t NP>
TriangleMesh<NG,NP>::~TriangleMesh() = default;


/****************************************************************
 * Estimate the maximum number of processors                     *
 ****************************************************************/
template<size_t NG, size_t NP>
size_t TriangleMesh<NG,NP>::maxProcs( const MeshParameters::shared_ptr &params )
{
    return estimateMeshSize( params );
}


/****************************************************************
 * Function to return the element given an ID                    *
 ****************************************************************/
template<size_t NG, size_t NP>
MeshElement TriangleMesh<NG,NP>::getElement( const MeshElementID &id ) const
{
    return TriangleMeshElement<NG,NP>( id, this );
}


/********************************************************
 * Function to return parents of an element              *
 ********************************************************/
template<size_t NG, size_t NP>
std::vector<MeshElement> TriangleMesh<NG,NP>::getElementParents( const MeshElement &meshelem,
                                                     const GeomType type ) const
{
    NULL_USE( meshelem );
    NULL_USE( type );
    AMP_ERROR("Not finished");
    return std::vector<MeshElement>();
}


/****************************************************************
 * Functions to return the number of elements                    *
 ****************************************************************/
template<size_t NG, size_t NP>
size_t TriangleMesh<NG,NP>::numLocalElements( const GeomType type ) const
{
    if ( type == GeomType::Vertex )
        return d_vert.size();
    if ( type == GeomType::Edge )
        return d_edge.size();
    if ( type == GeomType::Face )
        return d_tri.size();
    if ( type == GeomType::Volume )
        return d_tet.size();
    return 0;
}
template<size_t NG, size_t NP>
size_t TriangleMesh<NG,NP>::numGlobalElements( const GeomType type ) const
{
    if ( static_cast<uint8_t>( type ) <= NG )
        return d_N_global[ static_cast<uint8_t>( type ) ];
    return 0;
}
template<size_t NG, size_t NP>
size_t TriangleMesh<NG,NP>::numGhostElements( const GeomType type, int gcw ) const
{
    if ( gcw == 0 || d_comm.getSize() == 1 )
        return 0;
    int index1 = static_cast<int>( type );
    int index2 = static_cast<int>( type ) + gcw * ( NG + 1 );
    return d_iterators[index2].size() - d_iterators[index1].size();
}


/****************************************************************
 * Function to get an iterator                                   *
 ****************************************************************/
template<size_t NG, size_t NP>
MeshIterator TriangleMesh<NG,NP>::getIterator( const GeomType type, const int gcw ) const
{
    if ( static_cast<size_t>( type ) > NG || gcw > d_max_gcw )
        return MeshIterator();
    int gcw2 = d_comm.getSize() == 1 ? 0 : gcw;
    int index = static_cast<int>( type ) + gcw2 * ( NG + 1 );
    return d_iterators[index];
}


/****************************************************************
 * Function to get an iterator over the surface                  *
 ****************************************************************/
template<size_t NG, size_t NP>
MeshIterator TriangleMesh<NG,NP>::getSurfaceIterator( const GeomType type, const int gcw ) const
{
    if ( static_cast<size_t>( type ) > NG || gcw > d_max_gcw )
        return MeshIterator();
    int gcw2 = d_comm.getSize() == 1 ? 0 : gcw;
    int index = static_cast<int>( type ) + gcw2 * ( NG + 1 );
    if ( index > (int) d_surface_iterators.size() )
        return MeshIterator();
    return d_surface_iterators[index];
}


/****************************************************************
 * Functions to get the boundaries                               *
 ****************************************************************/
template<size_t NG, size_t NP>
std::vector<int> TriangleMesh<NG,NP>::getBoundaryIDs() const
{
    std::vector<int> ids( d_boundary_iterators.size() );
    for (size_t i=0; i<d_boundary_iterators.size(); i++)
        ids[i] = i;
    return ids;
}
template<size_t NG, size_t NP>
MeshIterator
TriangleMesh<NG,NP>::getBoundaryIDIterator( const GeomType type, const int id, const int gcw ) const
{
    if ( static_cast<size_t>( type ) > NG || gcw > d_max_gcw || id > (int) d_boundary_iterators.size() )
        return MeshIterator();
    int gcw2 = d_comm.getSize() == 1 ? 0 : gcw;
    int index = static_cast<int>( type ) + gcw2 * ( NG + 1 );
    if ( index > (int) d_boundary_iterators[id].size() )
        return MeshIterator();
    return d_boundary_iterators[id][index];
}
template<size_t NG, size_t NP>
std::vector<int> TriangleMesh<NG,NP>::getBlockIDs() const
{
    std::vector<int> ids( d_block_iterators.size() );
    for (size_t i=0; i<d_block_iterators.size(); i++)
        ids[i] = i;
    return ids;
}
template<size_t NG, size_t NP>
MeshIterator TriangleMesh<NG,NP>::getBlockIDIterator( const GeomType type, const int id, const int gcw ) const
{
    if ( static_cast<size_t>( type ) > NG || gcw > d_max_gcw || id > (int) d_block_iterators.size() )
        return MeshIterator();
    int gcw2 = d_comm.getSize() == 1 ? 0 : gcw;
    int index = static_cast<int>( type ) + gcw2 * ( NG + 1 );
    if ( index > (int) d_block_iterators[id].size() )
        return MeshIterator();
    return d_block_iterators[id][index];
}


/****************************************************************
 * Functions to dispace the mesh                                 *
 ****************************************************************/
template<size_t NG, size_t NP>
void TriangleMesh<NG,NP>::displaceMesh( const std::vector<double> &x )
{
    AMP_ASSERT( x.size() == NP );
    for ( auto& p : d_vert ) {
        for ( size_t d=0; d<NP; d++)
            p[d] += x[d];
    }
    for ( auto& p : d_remote_vert ) {
        for ( size_t d=0; d<NP; d++)
            p.second[d] += x[d];
    }
}
#ifdef USE_AMP_VECTORS
template<size_t NG, size_t NP>
void TriangleMesh<NG,NP>::displaceMesh( AMP::shared_ptr<const AMP::LinearAlgebra::Vector> x )
{
    NULL_USE( x );
    AMP_ERROR("Not finished");
}
#endif


/****************************************************************
 *  Get the coordinated of the given vertex or the centroid      *
 ****************************************************************/
template<size_t NG, size_t NP>
const std::array<double,NP>& TriangleMesh<NG,NP>::getPos( const MeshElementID& id ) const
{
    if ( id.is_local() )
        return d_vert[id.local_id()];
    auto it = d_remote_vert.find( id );
    return it->second;
}


/****************************************************************
 * Return the IDs of the elements composing the current element  *
 ****************************************************************/
template<size_t NG, size_t NP>
inline void TriangleMesh<NG,NP>::getVerticies( const MeshElementID& id, int& N, uint64_t* IDs ) const
{
    auto type = id.type();
    if ( type == GeomType::Vertex ) {
        N = 1;
        IDs[0] = id.elemID();
        return;
    }
    if ( type == GeomType::Edge ) {
        N = 2;
        Edge edge;
        if ( id.is_local() )
            edge = d_edge[id.local_id()];
        else
            edge = d_remote_edge.find( id )->second;
        IDs[0] = edge[0];
        IDs[1] = edge[1];
    } else if ( type == GeomType::Face ) {
        N = 3;
        Triangle tri;
        if ( id.is_local() )
            tri = d_tri[id.local_id()];
        else
            tri = d_remote_tri.find( id )->second;
        IDs[0] = tri[0];
        IDs[1] = tri[1];
        IDs[2] = tri[2];
    } else if ( type == GeomType::Volume ) {
        N = 4;
        Tetrahedron tet;
        if ( id.is_local() )
            tet = d_tet[id.local_id()];
        else
            tet = d_remote_tet.find( id )->second;
        IDs[0] = tet[0];
        IDs[1] = tet[1];
        IDs[2] = tet[2];
        IDs[3] = tet[3];
    } else {
        AMP_ERROR( "Not finished" );
    }
}
template<size_t NG, size_t NP>
void TriangleMesh<NG,NP>::getElementsIDs( const MeshElementID& id, const GeomType type, MeshElementID* IDs ) const
{
    int N_vertex = 0;
    uint64_t vertex_ids[NG+1];
    getVerticies( id, N_vertex, vertex_ids );
    if ( type == GeomType::Vertex ) {
        for ( int i=0; i<N_vertex; i++)
            IDs[i] = MeshElementID( d_meshID, vertex_ids[i] );
        return;
    }

    AMP_ERROR( "Not finished" );
}


/********************************************************
 *  Explicit instantiations of TriangleMesh              *
 ********************************************************/
template class TriangleMesh<1,1>;
template class TriangleMesh<1,2>;
template class TriangleMesh<1,3>;
template class TriangleMesh<2,2>;
template class TriangleMesh<2,3>;
template class TriangleMesh<3,3>;


} // namespace Mesh
} // namespace AMP
