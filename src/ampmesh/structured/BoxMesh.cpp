#include "AMP/ampmesh/structured/BoxMesh.h"
#include "AMP/ampmesh/MeshParameters.h"
#include "AMP/ampmesh/MultiIterator.h"
#include "AMP/ampmesh/structured/MovableBoxMesh.h"
#include "AMP/ampmesh/structured/PureLogicalMesh.h"
#include "AMP/ampmesh/structured/StructuredGeometryMesh.h"
#include "AMP/ampmesh/structured/structuredMeshElement.h"
#include "AMP/ampmesh/structured/structuredMeshIterator.h"
#include "AMP/utils/Utilities.h"

#ifdef USE_AMP_VECTORS
#include "AMP/discretization/simpleDOF_Manager.h"
#include "AMP/vectors/Variable.h"
#include "AMP/vectors/Vector.h"
#include "AMP/vectors/VectorBuilder.h"
#include "AMP/vectors/data/ArrayVectorData.h"
#include "AMP/vectors/operations/VectorOperationsDefault.h"
#endif

#include "ProfilerApp.h"

#include <algorithm>
#include <cstring>
#include <iostream>


namespace AMP {
namespace Mesh {


/****************************************************************
 * Generator                                                     *
 ****************************************************************/
std::shared_ptr<BoxMesh> BoxMesh::generate( std::shared_ptr<const MeshParameters> params )
{
    auto db        = params->getDatabase();
    auto generator = db->getWithDefault<std::string>( "Generator", "" );
    if ( generator == "logical" ) {
        return std::make_shared<PureLogicalMesh>( params );
    } else {
        std::shared_ptr<BoxMesh> mesh( new StructuredGeometryMesh( params ) );
        bool static_mesh = db->getWithDefault<bool>( "static", false );
        if ( !static_mesh )
            mesh.reset( new MovableBoxMesh( *mesh ) );
        return mesh;
    }
}


/****************************************************************
 * Estimate the mesh size                                        *
 ****************************************************************/
size_t BoxMesh::estimateMeshSize( std::shared_ptr<const MeshParameters> params )
{
    auto size = estimateLogicalMeshSize( params );
    size_t N  = 1;
    for ( auto s : size )
        N *= s;
    return N;
}

std::vector<size_t> BoxMesh::estimateLogicalMeshSize( std::shared_ptr<const MeshParameters> params )
{
    auto db    = params->getDatabase();
    auto geom  = AMP::Geometry::Geometry::buildGeometry( db );
    auto geom2 = std::dynamic_pointer_cast<AMP::Geometry::LogicalGeometry>( geom );
    auto size  = geom2->getLogicalGridSize( db->getVector<int>( "Size" ) );
    std::vector<size_t> N( size.size() );
    for ( size_t i = 0; i < N.size(); i++ )
        N[i] = size[i];
    return N;
}


/****************************************************************
 * Constructor                                                   *
 ****************************************************************/
BoxMesh::BoxMesh( std::shared_ptr<const MeshParameters> params_in ) : Mesh( params_in )
{
    // Check for valid inputs
    AMP_INSIST( d_params != nullptr, "Params must not be null" );
    AMP_INSIST( d_db.get(), "Database must exist" );
    d_isPeriodic.fill( false );
    d_globalSize.fill( 1 );
    d_indexSize.fill( 0 );
    d_localIndex.fill( 0 );
    d_numBlocks.fill( 1 );
    d_surfaceId.fill( -1 );
    d_rank = -1;
    d_size = 0;
    if ( !d_comm.isNull() ) {
        d_rank = d_comm.getRank();
        d_size = d_comm.getSize();
    }
}
BoxMesh::BoxMesh( const BoxMesh &mesh ) : Mesh( mesh )
{
    PhysicalDim  = mesh.PhysicalDim;
    GeomDim      = mesh.GeomDim;
    d_max_gcw    = mesh.d_max_gcw;
    d_comm       = mesh.d_comm;
    d_rank       = mesh.d_rank;
    d_size       = mesh.d_size;
    d_name       = mesh.d_name;
    d_box        = mesh.d_box;
    d_box_local  = mesh.d_box_local;
    d_isPeriodic = mesh.d_isPeriodic;
    d_globalSize = mesh.d_globalSize;
    d_numBlocks  = mesh.d_numBlocks;
    d_localIndex = mesh.d_localIndex;
    d_indexSize  = mesh.d_indexSize;
    for ( int d = 0; d < 3; d++ ) {
        d_startIndex[d] = mesh.d_startIndex[d];
        d_endIndex[d]   = mesh.d_endIndex[d];
    }
    d_surfaceId = mesh.d_surfaceId;
    for ( int d = 0; d < 4; d++ ) {
        for ( int i = 0; i < 6; i++ )
            d_globalSurfaceList[i][d] = mesh.d_globalSurfaceList[i][d];
    }
}


/****************************************************************
 * Perform the load balancing                                    *
 ****************************************************************/
void BoxMesh::loadBalance( std::array<int, 3> size,
                           int N_procs,
                           std::vector<int> *startIndex,
                           const AMP::Database *db )
{
    AMP_ASSERT( size[0] > 0 && size[1] > 0 && size[2] > 0 );
    // Check if we are dealing with a serial mesh
    if ( N_procs == 1 ) {
        startIndex[0] = std::vector<int>( 1, 0 );
        startIndex[1] = std::vector<int>( 1, 0 );
        startIndex[2] = std::vector<int>( 1, 0 );
        return;
    }
    // Get the minimum size / proc
    std::array<int, 3> minSize = { 1, 1, 1 };
    if ( db ) {
        auto tmp = db->getWithDefault<std::vector<int>>( "LoadBalanceMinSize", { 1, 1, 1 } );
        if ( tmp.size() == 1 )
            tmp.resize( 3, tmp[0] );
        tmp.resize( 3, -1 );
        minSize = { tmp[0], tmp[1], tmp[2] };
        for ( int d = 0; d < 3; d++ ) {
            if ( minSize[d] == -1 )
                minSize[d] = size[d];
        }
    }
    // Get the number of processors for each dimension
    int numBlocks[3] = { 1, 1, 1 };
    auto factors     = AMP::Utilities::factor( N_procs );
    while ( !factors.empty() ) {
        int d    = -1;
        double v = -1;
        for ( int i = 0; i < 3; i++ ) {
            double tmp = (double) size[i] / (double) numBlocks[i];
            if ( tmp > v && tmp > minSize[i] && minSize[i] >= 0 ) {
                d = i;
                v = tmp;
            }
        }
        if ( d == -1 )
            break;
        numBlocks[d] *= factors.back();
        factors.pop_back();
    }
    // Calculate the starting index for each dimension
    for ( int d = 0; d < 3; d++ ) {
        double n = size[d] / static_cast<double>( numBlocks[d] ) + 1e-12;
        startIndex[d].resize( numBlocks[d] );
        for ( int i = 0; i < numBlocks[d]; i++ )
            startIndex[d][i] = static_cast<int>( i * n );
    }
}


/****************************************************************
 * Initialize the mesh                                           *
 ****************************************************************/
void BoxMesh::initialize()
{
    PROFILE_SCOPED( timer, "initialize" );
    // Check some assumptions/variables
    AMP_INSIST( static_cast<int>( GeomDim ) <= 3, "Geometric dimension must be <= 3" );
    for ( int i = 2 * static_cast<int>( GeomDim ); i < 6; i++ )
        d_surfaceId[i] = -1;
    for ( int i = 0; i < 6; i++ ) {
        if ( d_isPeriodic[i / 2] )
            AMP_ASSERT( d_surfaceId[i] == -1 );
    }
    for ( int d = 0; d < static_cast<int>( GeomDim ); d++ )
        AMP_ASSERT( d_globalSize[d] > 0 );
    // Create the load balance
    AMP_INSIST( d_size > 0, "Communicator must be set" );
    loadBalance( d_globalSize, d_size, d_startIndex, d_db.get() );
    // Set some cached values
    for ( int d = 0; d < 3; d++ ) {
        AMP_ASSERT( !d_startIndex[d].empty() );
        d_numBlocks[d] = d_startIndex[d].size();
        d_endIndex[d].resize( d_numBlocks[d] );
        for ( int i = 1; i < d_numBlocks[d]; i++ )
            d_endIndex[d][i - 1] = d_startIndex[d][i];
        d_endIndex[d].back() = d_globalSize[d];
    }
    auto block   = getLocalBlock( d_rank );
    d_indexSize  = { block[1] - block[0] + 3, block[3] - block[2] + 3, block[5] - block[4] + 3 };
    d_localIndex = block;
    for ( int d = 0; d < 3; d++ ) {
        if ( d_localIndex[2 * d + 1] == d_globalSize[d] - 1 )
            d_localIndex[2 * d + 1] = d_globalSize[d];
        d_localIndex[2 * d + 1]++;
    }
    // Create the list of elements on each surface
    const std::array<int, 6> globalRange = { 0, std::max( d_globalSize[0] - 1, 0 ),
                                             0, std::max( d_globalSize[1] - 1, 0 ),
                                             0, std::max( d_globalSize[2] - 1, 0 ) };
    for ( int d = 0; d < static_cast<int>( GeomDim ); d++ ) {
        // Loop through the different geometry types;
        for ( int t = 0; t <= static_cast<int>( GeomDim ); t++ ) {
            auto type         = static_cast<GeomType>( t );
            auto range1       = globalRange;
            auto range2       = globalRange;
            range1[2 * d + 0] = 0;
            range1[2 * d + 1] = 0;
            range2[2 * d + 0] = d_globalSize[d] - 1;
            range2[2 * d + 1] = d_globalSize[d] - 1;
            auto set1         = getIteratorRange( range1, type, 0 );
            auto set2         = getIteratorRange( range2, type, 0 );
            // Create the surface list
            if ( type == GeomDim ) {
                AMP_ASSERT( set1.size() == 1u && set2.size() == 1u );
                d_globalSurfaceList[2 * d + 0][t].push_back( set1[0] );
                d_globalSurfaceList[2 * d + 1][t].push_back( set2[0] );
            } else if ( type == GeomType::Vertex ) {
                AMP_ASSERT( set1.size() == 1u && set2.size() == 1u );
                set1[0].first.index( d )  = 0;
                set1[0].second.index( d ) = 0;
                set2[0].first.index( d )  = d_globalSize[d];
                set2[0].second.index( d ) = d_globalSize[d];
                d_globalSurfaceList[2 * d + 0][t].push_back( set1[0] );
                d_globalSurfaceList[2 * d + 1][t].push_back( set2[0] );
            } else if ( type == GeomType::Edge && GeomDim == GeomType::Face ) {
                AMP_ASSERT( set1.size() == 2u && set2.size() == 2u );
                for ( size_t i = 0; i < set1.size(); i++ ) {
                    if ( set1[i].first.side() == d ) {
                        set1[i].first.index( d )  = 0;
                        set1[i].second.index( d ) = 0;
                        set2[i].first.index( d )  = d_globalSize[d];
                        set2[i].second.index( d ) = d_globalSize[d];
                        d_globalSurfaceList[2 * d + 0][t].push_back( set1[i] );
                        d_globalSurfaceList[2 * d + 1][t].push_back( set2[i] );
                    }
                }
            } else if ( type == GeomType::Edge && GeomDim == GeomType::Volume ) {
                AMP_ASSERT( set1.size() == 3u && set2.size() == 3u );
                for ( size_t i = 0; i < set1.size(); i++ ) {
                    if ( set1[i].first.side() == d ) {
                        set1[i].first.index( d )  = 0;
                        set1[i].second.index( d ) = 0;
                        set2[i].first.index( d )  = d_globalSize[d];
                        set2[i].second.index( d ) = d_globalSize[d];
                        d_globalSurfaceList[2 * d + 0][t].push_back( set1[i] );
                        d_globalSurfaceList[2 * d + 1][t].push_back( set2[i] );
                    }
                }
            } else if ( type == GeomType::Face && GeomDim == GeomType::Volume ) {
                AMP_ASSERT( set1.size() == 3u && set2.size() == 3u );
                for ( size_t i = 0; i < set1.size(); i++ ) {
                    if ( set1[i].first.side() == d ) {
                        set1[i].first.index( d )  = 0;
                        set1[i].second.index( d ) = 0;
                        set2[i].first.index( d )  = d_globalSize[d];
                        set2[i].second.index( d ) = d_globalSize[d];
                        d_globalSurfaceList[2 * d + 0][t].push_back( set1[i] );
                        d_globalSurfaceList[2 * d + 1][t].push_back( set2[i] );
                    }
                }
            } else {
                AMP_ERROR( "Unknown type" );
            }
        }
    }
    for ( int i = 0; i < 6; i++ ) {
        if ( d_surfaceId[i] == -1 ) {
            for ( int j = 0; j < static_cast<int>( GeomDim ); j++ )
                d_globalSurfaceList[i][j].clear();
        }
    }
    // Create the initial boundary info
}
void BoxMesh::createBoundingBox()
{
    // Fill the bounding box
    AMP_INSIST( !d_comm.isNull(), "Communicator must be set" );
    d_box_local = std::vector<double>( 2 * PhysicalDim );
    for ( int d = 0; d < PhysicalDim; d++ ) {
        d_box_local[2 * d + 0] = 1e100;
        d_box_local[2 * d + 1] = -1e100;
    }
    double x[3] = { 0, 0, 0 };
    for ( auto &node : getIterator( GeomType::Vertex, 0 ) ) {
        auto element = dynamic_cast<structuredMeshElement *>( node.getRawElement() );
        AMP_ASSERT( element != nullptr );
        coord( element->getIndex(), x );
        for ( int d = 0; d < PhysicalDim; d++ ) {
            if ( x[d] != x[d] )
                AMP_ERROR( "NaNs detected" );
            d_box_local[2 * d + 0] = std::min( d_box_local[2 * d + 0], x[d] );
            d_box_local[2 * d + 1] = std::max( d_box_local[2 * d + 1], x[d] );
        }
    }
    d_box = std::vector<double>( PhysicalDim * 2 );
    for ( int i = 0; i < PhysicalDim; i++ ) {
        d_box[2 * i + 0] = d_comm.minReduce( d_box_local[2 * i + 0] );
        d_box[2 * i + 1] = d_comm.maxReduce( d_box_local[2 * i + 1] );
    }
}
void BoxMesh::finalize()
{
    PROFILE_START( "finalize" );
    // Fill in the final info for the mesh
    AMP_INSIST( d_db->keyExists( "MeshName" ), "MeshName must exist in input database" );
    d_name = d_db->getString( "MeshName" );
    createBoundingBox();
    // Displace the mesh
    std::vector<double> displacement( PhysicalDim, 0.0 );
    if ( d_db->keyExists( "x_offset" ) && PhysicalDim >= 1 )
        displacement[0] = d_db->getScalar<double>( "x_offset" );
    if ( d_db->keyExists( "y_offset" ) && PhysicalDim >= 2 )
        displacement[1] = d_db->getScalar<double>( "y_offset" );
    if ( d_db->keyExists( "z_offset" ) && PhysicalDim >= 3 )
        displacement[2] = d_db->getScalar<double>( "z_offset" );
    bool test = false;
    for ( auto &elem : displacement ) {
        if ( elem != 0.0 )
            test = true;
    }
    if ( test )
        displaceMesh( displacement );
    PROFILE_STOP( "finalize" );
}


/****************************************************************
 * De-constructor                                                *
 ****************************************************************/
BoxMesh::~BoxMesh() = default;


/****************************************************************
 * Estimate the maximum number of processors                     *
 ****************************************************************/
size_t BoxMesh::maxProcs( std::shared_ptr<const MeshParameters> params )
{
    // Check for valid inputs
    AMP_INSIST( params.get(), "Params must not be null" );
    auto db = params->getDatabase();
    AMP_INSIST( db.get(), "Database must exist" );
    size_t maxProcs = 1;
    if ( db->keyExists( "LoadBalanceMinSize" ) ) {
        auto minSize  = db->getVector<int>( "LoadBalanceMinSize" );
        auto meshSize = estimateLogicalMeshSize( params );
        for ( size_t i = 0; i < meshSize.size(); i++ ) {
            if ( minSize[i] == 0 )
                minSize[i] = 1;
            if ( minSize[i] != -1 )
                maxProcs *= ( meshSize[i] / minSize[i] );
        }
    } else {
        maxProcs = estimateMeshSize( params );
    }
    return maxProcs;
}


/****************************************************************
 * Function to return the element given an ID                    *
 ****************************************************************/
MeshElement BoxMesh::getElement( const MeshElementID &id ) const
{
    // Get the index of the element
    MeshElementIndex index = convert( id );
    // Create the element
    structuredMeshElement elem( index, this );
    AMP_ASSERT( elem.globalID() == id );
    return std::move( elem );
}
MeshElement BoxMesh::getElement( const MeshElementIndex &index ) const
{
    return structuredMeshElement( index, this );
}


/****************************************************************
 * Find the mesh element index from a point                      *
 ****************************************************************/
static inline int to_nearest( double x ) { return static_cast<int>( floor( x + 0.5 ) ); }
BoxMesh::MeshElementIndex BoxMesh::getElementFromLogical( const AMP::Geometry::Point &x0,
                                                          GeomType type ) const
{
    // Correct x for periodic boundaries
    double x[3] = { x0.x(), x0.y(), x0.z() };
    for ( int d = 0; d < static_cast<int>( GeomDim ); d++ ) {
        if ( d_isPeriodic[d] ) {
            while ( x[d] < 0 )
                x[d] += 1.0;
            while ( x[d] >= 1.0 )
                x[d] -= 1.0;
        }
    }
    // Convert x to [0,size]
    x[0] = x[0] * d_globalSize[0];
    x[1] = x[1] * d_globalSize[1];
    x[2] = x[2] * d_globalSize[2];
    // Check if element is outside domain
    for ( int d = 0; d < static_cast<int>( GeomDim ); d++ ) {
        if ( fabs( x[d] ) < 1e-6 )
            x[d] = 0;
        if ( fabs( x[d] - d_globalSize[d] ) < 1e-6 )
            x[d] = d_globalSize[d];
        if ( x[d] < 0 || x[d] > d_globalSize[d] )
            return MeshElementIndex();
    }
    // Compute the index
    MeshElementIndex index;
    if ( type == GeomDim ) {
        index = MeshElementIndex( GeomDim, 0, x[0], x[1], x[2] );
    } else if ( type == GeomType::Vertex ) {
        int i     = to_nearest( x[0] );
        int j     = to_nearest( x[1] );
        int k     = to_nearest( x[2] );
        bool keep = fabs( x[0] - i ) < 1e-6 && fabs( x[1] - j ) < 1e-6 && fabs( x[2] - k ) < 1e-6;
        keep      = keep && i >= 0 && j >= 0 && k >= 0;
        keep      = keep && i <= d_globalSize[0] && j <= d_globalSize[1] && k <= d_globalSize[2];
        if ( keep )
            index = MeshElementIndex( GeomType::Vertex, 0, i, j, k );
    } else if ( type == GeomType::Edge ) {
        AMP_ERROR( "Not finished" );
    } else if ( type == GeomType::Face ) {
        int i      = to_nearest( x[0] );
        int j      = to_nearest( x[1] );
        int k      = to_nearest( x[2] );
        int ijk    = 0;
        double min = fabs( x[0] - i );
        if ( fabs( x[1] - j ) < min && static_cast<int>( GeomDim ) >= 2 ) {
            min = fabs( x[1] - j );
            ijk = 1;
        }
        if ( fabs( x[2] - k ) < min && static_cast<int>( GeomDim ) >= 3 ) {
            min = fabs( x[2] - k );
            ijk = 2;
        }
        if ( min > 1e-6 ) {
            // Point is not on any face
        } else if ( ijk == 0 ) {
            index = MeshElementIndex( GeomType::Face, 0, i, x[1], x[2] );
        } else if ( ijk == 1 ) {
            index = MeshElementIndex( GeomType::Face, 1, x[0], j, x[2] );
        } else if ( ijk == 2 ) {
            index = MeshElementIndex( GeomType::Face, 2, x[0], x[1], k );
        }
    } else if ( type == GeomType::Volume ) {
        AMP_ERROR( "Not finished" );
    } else {
        AMP_ERROR( "Unknown mesh element type" );
    }
    return index;
}
BoxMesh::MeshElementIndex BoxMesh::getElementFromPhysical( const AMP::Geometry::Point &x,
                                                           GeomType type ) const
{
    auto logical = physicalToLogical( x );
    auto index   = getElementFromLogical( logical, type );
    return index;
}


/********************************************************
 * Function to return parents of an element              *
 ********************************************************/
std::vector<MeshElement> BoxMesh::getElementParents( const MeshElement &meshelem,
                                                     const GeomType type ) const
{
    auto id = meshelem.globalID();
    if ( type == id.type() )
        return std::vector<MeshElement>( 1, meshelem );
    AMP_INSIST( id.meshID() == d_meshID, "MeshElement is not from the given mesh" );
    // AMP_INSIST( type >= id.type() && type <= GeomDim,
    //            "Cannot get the parents of the given type for the current element" );
    // Get the element of interest
    const auto *elem = dynamic_cast<const structuredMeshElement *>( meshelem.getRawElement() );
    AMP_ASSERT( elem != nullptr );
    return elem->getParents( type );
}


/****************************************************************
 * Functions to return the number of elements                    *
 ****************************************************************/
size_t BoxMesh::numLocalElements( const GeomType type ) const
{
    if ( type > GeomDim )
        return 0;
    auto box   = getLocalBlock( d_rank );
    auto range = getIteratorRange( box, type, 0 );
    size_t N   = 0;
    for ( const auto &tmp : range )
        N += BoxMesh::MeshElementIndex::numElements( tmp.first, tmp.second );
    return N;
}
size_t BoxMesh::numGlobalElements( const GeomType type ) const
{
    if ( type > GeomDim )
        return 0;
    std::array<int, 6> box = {
        { 0, d_globalSize[0] - 1, 0, d_globalSize[1] - 1, 0, d_globalSize[2] - 1 }
    };
    auto range = getIteratorRange( box, type, 0 );
    size_t N   = 0;
    for ( auto tmp : range )
        N += MeshElementIndex::numElements( tmp.first, tmp.second );
    return N;
}
size_t BoxMesh::numGhostElements( const GeomType type, int gcw ) const
{
    if ( type > GeomDim )
        return 0;
    auto box    = getLocalBlock( d_rank );
    auto range1 = getIteratorRange( box, type, 0 );
    auto range2 = getIteratorRange( box, type, gcw );
    size_t N    = 0;
    for ( size_t i = 0; i < range1.size(); i++ ) {
        size_t N1 = BoxMesh::MeshElementIndex::numElements( range1[i].first, range1[i].second );
        size_t N2 = BoxMesh::MeshElementIndex::numElements( range2[i].first, range2[i].second );
        N += N2 - N1;
    }
    return N;
}


/****************************************************************
 * Function to get an iterator                                   *
 ****************************************************************/
BoxMesh::ElementBlocks
BoxMesh::getIteratorRange( std::array<int, 6> range, const GeomType type, const int gcw ) const
{
    AMP_ASSERT( type <= GeomDim );
    // Get the range of cells we care about
    for ( int d = 0; d < static_cast<int>( GeomDim ); d++ ) {
        range[2 * d + 0] -= gcw;
        range[2 * d + 1] += gcw;
        if ( !d_isPeriodic[d] ) {
            range[2 * d + 0] = std::max( range[2 * d + 0], 0 );
            range[2 * d + 1] = std::min( range[2 * d + 1], d_globalSize[d] - 1 );
        }
    }
    // Get the element blocks we want to process
    ElementBlocks blocks;
    blocks.reserve( 3 );
    if ( type == GeomDim ) {
        blocks.emplace_back( MeshElementIndex( type, 0, range[0], range[2], range[4] ),
                             MeshElementIndex( type, 0, range[1], range[3], range[5] ) );
    } else if ( type == GeomType::Vertex ) {
        for ( int d = 0; d < static_cast<int>( GeomDim ); d++ ) {
            if ( gcw != 0 )
                range[2 * d + 1]++;
            else if ( !d_isPeriodic[d] && range[2 * d + 1] == d_globalSize[d] - 1 )
                range[2 * d + 1]++;
        }
        blocks.emplace_back( MeshElementIndex( type, 0, range[0], range[2], range[4] ),
                             MeshElementIndex( type, 0, range[1], range[3], range[5] ) );
    } else if ( type == GeomType::Edge && GeomDim == GeomType::Face ) {
        auto range1 = range;
        auto range2 = range;
        if ( gcw != 0 || ( gcw == 0 && !d_isPeriodic[0] && range[1] == d_globalSize[0] - 1 ) )
            range2[1]++;
        if ( gcw != 0 || ( gcw == 0 && !d_isPeriodic[1] && range[3] == d_globalSize[1] - 1 ) )
            range1[3]++;
        blocks.emplace_back( MeshElementIndex( type, 0, range1[0], range1[2], range1[4] ),
                             MeshElementIndex( type, 0, range1[1], range1[3], range1[5] ) );
        blocks.emplace_back( MeshElementIndex( type, 1, range2[0], range2[2], range2[4] ),
                             MeshElementIndex( type, 1, range2[1], range2[3], range2[5] ) );
    } else if ( type == GeomType::Edge && GeomDim == GeomType::Volume ) {
        auto range1 = range;
        auto range2 = range;
        auto range3 = range;
        if ( gcw != 0 || ( gcw == 0 && !d_isPeriodic[0] && range[1] == d_globalSize[0] - 1 ) ) {
            range2[1]++;
            range3[1]++;
        }
        if ( gcw != 0 || ( gcw == 0 && !d_isPeriodic[1] && range[3] == d_globalSize[1] - 1 ) ) {
            range1[3]++;
            range3[3]++;
        }
        if ( gcw != 0 || ( gcw == 0 && !d_isPeriodic[2] && range[5] == d_globalSize[2] - 1 ) ) {
            range1[5]++;
            range2[5]++;
        }
        blocks.emplace_back( MeshElementIndex( type, 0, range1[0], range1[2], range1[4] ),
                             MeshElementIndex( type, 0, range1[1], range1[3], range1[5] ) );
        blocks.emplace_back( MeshElementIndex( type, 1, range2[0], range2[2], range2[4] ),
                             MeshElementIndex( type, 1, range2[1], range2[3], range2[5] ) );
        blocks.emplace_back( MeshElementIndex( type, 2, range3[0], range3[2], range3[4] ),
                             MeshElementIndex( type, 2, range3[1], range3[3], range3[5] ) );
    } else if ( type == GeomType::Face && GeomDim == GeomType::Volume ) {
        auto range1 = range;
        auto range2 = range;
        auto range3 = range;
        if ( gcw != 0 || ( gcw == 0 && !d_isPeriodic[0] && range[1] == d_globalSize[0] - 1 ) )
            range1[1]++;
        if ( gcw != 0 || ( gcw == 0 && !d_isPeriodic[1] && range[3] == d_globalSize[1] - 1 ) )
            range2[3]++;
        if ( gcw != 0 || ( gcw == 0 && !d_isPeriodic[2] && range[5] == d_globalSize[2] - 1 ) )
            range3[5]++;
        blocks.emplace_back( MeshElementIndex( type, 0, range1[0], range1[2], range1[4] ),
                             MeshElementIndex( type, 0, range1[1], range1[3], range1[5] ) );
        blocks.emplace_back( MeshElementIndex( type, 1, range2[0], range2[2], range2[4] ),
                             MeshElementIndex( type, 1, range2[1], range2[3], range2[5] ) );
        blocks.emplace_back( MeshElementIndex( type, 2, range3[0], range3[2], range3[4] ),
                             MeshElementIndex( type, 2, range3[1], range3[3], range3[5] ) );
    } else {
        AMP_ERROR( "Unknown case" );
    }
    // Check that each block does not have duplicate elements
    for ( auto &block : blocks ) {
        for ( int d = 0; d < static_cast<int>( GeomDim ); d++ ) {
            if ( d_isPeriodic[d] ) {
                auto &first = block.first;
                auto &last  = block.second;
                if ( first.index( d ) + d_globalSize[d] <= last.index( d ) ) {
                    first.index( d ) = 0;
                    last.index( d )  = d_globalSize[d] - 1;
                }
            }
        }
    }
    return blocks;
}
BoxMesh::ElementBlocks BoxMesh::intersect( const ElementBlocks &set1, const ElementBlocks &set2 )
{
    ElementBlocks set;
    set.reserve( set1.size() * set2.size() );
    for ( const auto &v1 : set1 ) {
        for ( const auto &v2 : set2 ) {
            if ( v1.first.type() != v2.first.type() || v1.first.side() != v2.first.side() )
                continue;
            auto v              = v1;
            v.first.index( 0 )  = std::max( v1.first.index( 0 ), v2.first.index( 0 ) );
            v.first.index( 1 )  = std::max( v1.first.index( 1 ), v2.first.index( 1 ) );
            v.first.index( 2 )  = std::max( v1.first.index( 2 ), v2.first.index( 2 ) );
            v.second.index( 0 ) = std::min( v1.second.index( 0 ), v2.second.index( 0 ) );
            v.second.index( 1 ) = std::min( v1.second.index( 1 ), v2.second.index( 1 ) );
            v.second.index( 2 ) = std::min( v1.second.index( 2 ), v2.second.index( 2 ) );
            if ( MeshElementIndex::numElements( v.first, v.second ) > 0 )
                set.push_back( v );
        }
    }
    return set;
}
inline MeshIterator BoxMesh::createIterator( const ElementBlocks &list ) const
{
    if ( list.empty() ) {
        return MeshIterator();
    } else if ( list.size() == 1 ) {
        return structuredMeshIterator( list[0].first, list[0].second, this, 0 );
    } else {
        std::vector<MeshIterator> iterator_list;
        iterator_list.reserve( list.size() );
        for ( const auto &item : list ) {
            if ( MeshElementIndex::numElements( item.first, item.second ) ) {
                structuredMeshIterator it( item.first, item.second, this, 0 );
                iterator_list.push_back( it );
            }
        }
        return MultiIterator( iterator_list, 0 );
    }
    return MeshIterator();
}
MeshIterator BoxMesh::getIterator( const GeomType type, const int gcw ) const
{
    if ( type > GeomDim )
        return MeshIterator();
    auto box   = getLocalBlock( d_rank );
    auto range = getIteratorRange( box, type, gcw );
    return createIterator( range );
}


/****************************************************************
 * Function to get an iterator over the surface                  *
 ****************************************************************/
MeshIterator BoxMesh::getSurfaceIterator( const GeomType type, const int gcw ) const
{
    if ( type > GeomDim )
        return MeshIterator();
    // Include each surface as needed
    ElementBlocks sufaceSet;
    for ( int i = 0; i < 2 * static_cast<int>( GeomDim ); i++ ) {
        if ( d_surfaceId[i] != -1 ) {
            for ( const auto &item : d_globalSurfaceList[i][static_cast<int>( type )] )
                sufaceSet.emplace_back( item );
        }
    }
    // Intersect with the local ghost box
    auto box          = getLocalBlock( d_rank );
    auto range        = getIteratorRange( box, type, gcw );
    auto intersection = intersect( sufaceSet, range );
    // Create a list if elements removing any duplicate elements
    std::set<MeshElementIndex> set;
    for ( auto block : intersection ) {
        auto itype = block.first.type();
        auto side  = block.first.side();
        for ( int k = block.second.index( 2 ); k >= block.first.index( 2 ); k-- ) {
            for ( int j = block.second.index( 1 ); j >= block.first.index( 1 ); j-- ) {
                for ( int i = block.second.index( 0 ); i >= block.first.index( 0 ); i-- ) {
                    set.emplace( MeshElementIndex( itype, side, i, j, k ) );
                }
            }
        }
    }
    // Create the iterator
    std::shared_ptr<std::vector<MeshElementIndex>> elements(
        new std::vector<MeshElementIndex>( set.begin(), set.end() ) );
    return structuredMeshIterator( elements, this, 0 );
}


/****************************************************************
 * Functions to get the boundaries                               *
 ****************************************************************/
std::vector<int> BoxMesh::getBoundaryIDs() const
{
    std::set<int> ids;
    for ( int i = 0; i < 2 * static_cast<int>( GeomDim ); i++ ) {
        if ( !d_isPeriodic[i / 2] && d_surfaceId[i] != -1 )
            ids.insert( d_surfaceId[i] );
    }
    return std::vector<int>( ids.begin(), ids.end() );
}
MeshIterator
BoxMesh::getBoundaryIDIterator( const GeomType type, const int id, const int gcw ) const
{
    if ( type > GeomDim )
        return MeshIterator();
    // Include each surface as needed
    ElementBlocks sufaceSet;
    for ( int i = 0; i < 2 * static_cast<int>( GeomDim ); i++ ) {
        if ( d_surfaceId[i] == id ) {
            for ( auto item : d_globalSurfaceList[i][static_cast<int>( type )] )
                sufaceSet.emplace_back( item );
        }
    }
    // Intersect with the local ghost box
    auto box          = getLocalBlock( d_rank );
    auto range        = getIteratorRange( box, type, gcw );
    auto intersection = intersect( sufaceSet, range );
    // Create a list if elements removing any duplicate elements
    std::set<MeshElementIndex> set;
    for ( auto block : intersection ) {
        auto itype = block.first.type();
        auto side  = block.first.side();
        for ( int k = block.second.index( 2 ); k >= block.first.index( 2 ); k-- ) {
            for ( int j = block.second.index( 1 ); j >= block.first.index( 1 ); j-- ) {
                for ( int i = block.second.index( 0 ); i >= block.first.index( 0 ); i-- ) {
                    set.emplace( MeshElementIndex( itype, side, i, j, k ) );
                }
            }
        }
    }
    // Create the iterator
    std::shared_ptr<std::vector<MeshElementIndex>> elements(
        new std::vector<MeshElementIndex>( set.begin(), set.end() ) );
    return structuredMeshIterator( elements, this, 0 );
}
std::vector<int> BoxMesh::getBlockIDs() const { return std::vector<int>( 1, 0 ); }
MeshIterator BoxMesh::getBlockIDIterator( const GeomType type, const int id, const int gcw ) const
{
    if ( id == 0 )
        return getIterator( type, gcw );
    return MeshIterator();
}
bool BoxMesh::isOnBoundary( const MeshElementIndex &index, int id ) const
{
    bool test = false;
    for ( int i = 0; i < 6; i++ ) {
        if ( d_surfaceId[i] == id ) {
            int d = i / 2;
            int s = i % 2;
            if ( index.type() == GeomDim ) {
                test |= ( s == 0 && index.index( d ) == 0 );
                test |= ( s == 1 && index.index( d ) == d_globalSize[d] - 1 );
            } else if ( index.type() == GeomType::Vertex ) {
                test |= ( s == 0 && index.index( d ) == 0 );
                test |= ( s == 1 && index.index( d ) == d_globalSize[d] );
            } else if ( index.type() == GeomType::Edge ) {
                test |= ( s == 0 && index.side() == d && index.index( d ) == 0 );
                test |= ( s == 1 && index.side() == d && index.index( d ) == d_globalSize[d] );
            } else if ( index.type() == GeomType::Face ) {
                test |= ( s == 0 && index.side() == d && index.index( d ) == 0 );
                test |= ( s == 1 && index.side() == d && index.index( d ) == d_globalSize[d] );
            } else {
                AMP_ERROR( "Unknown type" );
            }
        }
    }
    return test;
}


/****************************************************************
 * Create an ArrayVector over the mesh                           *
 ****************************************************************/
#ifdef USE_AMP_VECTORS
std::shared_ptr<AMP::LinearAlgebra::Vector> BoxMesh::createVector( const std::string &name,
                                                                   int gcw )
{
    AMP_ASSERT( getComm().getSize() == 1 );
    auto mesh = shared_from_this();
    auto size = getLocalBox().size();
    auto type = GeomDim;
    auto var  = std::make_shared<AMP::LinearAlgebra::Variable>( name );
    auto DOFs = AMP::Discretization::simpleDOFManager::create( mesh, type, gcw, 1, true );
    auto ops  = std::make_shared<AMP::LinearAlgebra::VectorOperationsDefault<double>>();
    auto data = AMP::LinearAlgebra::ArrayVectorData<double>::create( size );
    auto vec  = std::make_shared<AMP::LinearAlgebra::Vector>( data, ops, var, DOFs );
    AMP_ASSERT( vec->getLocalSize() == getIterator( type ).size() );
    return vec;
}
#endif


/****************************************************************
 * Check if two meshes are equal                                 *
 ****************************************************************/
bool BoxMesh::operator==( const Mesh &rhs ) const
{
    // Check if we are comparing to *this
    if ( &rhs == this )
        return true;
    // Check if we can cast to a BoxMesh
    auto mesh = dynamic_cast<const BoxMesh *>( &rhs );
    if ( !mesh )
        return false;
    // Perform basic comparison
    if ( d_isPeriodic != mesh->d_isPeriodic || d_globalSize != mesh->d_globalSize )
        return false;
    if ( d_numBlocks != mesh->d_numBlocks || d_indexSize != mesh->d_indexSize ||
         d_localIndex != mesh->d_localIndex )
        return false;
    if ( d_surfaceId != mesh->d_surfaceId || *d_geometry != *mesh->d_geometry )
        return false;
    for ( int d = 0; d < 3; d++ ) {
        if ( d_startIndex[d] != mesh->d_startIndex[d] || d_endIndex[d] != mesh->d_endIndex[d] )
            return false;
    }
    return true;
}


/****************************************************************
 * Print the index                                               *
 ****************************************************************/
std::ostream &operator<<( std::ostream &out, const BoxMesh::MeshElementIndex &x )
{
    const char *type[] = { "Vertex", "Edge", "Face", "Volume" };
    char tmp[128];
    sprintf( tmp,
             "(%i,%i,%i,%s,%i)",
             x.index( 0 ),
             x.index( 1 ),
             x.index( 2 ),
             type[static_cast<int>( x.type() )],
             x.side() );
    out << tmp;
    return out;
}


} // namespace Mesh
} // namespace AMP
