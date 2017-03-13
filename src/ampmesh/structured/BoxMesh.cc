#include "ampmesh/structured/BoxMesh.h"
#include "ampmesh/structured/CubeMesh.h"
#include "ampmesh/structured/TubeMesh.h"
#include "ampmesh/structured/CircleMesh.h"
#include "ampmesh/structured/CylinderMesh.h"
#include "ampmesh/structured/ShellMesh.h"
#include "ampmesh/structured/SphereMesh.h"
#include "ampmesh/structured/MovableBoxMesh.h"

#include "ampmesh/MultiIterator.h"
#include "ampmesh/structured/structuredMeshElement.h"
#include "ampmesh/structured/structuredMeshIterator.h"


#include "utils/Utilities.h"
#ifdef USE_AMP_VECTORS
#include "vectors/Variable.h"
#include "vectors/Vector.h"
#include "vectors/VectorBuilder.h"
#endif

#include "ProfilerApp.h"

#include <algorithm>
#include <iostream>
#include <string.h>


namespace AMP {
namespace Mesh {

static inline double round( double x ) { return x < 0.0 ? ceil( x - 0.5 ) : floor( x + 0.5 ); }


/****************************************************************
* Generator                                                     *
****************************************************************/
AMP::shared_ptr<BoxMesh> BoxMesh::generate( MeshParameters::shared_ptr params )
{
    auto db = params->getDatabase();
    std::string generator = db->getString( "Generator" );
    bool static_mesh = db->getBoolWithDefault( "static", false );
    AMP::shared_ptr<BoxMesh> mesh;
    if ( generator.compare( "cube" ) == 0 ) {
        mesh.reset( new CubeMesh( params ) );
    } else if ( generator.compare( "tube" ) == 0 ) {
        mesh.reset( new TubeMesh( params ) );
    } else if ( generator.compare( "circle" ) == 0 ) {
        mesh.reset( new CircleMesh( params ) );
    } else if ( generator.compare( "cylinder" ) == 0 ) {
        mesh.reset( new CylinderMesh( params ) );
    } else if ( generator.compare( "shell" ) == 0 ) {
        mesh.reset( new ShellMesh( params ) );
    } else if ( generator.compare( "sphere" ) == 0 ) {
        mesh.reset( new SphereMesh( params ) );
    } else {
        AMP_ERROR( "Unknown generator" );
    }
    if ( !static_mesh )
        mesh.reset( new MovableBoxMesh( *mesh ) );
    return mesh;
}


/****************************************************************
* Estimate the mesh size                                        *
****************************************************************/
size_t BoxMesh::estimateMeshSize( const MeshParameters::shared_ptr &params )
{
    auto size = estimateLogicalMeshSize( params );
    size_t N = 1;
    for ( auto s : size )
        N *= s;
    return N;
}
std::vector<size_t> BoxMesh::estimateLogicalMeshSize( const MeshParameters::shared_ptr &params )
{
    auto db = params->getDatabase();
    std::string generator = db->getString( "Generator" );
    std::vector<size_t> N;
    if ( generator.compare( "cube" ) == 0 ) {
        N = CubeMesh::estimateLogicalMeshSize( params );
    } else if ( generator.compare( "tube" ) == 0 ) {
        N = TubeMesh::estimateLogicalMeshSize( params );
    } else if ( generator.compare( "circle" ) == 0 ) {
        N = CircleMesh::estimateLogicalMeshSize( params );
    } else if ( generator.compare( "cylinder" ) == 0 ) {
        N = CylinderMesh::estimateLogicalMeshSize( params );
    } else if ( generator.compare( "shell" ) == 0 ) {
        N = ShellMesh::estimateLogicalMeshSize( params );
    } else if ( generator.compare( "sphere" ) == 0 ) {
        N = SphereMesh::estimateLogicalMeshSize( params );
    } else {
        AMP_ERROR( "Unknown generator" );
    }
    return N;
}


/****************************************************************
* Constructor                                                   *
****************************************************************/
BoxMesh::BoxMesh( MeshParameters::shared_ptr params_in ):
    Mesh( params_in )
{
    // Check for valid inputs
    AMP_INSIST( d_params != nullptr, "Params must not be null" );
    AMP_INSIST( !d_comm.isNull(), "Communicator must be set" );
    AMP_INSIST( d_db.get(), "Database must exist" );
    for (int i=0; i<3; i++) {
        d_isPeriodic[i] = false;
        d_globalSize[i] = 0;
        d_numBlocks[i] = 0;
    }
    for (int i=0; i<6; i++) {
        d_surfaceId[i] = i;
        d_onSurface[i] = true;
    }
}
BoxMesh::BoxMesh( const BoxMesh &mesh ):
    Mesh( mesh.d_params )
{
    PhysicalDim = mesh.PhysicalDim;
    GeomDim = mesh.GeomDim;
    d_max_gcw = mesh.d_max_gcw;
    d_comm = mesh.d_comm;
    d_name = mesh.d_name;
    d_box = mesh.d_box;
    d_box_local = mesh.d_box_local;
    for (int d=0; d<3; d++) {
        d_isPeriodic[d] = mesh.d_isPeriodic[d];
        d_globalSize[d] = mesh.d_globalSize[d];
        d_numBlocks[d] = mesh.d_numBlocks[d];
    }
    for (int i=0; i<6; i++) {
        d_surfaceId[i] = mesh.d_surfaceId[i];
        d_onSurface[i] = mesh.d_onSurface[i];
    }
    for (int d=0; d<4; d++) {
        for (int i=0; i<6; i++)
            d_globalSurfaceList[i][d] = mesh.d_globalSurfaceList[i][d];
    }
}


/****************************************************************
* Initialize the mesh                                           *
****************************************************************/
void BoxMesh::initialize()
{
    PROFILE_START( "initialize" );
    // Check some assumptions/variables
    AMP_INSIST( GeomDim <= 3, "Geometric dimension must be <= 3" );
    for (int i=2*GeomDim; i<6; i++) {
        d_onSurface[i] = false;
        d_surfaceId[i] = -1;
    }
    for (int i=0; i<6; i++) {
        if ( d_isPeriodic[i/2] ) {
            AMP_ASSERT(d_onSurface[i]==false);
            AMP_ASSERT(d_surfaceId[i]==-1);
        } else {
            if ( !d_onSurface[i] )
                AMP_ASSERT(d_surfaceId[i]==-1);
        }
    }
    for (int d=0; d<GeomDim; d++)
        AMP_ASSERT(d_globalSize[d]>0);
    // Get the minimum mesh size
    std::vector<int> minSize( GeomDim, 1 );
    if ( d_db->keyExists( "LoadBalanceMinSize" ) ) {
        minSize = d_db->getIntegerArray( "LoadBalanceMinSize" );
        AMP_ASSERT( (int) minSize.size() == (int) GeomDim );
        for (int d=0; d<GeomDim; d++) {
            if ( minSize[d] == -1 )
                minSize[d] = d_globalSize[d];
            if ( minSize[d] == 0 )
                minSize[d] = 1;
        }
    }
    // Create the load balance
    if ( d_comm.getSize() == 1 ) {
        // We are dealing with a serial mesh (do nothing to change the local box sizes)
        for ( int d = 0; d < GeomDim; d++ )
            d_numBlocks[d] = 1;
    } else {
        // We are dealing with a parallel mesh
        // First, get the prime factors for number of processors and divide the dimensions
        std::vector<int> factors = AMP::Utilities::factor( d_comm.getSize() );
        for (int d=0; d<3; d++)
            d_numBlocks[d] = 1;
        while ( !factors.empty() ) {
            int d    = -1;
            double v = -1;
            for ( int i = 0; i < GeomDim; i++ ) {
                double tmp = (double) d_globalSize[i] / (double) d_numBlocks[i];
                if ( tmp > v && tmp > minSize[i] && minSize[i] >= 0 ) {
                    d = i;
                    v = tmp;
                }
            }
            if ( d == -1 )
                break;
            d_numBlocks[d] *= factors[factors.size() - 1];
            factors.resize( factors.size() - 1 );
        }
    }
    d_localBlock = getLocalBlock( d_comm.getRank() );
    // Create the list of elements on each surface
    const std::array<int,6> globalRange = { 0, std::max(d_globalSize[0]-1,0),
                                            0, std::max(d_globalSize[1]-1,0),
                                            0, std::max(d_globalSize[2]-1,0) };
    for ( int d = 0; d < GeomDim; d++ ) {
        // Loop through the different geometry types;
        for ( int t = 0; t <= GeomDim; t++ ) {
            GeomType type = static_cast<GeomType>(t);
            auto range1 = globalRange;
            auto range2 = globalRange;
            range1[2*d+0] = 0;
            range1[2*d+1] = 0;
            range2[2*d+0] = d_globalSize[d]-1;
            range2[2*d+1] = d_globalSize[d]-1;
            auto set1 = getIteratorRange( range1, type, 0 );
            auto set2 = getIteratorRange( range2, type, 0 );
            // Create the surface list
            if ( type == GeomDim ) {
                AMP_ASSERT( set1.size()==1u && set2.size()==1u );
                d_globalSurfaceList[2*d+0][t].push_back( set1[0] );
                d_globalSurfaceList[2*d+1][t].push_back( set2[0] );
            } else if ( type == Vertex ) {
                AMP_ASSERT( set1.size()==1u && set2.size()==1u );
                set1[0].first.index(d)  = 0;
                set1[0].second.index(d) = 0;
                set2[0].first.index(d)  = d_globalSize[d];
                set2[0].second.index(d) = d_globalSize[d];
                d_globalSurfaceList[2*d+0][t].push_back( set1[0] );
                d_globalSurfaceList[2*d+1][t].push_back( set2[0] );
            } else if ( type==Edge && GeomDim==Face ) {
                AMP_ASSERT( set1.size()==2u && set2.size()==2u );
                for (size_t i=0; i<set1.size(); i++) {
                    if ( set1[i].first.side() == d ) {
                        set1[i].first.index(d)  = 0;
                        set1[i].second.index(d) = 0;
                        set2[i].first.index(d)  = d_globalSize[d];
                        set2[i].second.index(d) = d_globalSize[d];
                        d_globalSurfaceList[2*d+0][t].push_back( set1[i] );
                        d_globalSurfaceList[2*d+1][t].push_back( set2[i] );
                    }
                }
            } else if ( type==Edge && GeomDim==Volume ) {
                AMP_ASSERT( set1.size()==3u && set2.size()==3u );
                for (size_t i=0; i<set1.size(); i++) {
                    if ( set1[i].first.side() == d ) {
                        set1[i].first.index(d)  = 0;
                        set1[i].second.index(d) = 0;
                        set2[i].first.index(d)  = d_globalSize[d];
                        set2[i].second.index(d) = d_globalSize[d];
                        d_globalSurfaceList[2*d+0][t].push_back( set1[i] );
                        d_globalSurfaceList[2*d+1][t].push_back( set2[i] );
                    }
                }
            } else if ( type==Face && GeomDim==Volume ) {
                AMP_ASSERT( set1.size()==3u && set2.size()==3u );
                for (size_t i=0; i<set1.size(); i++) {
                    if ( set1[i].first.side() == d ) {
                        set1[i].first.index(d)  = 0;
                        set1[i].second.index(d) = 0;
                        set2[i].first.index(d)  = d_globalSize[d];
                        set2[i].second.index(d) = d_globalSize[d];
                        d_globalSurfaceList[2*d+0][t].push_back( set1[i] );
                        d_globalSurfaceList[2*d+1][t].push_back( set2[i] );
                    }
                }
            } else {
                AMP_ERROR("Unknown type");
            }
        }
    }
    for (int i=0; i<6; i++) {
        if ( !d_onSurface[i] ) {
            for (int j=0; j<GeomDim; j++)
                d_globalSurfaceList[i][j].clear();
        }
    }
    // Create the initial boundary info
    PROFILE_STOP( "initialize" );
}
void BoxMesh::finalize()
{
    PROFILE_START( "finalize" );
    // Fill in the final info for the mesh
    AMP_INSIST( d_db->keyExists( "MeshName" ), "MeshName must exist in input database" );
    d_name      = d_db->getString( "MeshName" );
    d_box_local = std::vector<double>( 2 * PhysicalDim );
    for ( int d = 0; d < PhysicalDim; d++ ) {
        d_box_local[2 * d + 0] = 1e100;
        d_box_local[2 * d + 1] = -1e100;
    }
    double x[3] = {0,0,0};
    for ( auto& node : getIterator(Vertex,0) ) {
        auto element = dynamic_cast<structuredMeshElement*>( node.getRawElement() );
        AMP_ASSERT( element != nullptr );
        coord( element->getIndex(), x );
        for ( int d = 0; d < PhysicalDim; d++ ) {
            if ( x[d] != x[d] )
                AMP_ERROR( "NaNs detected" );
            d_box_local[2*d+0] = std::min(d_box_local[2*d+0],x[d]);
            d_box_local[2*d+1] = std::max(d_box_local[2*d+1],x[d]);
        }
    }
    d_box = std::vector<double>( PhysicalDim * 2 );
    for ( int i = 0; i < PhysicalDim; i++ ) {
        d_box[2 * i + 0] = d_comm.minReduce( d_box_local[2 * i + 0] );
        d_box[2 * i + 1] = d_comm.maxReduce( d_box_local[2 * i + 1] );
    }
    // Displace the mesh
    std::vector<double> displacement( PhysicalDim, 0.0 );
    if ( d_db->keyExists( "x_offset" ) && PhysicalDim >= 1 )
        displacement[0] = d_db->getDouble( "x_offset" );
    if ( d_db->keyExists( "y_offset" ) && PhysicalDim >= 2 )
        displacement[1] = d_db->getDouble( "y_offset" );
    if ( d_db->keyExists( "z_offset" ) && PhysicalDim >= 3 )
        displacement[2] = d_db->getDouble( "z_offset" );
    bool test           = false;
    for ( auto &elem : displacement ) {
        if ( elem != 0.0 )
            test = true;
    }
    if ( test )
        displaceMesh( displacement );
    // Get the global ranks for the comm to make sure it is set
    auto globalRanks = getComm().globalRanks();
    AMP_ASSERT( !globalRanks.empty() );
    PROFILE_STOP( "finalize" );
}


/****************************************************************
* De-constructor                                                *
****************************************************************/
BoxMesh::~BoxMesh() {}


/****************************************************************
* Estimate the maximum number of processors                     *
****************************************************************/
size_t BoxMesh::maxProcs( const MeshParameters::shared_ptr &params )
{
    // Check for valid inputs
    AMP_INSIST( params.get(), "Params must not be null" );
    AMP::shared_ptr<AMP::Database> db = params->getDatabase();
    AMP_INSIST( db.get(), "Database must exist" );
    size_t maxProcs = 1;
    if ( db->keyExists( "LoadBalanceMinSize" ) ) {
        auto minSize  = db->getIntegerArray( "LoadBalanceMinSize" );
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
MeshElement BoxMesh::getElement( const MeshElementID &elem_id ) const
{
    auto range = getLocalBlock( elem_id.owner_rank() );
    // Increase the index range for the boxes on the boundary for all elements except the current
    // dimension
    if ( elem_id.type() != GeomDim ) {
        for ( int d = 0; d < GeomDim; d++ ) {
            if ( range[2 * d + 1] == d_globalSize[d]-1 && !d_isPeriodic[d] )
                range[2 * d + 1]++;
        }
    }
    // Get the 3-index from the local id
    size_t myBoxSize[3] = { 1, 1, 1 };
    for ( int d      = 0; d < GeomDim; d++ )
        myBoxSize[d] = range[2 * d + 1] - range[2 * d + 0] + 1;
    size_t local_id = elem_id.local_id();
    auto type = elem_id.type();
    auto side = local_id / ( myBoxSize[0] * myBoxSize[1] * myBoxSize[2] );
    int x = range[0] + (int) local_id % myBoxSize[0];
    int y = range[2] + (int) ( local_id / myBoxSize[0] ) % myBoxSize[1];
    int z = range[4] + (int) ( local_id / ( myBoxSize[0] * myBoxSize[1] ) ) % myBoxSize[2];
    MeshElementIndex index( type, side, x, y, z );
    // Create the element
    structuredMeshElement elem( index, this );
    AMP_ASSERT( elem.globalID() == elem_id );
    return elem;
}
MeshElement BoxMesh::getElement( const MeshElementIndex &index ) const
{
    return structuredMeshElement( index, this );
}


/********************************************************
* Function to return parents of an element              *
********************************************************/
std::vector<MeshElement> BoxMesh::getElementParents( const MeshElement &meshelem,
                                                     const GeomType type ) const
{
    AMP_INSIST( meshelem.globalID().meshID() == d_meshID,
                "MeshElement is not from the given mesh" );
    AMP_INSIST( type >= meshelem.globalID().type() && type <= GeomDim,
                "Cannot get the parents of the given type for the current element" );
    if ( type == meshelem.globalID().type() )
        return std::vector<MeshElement>( 1, meshelem );
    // Get the element of interest
    const structuredMeshElement *elem =
        dynamic_cast<const structuredMeshElement *>( meshelem.getRawElement() );
    AMP_ASSERT( elem != nullptr );
    return elem->getParents( type );
}


/****************************************************************
* Functions to return the number of elements                    *
****************************************************************/
size_t BoxMesh::numLocalElements( const GeomType type ) const
{
    auto box = getLocalBlock( d_comm.getRank() );
    auto range = getIteratorRange( box, type, 0 );
    size_t N = 0;
    for ( const auto& tmp : range )
        N += BoxMesh::MeshElementIndex::numElements( tmp.first, tmp.second );
    return N;
}
size_t BoxMesh::numGlobalElements( const GeomType type ) const
{
    std::array<int,6> box = { 0, d_globalSize[0]-1, 0, d_globalSize[1]-1, 0, d_globalSize[2]-1 };
    auto range = getIteratorRange( box, type, 0 );
    size_t N = 0;
    for ( auto tmp : range )
        N += MeshElementIndex::numElements( tmp.first, tmp.second );
    return N;
}
size_t BoxMesh::numGhostElements( const GeomType type, int gcw ) const
{
    auto box = getLocalBlock( d_comm.getRank() );
    auto range1 = getIteratorRange( box, type, 0 );
    auto range2 = getIteratorRange( box, type, gcw );
    size_t N = 0;
    for (size_t i=0; i<range1.size(); i++) {
        size_t N1 = BoxMesh::MeshElementIndex::numElements( range1[i].first, range1[i].second );
        size_t N2 = BoxMesh::MeshElementIndex::numElements( range2[i].first, range2[i].second );
        N += N2-N1;
    }
    return N;
}


/****************************************************************
* Function to get an iterator                                   *
****************************************************************/
BoxMesh::ElementBlocks
BoxMesh::getIteratorRange( std::array<int,6> range, const GeomType type, const int gcw ) const
{
    AMP_ASSERT( type <= GeomDim );
    // Get the range of cells we care about
    for (int d=0; d<GeomDim; d++) {
        range[2*d+0] -= gcw;
        range[2*d+1] += gcw;
        if ( !d_isPeriodic[d] ) {
            range[2*d+0] = std::max(range[2*d+0],0);
            range[2*d+1] = std::min(range[2*d+1],d_globalSize[d]-1);
        }
    }
    // Get the element blocks we want to process
    ElementBlocks blocks;
    blocks.reserve(3);
    if ( type == GeomDim ) {
        blocks.push_back( std::make_pair(
            MeshElementIndex( type, 0, range[0], range[2], range[4] ),
            MeshElementIndex( type, 0, range[1], range[3], range[5] ) ) );
    } else if ( type == Vertex ) {
        for (int d=0; d<GeomDim; d++) {
            if ( gcw != 0 )
                range[2*d+1]++;
            else if ( !d_isPeriodic[d] && range[2*d+1]==d_globalSize[d]-1 )
                range[2*d+1]++;
        }
        blocks.push_back( std::make_pair(
            MeshElementIndex( type, 0, range[0], range[2], range[4] ),
            MeshElementIndex( type, 0, range[1], range[3], range[5] ) ) );
    } else if ( type==Edge && GeomDim==Face ) {
        auto range1 = range;
        auto range2 = range;
        if ( gcw!=0 || ( gcw==0 && !d_isPeriodic[0] && range[1]==d_globalSize[0]-1 ) )
            range2[1]++;
        if ( gcw!=0 || ( gcw==0 && !d_isPeriodic[1] && range[3]==d_globalSize[1]-1 ) )
            range1[3]++;
        blocks.push_back( std::make_pair(
            MeshElementIndex( type, 0, range1[0], range1[2], range1[4] ),
            MeshElementIndex( type, 0, range1[1], range1[3], range1[5] ) ) );
        blocks.push_back( std::make_pair(
            MeshElementIndex( type, 1, range2[0], range2[2], range2[4] ),
            MeshElementIndex( type, 1, range2[1], range2[3], range2[5] ) ) );
    } else if ( type==Edge && GeomDim==Volume ) {
        auto range1 = range;
        auto range2 = range;
        auto range3 = range;
        if ( gcw!=0 || ( gcw==0 && !d_isPeriodic[0] && range[1]==d_globalSize[0]-1 ) ) {
            range2[1]++;
            range3[1]++;
        }
        if ( gcw!=0 || ( gcw==0 && !d_isPeriodic[1] && range[3]==d_globalSize[1]-1 ) ) {
            range1[3]++;
            range3[3]++;
        }
        if ( gcw!=0 || ( gcw==0 && !d_isPeriodic[2] && range[5]==d_globalSize[2]-1 ) ) {
            range1[5]++;
            range2[5]++;
        }
        blocks.push_back( std::make_pair(
            MeshElementIndex( type, 0, range1[0], range1[2], range1[4] ),
            MeshElementIndex( type, 0, range1[1], range1[3], range1[5] ) ) );
        blocks.push_back( std::make_pair(
            MeshElementIndex( type, 1, range2[0], range2[2], range2[4] ),
            MeshElementIndex( type, 1, range2[1], range2[3], range2[5] ) ) );
        blocks.push_back( std::make_pair(
            MeshElementIndex( type, 2, range3[0], range3[2], range3[4] ),
            MeshElementIndex( type, 2, range3[1], range3[3], range3[5] ) ) );
    } else if ( type==Face && GeomDim==Volume ) {
        auto range1 = range;
        auto range2 = range;
        auto range3 = range;
        if ( gcw!=0 || ( gcw==0 && !d_isPeriodic[0] && range[1]==d_globalSize[0]-1 ) )
            range1[1]++;
        if ( gcw!=0 || ( gcw==0 && !d_isPeriodic[1] && range[3]==d_globalSize[1]-1 ) )
            range2[3]++;
        if ( gcw!=0 || ( gcw==0 && !d_isPeriodic[2] && range[5]==d_globalSize[2]-1 ) )
            range3[5]++;
        blocks.push_back( std::make_pair(
            MeshElementIndex( type, 0, range1[0], range1[2], range1[4] ),
            MeshElementIndex( type, 0, range1[1], range1[3], range1[5] ) ) );
        blocks.push_back( std::make_pair(
            MeshElementIndex( type, 1, range2[0], range2[2], range2[4] ),
            MeshElementIndex( type, 1, range2[1], range2[3], range2[5] ) ) );
        blocks.push_back( std::make_pair(
            MeshElementIndex( type, 2, range3[0], range3[2], range3[4] ),
            MeshElementIndex( type, 2, range3[1], range3[3], range3[5] ) ) );
    } else {
        AMP_ERROR("Unknown case");
    }
    // Check that each block does not have duplicate elements
    for ( auto& block : blocks ) {
        for (int d=0; d<GeomDim; d++) {
            if ( d_isPeriodic[d] ) {
                auto& first = block.first;
                auto& last = block.second;
                if ( first.index(d)+d_globalSize[d] <= last.index(d) ) {
                    first.index(d) = 0;
                    last.index(d) = d_globalSize[d]-1;
                }
            }
        }
    }
    return blocks;
}
BoxMesh::ElementBlocks BoxMesh::intersect( const ElementBlocks& set1, const ElementBlocks& set2 )
{
    ElementBlocks set;
    set.reserve(set1.size()*set2.size());
    for ( const auto v1 : set1 ) {
        for ( const auto v2 : set2 ) {
            if ( v1.first.type()!=v2.first.type() || v1.first.side() !=v2.first.side() )
                continue;
            auto v = v1;
            v.first.index(0)  = std::max( v1.first.index(0),  v2.first.index(0)  );
            v.first.index(1)  = std::max( v1.first.index(1),  v2.first.index(1)  );
            v.first.index(2)  = std::max( v1.first.index(2),  v2.first.index(2)  );
            v.second.index(0) = std::min( v1.second.index(0), v2.second.index(0) );
            v.second.index(1) = std::min( v1.second.index(1), v2.second.index(1) );
            v.second.index(2) = std::min( v1.second.index(2), v2.second.index(2) );
            if ( MeshElementIndex::numElements(v.first,v.second) > 0 )
                set.push_back(v);
        }
    }
    return set;
}
inline MeshIterator BoxMesh::createIterator( const ElementBlocks& list ) const
{
    if ( list.empty() ) {
        return MeshIterator();
    } else if ( list.size()==1 ) {
        return structuredMeshIterator( list[0].first, list[0].second, this, 0 );
    } else {
        std::vector<AMP::shared_ptr<MeshIterator>> iterator_list;
        iterator_list.reserve( list.size() );
        for ( const auto item : list ) {
            if ( MeshElementIndex::numElements(item.first,item.second) ) {
                AMP::shared_ptr<structuredMeshIterator> it( 
                    new structuredMeshIterator( item.first, item.second, this, 0 ) );
                iterator_list.push_back( it );
            }
        }
        return MultiIterator( iterator_list, 0 );
    }
    return MeshIterator();
}
MeshIterator BoxMesh::getIterator( const GeomType type, const int gcw ) const
{
    auto box = getLocalBlock( d_comm.getRank() );
    auto range = getIteratorRange( box, type, gcw );
    return createIterator( range );
}


/****************************************************************
* Function to get an iterator over the surface                  *
****************************************************************/
MeshIterator BoxMesh::getSurfaceIterator( const GeomType type, const int gcw ) const
{
    // Include each surface as needed
    ElementBlocks sufaceSet;
    for (int i=0; i<2*GeomDim; i++) {
        if ( d_onSurface[i] ) {
            for ( const auto& item : d_globalSurfaceList[i][type] )
                sufaceSet.emplace_back( item );
        }
    }
    // Intersect with the local ghost box
    auto box = getLocalBlock( d_comm.getRank() );
    auto range = getIteratorRange( box, type, gcw );
    auto intersection = intersect( sufaceSet, range );
    // Create a list if elements removing any duplicate elements
    std::set<MeshElementIndex> set;
    for ( auto block : intersection ) {
        auto type = block.first.type();
        auto side = block.first.side();
        for (int k=block.second.index(2); k>=block.first.index(2); k--) {
            for (int j=block.second.index(1); j>=block.first.index(1); j--) {
                for (int i=block.second.index(0); i>=block.first.index(0); i--) {
                    set.emplace( MeshElementIndex( type, side, i, j, k ) );
                }
            }
        }
    }
    // Create the iterator
    AMP::shared_ptr<std::vector<MeshElementIndex>> elements(
        new std::vector<MeshElementIndex>( set.begin(), set.end() ) );
    return structuredMeshIterator( elements, this, 0 );
}


/****************************************************************
* Functions to get the boundaries                               *
****************************************************************/
std::vector<int> BoxMesh::getBoundaryIDs() const
{
    std::set<int> ids;
    for (int i=0; i<2*GeomDim; i++) {
        if ( !d_isPeriodic[i/2] && d_surfaceId[i]!=-1 )
            ids.insert(d_surfaceId[i]);
    }
    return std::vector<int>(ids.begin(),ids.end());
}
MeshIterator
BoxMesh::getBoundaryIDIterator( const GeomType type, const int id, const int gcw ) const
{
    // Include each surface as needed
    ElementBlocks sufaceSet;
    for (int i=0; i<2*GeomDim; i++) {
        if ( d_surfaceId[i] == id ) {
            for ( auto item : d_globalSurfaceList[i][type] )
                sufaceSet.emplace_back( item );
        }
    }
    // Intersect with the local ghost box
    auto box = getLocalBlock( d_comm.getRank() );
    auto range = getIteratorRange( box, type, gcw );
    auto intersection = intersect( sufaceSet, range );
    // Create a list if elements removing any duplicate elements
    std::set<MeshElementIndex> set;
    for ( auto block : intersection ) {
        auto type = block.first.type();
        auto side = block.first.side();
        for (int k=block.second.index(2); k>=block.first.index(2); k--) {
            for (int j=block.second.index(1); j>=block.first.index(1); j--) {
                for (int i=block.second.index(0); i>=block.first.index(0); i--) {
                    set.emplace( MeshElementIndex( type, side, i, j, k ) );
                }
            }
        }
    }
    // Create the iterator
    AMP::shared_ptr<std::vector<MeshElementIndex>> elements(
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
    bool on_boundary = false;
    for (int i=0; i<6; i++) {
        if ( d_surfaceId[i]==id ) {
            int d = i/2;
            int s = i%2;
            if ( index.type() == GeomDim ) {
                on_boundary = on_boundary || ( s==0 && index.index(d)==0 );
                on_boundary = on_boundary || ( s==1 && index.index(d)==d_globalSize[d]-1 );
            } else if ( index.type() == Vertex ) {
                on_boundary = on_boundary || ( s==0 && index.index(d)==0 );
                on_boundary = on_boundary || ( s==1 && index.index(d)==d_globalSize[d] );
            } else if ( index.type() == Edge ) {
                on_boundary = on_boundary || ( s==0 && index.side()==d && index.index(d)==0 );
                on_boundary = on_boundary || ( s==1 && index.side()==d && index.index(d)==d_globalSize[d] );
            } else if ( index.type() == Face ) {
                on_boundary = on_boundary || ( s==0 && index.side()==d && index.index(d)==0 );
                on_boundary = on_boundary || ( s==1 && index.side()==d && index.index(d)==d_globalSize[d] );
            } else {
                AMP_ERROR("Unknown type");
            }
        }
    }
    return on_boundary;
}


/****************************************************************
* Helper function to return the indices and rank of the owning  *
* block for a given MeshElementIndex                            *
****************************************************************/
void BoxMesh::getOwnerBlock( const MeshElementIndex &index, unsigned int &rank, int *range ) const
{
    int myBoxIndex[3] = { 0, 0, 0 };
    for ( int d = 0; d < GeomDim; d++ ) {
        size_t size     = (size_t) d_globalSize[d];
        size_t N_blocks = (size_t) d_numBlocks[d];
        if ( index.index(d) == d_globalSize[d] ) {
            // The element lies on the physical bounadry
            AMP_ASSERT( index.type() < GeomDim );
            myBoxIndex[d] = d_numBlocks[d] - 1;
        } else {
            // Find the owning box
            myBoxIndex[d] = (int) ( ( ( (size_t) index.index(d) + 1 ) * N_blocks - 1 ) / size );
        }
        range[2 * d + 0] = (int) ( ( size * ( (size_t) myBoxIndex[d] ) ) / N_blocks );
        range[2 * d + 1] = (int) ( ( size * ( (size_t) myBoxIndex[d] + 1 ) ) / N_blocks );
        range[2 * d + 1] = std::min( range[2 * d + 1], d_globalSize[d] );
    }
    // Increase the index range for the boxes on the boundary for all elements except the current
    // dimension
    if ( index.type() != GeomDim ) {
        for ( int d = 0; d < GeomDim; d++ ) {
            if ( range[2 * d + 1] == d_globalSize[d] && !d_isPeriodic[d] )
                range[2 * d + 1]++;
        }
    }
    rank = (unsigned int) ( myBoxIndex[0] + myBoxIndex[1] * d_numBlocks[0] +
                            myBoxIndex[2] * d_numBlocks[0] * d_numBlocks[1] );
}


} // Mesh namespace
} // AMP namespace
