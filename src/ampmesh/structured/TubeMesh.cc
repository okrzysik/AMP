#include "AMP/ampmesh/structured/TubeMesh.h"
#include "AMP/ampmesh/MultiIterator.h"
#include "AMP/ampmesh/shapes/Tube.h"

namespace AMP {
namespace Mesh {


/****************************************************************
 * Constructors                                                  *
 ****************************************************************/
TubeMesh::TubeMesh( MeshParameters::shared_ptr params ) : StructuredGeometryMesh( params )
{
    for ( int d = 0; d < 3; d++ ) {
        d_globalSize[d] = 1;
        d_isPeriodic[d] = false;
        d_numBlocks[d]  = 0;
    }
    // Check for valid inputs
    AMP_INSIST( params.get(), "Params must not be null" );
    AMP_INSIST( d_comm != AMP_MPI( AMP_COMM_NULL ), "Communicator must be set" );
    AMP_INSIST( d_db.get(), "Database must exist" );
    // Input options from the database
    PhysicalDim = d_db->getInteger( "dim" );
    GeomDim     = (GeomType) PhysicalDim;
    auto size   = d_db->getIntegerArray( "Size" );
    auto range  = d_db->getDoubleArray( "Range" );
    d_max_gcw   = d_db->getIntegerWithDefault( "GCW", 2 );
    std::vector<unsigned char> per( 1, false );
    if ( d_db->keyExists( "Periodic" ) )
        per = d_db->getBoolArray( "Periodic" );
    AMP_INSIST( size.size() == 3u, "Size must be an array of length 3" );
    AMP_INSIST( per.size() == 1u, "Periodic must be an array of length 1" );
    AMP_INSIST( range.size() == 4u, "Range must be an array of length 4" );
    AMP_INSIST( PhysicalDim == 3, "dim must be 3" );
    for ( int i = 0; i < 3; i++ )
        d_globalSize[i] = size[i];
    // Change the surface ids to match the standard ids
    // 0 - 8: Inner surface
    // 1 - 4: Outer surface
    // 4 - 2: Bottom surface
    // 5 - 1: Top surface
    d_surfaceId[0] = 8;
    d_surfaceId[1] = 4;
    d_surfaceId[2] = -1;
    d_surfaceId[3] = -1;
    d_surfaceId[4] = 2;
    d_surfaceId[5] = 1;
    d_onSurface[0] = true;
    d_onSurface[1] = true;
    d_onSurface[2] = false;
    d_onSurface[3] = false;
    d_onSurface[4] = true;
    d_onSurface[5] = true;
    // Initialize the logical mesh
    BoxMesh::initialize();
    // Set the geometry
    d_geometry.reset( new Geometry::Tube( range[0], range[1], range[2], range[3] ) );
    // Finalize the logical mesh
    BoxMesh::finalize();
}


/****************************************************************
 * Estimate the mesh size                                        *
 ****************************************************************/
std::vector<size_t> TubeMesh::estimateLogicalMeshSize( const MeshParameters::shared_ptr &params )
{
    auto db               = params->getDatabase();
    std::vector<int> size = db->getIntegerArray( "Size" );
    AMP_ASSERT( size.size() == 3u );
    std::vector<size_t> size2( size.size() );
    for ( size_t d = 0; d < size.size(); d++ )
        size2[d] = size[d];
    return size2;
}


/****************************************************************
 * Copy the mesh                                                 *
 ****************************************************************/
AMP::shared_ptr<Mesh> TubeMesh::clone() const { return AMP::make_shared<TubeMesh>( *this ); }


} // namespace Mesh
} // namespace AMP
