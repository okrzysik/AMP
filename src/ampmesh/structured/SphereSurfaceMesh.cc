#include "AMP/ampmesh/structured/SphereSurfaceMesh.h"
#include "AMP/ampmesh/shapes/SphereSurface.h"
#include "AMP/ampmesh/structured/BoxMesh.h"
#include "AMP/ampmesh/structured/BoxMeshHelpers.h"

namespace AMP {
namespace Mesh {


/****************************************************************
 * Constructors                                                  *
 ****************************************************************/
SphereSurfaceMesh::SphereSurfaceMesh( MeshParameters::shared_ptr params )
    : StructuredGeometryMesh( params )
{
    // Input options from the database
    PhysicalDim = 3;
    GeomDim     = GeomType::Face;
    auto size   = d_db->getIntegerArray( "Size" );
    auto range  = d_db->getDoubleArray( "Range" );
    d_max_gcw   = d_db->getIntegerWithDefault( "GCW", 2 );
    AMP_INSIST( size.size() == 1u, "Size must be an array of length 1" );
    AMP_INSIST( range.size() == 1u, "Range must be an array of length 1" );
    AMP_INSIST( (int) PhysicalDim == 3, "dim must be 3" );
    AMP_ASSERT( range[0] >= 0 );
    d_isPeriodic[0] = true;
    d_isPeriodic[1] = false;
    d_isPeriodic[2] = false;
    d_globalSize[0] = size[0];
    d_globalSize[1] = size[0] / 2;
    d_globalSize[2] = 1;
    // Change the surface ids to match the standard ids
    // 0,1,2,3 - 4: Outer surface
    // 4 - 2: Bottom surface
    // 5 - 1: Top surface
    d_surfaceId[0] = -1;
    d_surfaceId[1] = -1;
    d_surfaceId[2] = -1;
    d_surfaceId[3] = -1;
    d_surfaceId[4] = -1;
    d_surfaceId[5] = -1;
    d_onSurface[0] = false;
    d_onSurface[1] = false;
    d_onSurface[2] = false;
    d_onSurface[3] = false;
    d_onSurface[4] = false;
    d_onSurface[5] = false;
    // Initialize the logical mesh
    BoxMesh::initialize();
    // Set the geometry
    d_geometry.reset( new Geometry::SphereSurface( range[0] ) );
    // Finalize the logical mesh
    BoxMesh::finalize();
}


/****************************************************************
 * Estimate the mesh size                                        *
 ****************************************************************/
std::vector<size_t>
SphereSurfaceMesh::estimateLogicalMeshSize( const MeshParameters::shared_ptr &params )
{
    auto db               = params->getDatabase();
    std::vector<int> size = db->getIntegerArray( "Size" );
    AMP_ASSERT( size.size() == 1u );
    std::vector<size_t> size2( 2 );
    size2[0] = size[0];
    size2[1] = size[0] / 2;
    return size2;
}


/****************************************************************
 * Copy the mesh                                                 *
 ****************************************************************/
AMP::shared_ptr<Mesh> SphereSurfaceMesh::clone() const
{
    return AMP::make_shared<SphereSurfaceMesh>( *this );
}


} // namespace Mesh
} // namespace AMP
