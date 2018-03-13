#include "AMP/ampmesh/structured/SphereMesh.h"
#include "AMP/ampmesh/shapes/Sphere.h"
#include "AMP/ampmesh/structured/BoxMesh.h"
#include "AMP/ampmesh/structured/BoxMeshHelpers.h"


namespace AMP {
namespace Mesh {


/****************************************************************
 * Constructors                                                  *
 ****************************************************************/
SphereMesh::SphereMesh( MeshParameters::shared_ptr params ) : StructuredGeometryMesh( params )
{
    // Input options from the database
    PhysicalDim = d_db->getInteger( "dim" );
    GeomDim     = (GeomType) PhysicalDim;
    auto size   = d_db->getIntegerArray( "Size" );
    auto range  = d_db->getDoubleArray( "Range" );
    d_max_gcw   = d_db->getIntegerWithDefault( "GCW", 2 );
    AMP_INSIST( size.size() == 1u, "Size must be an array of length 1" );
    AMP_INSIST( range.size() == 1u, "Range must be an array of length 3" );
    AMP_INSIST( (int) PhysicalDim == 3, "dim must be 3" );
    d_globalSize[0] = 2 * size[0];
    d_globalSize[1] = 2 * size[0];
    d_globalSize[2] = 2 * size[0];
    // Change the surface ids to match the standard ids
    // 0,1,2,3 - 4: Outer surface
    // 4 - 2: Bottom surface
    // 5 - 1: Top surface
    d_surfaceId[0] = 4;
    d_surfaceId[1] = 4;
    d_surfaceId[2] = 4;
    d_surfaceId[3] = 4;
    d_surfaceId[4] = 2;
    d_surfaceId[5] = 1;
    d_onSurface[0] = true;
    d_onSurface[1] = true;
    d_onSurface[2] = true;
    d_onSurface[3] = true;
    d_onSurface[4] = true;
    d_onSurface[5] = true;
    // Initialize the logical mesh
    BoxMesh::initialize();
    // Set the geometry
    d_geometry.reset( new Geometry::Sphere( range[0] ) );
    // Finalize the logical mesh
    BoxMesh::finalize();
}


/****************************************************************
 * Estimate the mesh size                                        *
 ****************************************************************/
std::vector<size_t> SphereMesh::estimateLogicalMeshSize( const MeshParameters::shared_ptr &params )
{
    auto db               = params->getDatabase();
    std::vector<int> size = db->getIntegerArray( "Size" );
    AMP_ASSERT( size.size() == 1u );
    return std::vector<size_t>( 3, 2 * size[0] );
}


/****************************************************************
 * Copy the mesh                                                 *
 ****************************************************************/
AMP::shared_ptr<Mesh> SphereMesh::clone() const { return AMP::make_shared<SphereMesh>( *this ); }


} // namespace Mesh
} // namespace AMP
