#include "AMP/ampmesh/structured/ShellMesh.h"
#include "AMP/ampmesh/shapes/Shell.h"
#include "AMP/ampmesh/structured/BoxMesh.h"

namespace AMP {
namespace Mesh {


/****************************************************************
 * Constructors                                                  *
 ****************************************************************/
ShellMesh::ShellMesh( MeshParameters::shared_ptr params ) : StructuredGeometryMesh( params )
{
    // Input options from the database
    PhysicalDim = d_db->getInteger( "dim" );
    GeomDim     = AMP::Mesh::GeomType::Face;
    auto size   = d_db->getIntegerArray( "Size" );
    auto range  = d_db->getDoubleArray( "Range" );
    d_max_gcw   = d_db->getIntegerWithDefault( "GCW", 2 );
    AMP_INSIST( size.size() == 2u, "Size must be an array of length 2" );
    AMP_INSIST( range.size() == 2u, "Range must be an array of length 2" );
    AMP_INSIST( (int) PhysicalDim == 3, "dim must be 3" );
    AMP_ASSERT( range[0] >= 0 && range[1] > 0 && ( range[1] - range[0] ) > 0 );
    d_isPeriodic[0] = true;
    d_isPeriodic[1] = false;
    d_isPeriodic[2] = false;
    d_globalSize[0] = size[1];
    d_globalSize[1] = size[1] / 2;
    d_globalSize[2] = size[0];
    // Change the surface ids to match the standard ids
    // 0,1,2,3 - 4: Outer surface
    // 4 - 2: Bottom surface
    // 5 - 1: Top surface
    d_surfaceId[0] = -1;
    d_surfaceId[1] = -1;
    d_surfaceId[2] = 1;
    d_surfaceId[3] = 2;
    d_surfaceId[4] = 3;
    d_surfaceId[5] = 4;
    d_onSurface[0] = false;
    d_onSurface[1] = false;
    d_onSurface[2] = true;
    d_onSurface[3] = true;
    d_onSurface[4] = true;
    d_onSurface[5] = true;
    // Initialize the logical mesh
    BoxMesh::initialize();
    // Set the geometry
    d_geometry.reset( new Geometry::Shell( range[0], range[1] ) );
    // Finalize the logical mesh
    BoxMesh::finalize();
}


/****************************************************************
 * Estimate the mesh size                                        *
 ****************************************************************/
std::vector<size_t> ShellMesh::estimateLogicalMeshSize( const MeshParameters::shared_ptr &params )
{
    auto db   = params->getDatabase();
    auto size = db->getIntegerArray( "Size" );
    AMP_ASSERT( size.size() == 2u );
    std::vector<size_t> size2( 3 );
    size2[0] = size[1];
    size2[1] = size[1] / 2;
    size2[2] = size[0];
    return size2;
}


/****************************************************************
 * Copy the mesh                                                 *
 ****************************************************************/
AMP::shared_ptr<Mesh> ShellMesh::clone() const { return AMP::make_shared<ShellMesh>( *this ); }


} // namespace Mesh
} // namespace AMP
