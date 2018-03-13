#include "AMP/ampmesh/structured/CylinderMesh.h"
#include "AMP/ampmesh/shapes/Cylinder.h"
#include "AMP/ampmesh/structured/BoxMesh.h"


namespace AMP {
namespace Mesh {


/****************************************************************
 * Constructors                                                  *
 ****************************************************************/
CylinderMesh::CylinderMesh( MeshParameters::shared_ptr params ) : StructuredGeometryMesh( params )
{
    // Input options from the database
    PhysicalDim = d_db->getInteger( "dim" );
    GeomDim     = (GeomType) PhysicalDim;
    auto size   = d_db->getIntegerArray( "Size" );
    auto range  = d_db->getDoubleArray( "Range" );
    d_max_gcw   = d_db->getIntegerWithDefault( "GCW", 2 );
    std::vector<unsigned char> per( 1, false );
    if ( d_db->keyExists( "Periodic" ) )
        per = d_db->getBoolArray( "Periodic" );
    AMP_INSIST( per.size() == 1u, "Periodic must be an array of length 1" );
    AMP_INSIST( size.size() == 2u, "Size must be an array of length 2" );
    AMP_INSIST( range.size() == 3u, "Range must be an array of length 3" );
    AMP_INSIST( (int) PhysicalDim == 3, "dim must be 3" );
    d_globalSize[0] = 2 * size[0];
    d_globalSize[1] = 2 * size[0];
    d_globalSize[2] = size[1];
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
    d_geometry.reset( new Geometry::Cylinder( range[0], range[1], range[2] ) );
    // Finalize the logical mesh
    BoxMesh::finalize();
}


/****************************************************************
 * Estimate the mesh size                                        *
 ****************************************************************/
std::vector<size_t>
CylinderMesh::estimateLogicalMeshSize( const MeshParameters::shared_ptr &params )
{
    auto db               = params->getDatabase();
    std::vector<int> size = db->getIntegerArray( "Size" );
    AMP_ASSERT( size.size() == 2u );
    std::vector<size_t> size2( 3 );
    size2[0] = 2 * size[0];
    size2[1] = 2 * size[0];
    size2[2] = size[1];
    return size2;
}


/****************************************************************
 * Copy the mesh                                                 *
 ****************************************************************/
AMP::shared_ptr<Mesh> CylinderMesh::clone() const
{
    return AMP::make_shared<CylinderMesh>( *this );
}


} // namespace Mesh
} // namespace AMP
