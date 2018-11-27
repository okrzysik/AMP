#include "AMP/ampmesh/structured/SquareFrustumMesh.h"
#include "AMP/ampmesh/shapes/SquareFrustum.h"
#include "AMP/ampmesh/structured/BoxMesh.h"
#include "AMP/ampmesh/structured/BoxMeshHelpers.h"


namespace AMP {
namespace Mesh {


/****************************************************************
 * Constructors                                                  *
 ****************************************************************/
SquareFrustumMesh::SquareFrustumMesh( MeshParameters::shared_ptr params )
    : StructuredGeometryMesh( params )
{
    // Input options from the database
    PhysicalDim   = 3;
    GeomDim       = GeomType::Volume;
    auto size     = d_db->getIntegerArray( "Size" );
    auto range    = d_db->getDoubleArray( "Range" );
    auto dir      = d_db->getString( "Dir" );
    double height = d_db->getDouble( "Height" );
    d_max_gcw     = d_db->getIntegerWithDefault( "GCW", 2 );
    AMP_INSIST( size.size() == 3u, "Size must be an array of length 3" );
    AMP_INSIST( range.size() == 6u, "Range must be an array of length 6" );
    d_globalSize[0] = size[0];
    d_globalSize[1] = size[1];
    d_globalSize[2] = size[2];
    // Change the surface ids to match the standard ids
    for ( int i = 0; i < 6; i++ ) {
        d_surfaceId[i] = i;
        d_onSurface[i] = true;
    }
    // Initialize the logical mesh
    BoxMesh::initialize();
    // Set the geometry
    int dir2 = 0;
    if ( dir == "-x" )
        dir2 = 0;
    else if ( dir == "+x" )
        dir2 = 1;
    else if ( dir == "-y" )
        dir2 = 2;
    else if ( dir == "+y" )
        dir2 = 3;
    else if ( dir == "-z" )
        dir2 = 4;
    else if ( dir == "+z" )
        dir2 = 5;
    else
        AMP_ERROR( "Invalid value for Dir" );
    d_geometry.reset( new Geometry::SquareFrustum( range, dir2, height ) );
    // Finalize the logical mesh
    BoxMesh::finalize();
}


/****************************************************************
 * Estimate the mesh size                                        *
 ****************************************************************/
std::vector<size_t>
SquareFrustumMesh::estimateLogicalMeshSize( const MeshParameters::shared_ptr &params )
{
    auto db   = params->getDatabase();
    auto size = db->getIntegerArray( "Size" );
    AMP_INSIST( size.size() == 3u, "Size must be an array of length 3" );
    return { (size_t) size[0], (size_t) size[1], (size_t) size[2] };
}


/****************************************************************
 * Copy the mesh                                                 *
 ****************************************************************/
AMP::shared_ptr<Mesh> SquareFrustumMesh::clone() const
{
    return AMP::make_shared<SquareFrustumMesh>( *this );
}


} // namespace Mesh
} // namespace AMP
