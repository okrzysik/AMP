#include "AMP/ampmesh/structured/CircleFrustumMesh.h"
#include "AMP/ampmesh/shapes/CircleFrustum.h"
#include "AMP/ampmesh/structured/BoxMesh.h"


namespace AMP {
namespace Mesh {


/****************************************************************
 * Constructors                                                  *
 ****************************************************************/
CircleFrustumMesh::CircleFrustumMesh( MeshParameters::shared_ptr params )
    : StructuredGeometryMesh( params )
{
    PhysicalDim = 3;
    GeomDim     = GeomType::Volume;
    // Input options from the database
    auto size = d_db->getIntegerArray( "Size" );
    double r1 = d_db->getDouble( "BaseRadius" );
    double r2 = d_db->getDouble( "TopRadius" );
    double h  = d_db->getDouble( "Height" );
    auto dir  = d_db->getString( "Dir" );
    d_max_gcw = d_db->getIntegerWithDefault( "GCW", 2 );
    AMP_INSIST( size.size() == 2u, "Size must be an array of length 2" );
    d_globalSize[0] = 2 * size[0];
    d_globalSize[1] = 2 * size[0];
    d_globalSize[2] = size[1];
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
    d_geometry.reset( new Geometry::CircleFrustum( { r1, r2 }, dir2, h ) );
    // Finalize the logical mesh
    BoxMesh::finalize();
}


/****************************************************************
 * Estimate the mesh size                                        *
 ****************************************************************/
std::vector<size_t>
CircleFrustumMesh::estimateLogicalMeshSize( const MeshParameters::shared_ptr &params )
{
    auto db   = params->getDatabase();
    auto size = db->getIntegerArray( "Size" );
    AMP_INSIST( size.size() == 2, "Size must be an array of length 2" );
    return { (size_t) size[0], (size_t) size[0], (size_t) size[1] };
}


/****************************************************************
 * Copy the mesh                                                 *
 ****************************************************************/
AMP::shared_ptr<Mesh> CircleFrustumMesh::clone() const
{
    return AMP::make_shared<CircleFrustumMesh>( *this );
}


} // namespace Mesh
} // namespace AMP
