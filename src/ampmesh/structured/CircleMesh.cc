#include "AMP/ampmesh/structured/CircleMesh.h"
#include "AMP/ampmesh/shapes/Circle.h"
#include "AMP/ampmesh/structured/BoxMesh.h"


namespace AMP {
namespace Mesh {


/****************************************************************
 * Constructors                                                  *
 ****************************************************************/
CircleMesh::CircleMesh( MeshParameters::shared_ptr params ) : StructuredGeometryMesh( params )
{
    // Check for valid inputs
    AMP_INSIST( params.get(), "Params must not be null" );
    AMP_INSIST( d_comm != AMP_MPI( AMP_COMM_NULL ), "Communicator must be set" );
    AMP_INSIST( d_db.get(), "Database must exist" );
    // Input options from the database
    PhysicalDim           = d_db->getInteger( "dim" );
    GeomDim               = (GeomType) PhysicalDim;
    std::vector<int> size = d_db->getIntegerArray( "Size" );
    auto range            = d_db->getDoubleArray( "Range" );
    d_max_gcw             = d_db->getIntegerWithDefault( "GCW", 2 );
    AMP_INSIST( (int) PhysicalDim == 2, "dim must be size 2" );
    AMP_INSIST( size.size() == 1u, "Size must be an array of length 1" );
    AMP_INSIST( range.size() == 1u, "Size must be an array of length 1" );
    d_globalSize[0] = 2 * size[0];
    d_globalSize[1] = 2 * size[0];
    d_isPeriodic[0] = false;
    d_isPeriodic[1] = false;
    double R        = range[0];
    // Change the surface ids to match the standard ids
    d_surfaceId[0] = 1;
    d_surfaceId[1] = 1;
    d_surfaceId[2] = 1;
    d_surfaceId[3] = 1;
    d_onSurface[0] = true;
    d_onSurface[1] = true;
    d_onSurface[2] = true;
    d_onSurface[3] = true;
    // Initialize the logical mesh
    BoxMesh::initialize();
    // Set the geometry
    d_geometry.reset( new Geometry::Circle( R ) );
    // Finalize the logical mesh
    BoxMesh::finalize();
}


/****************************************************************
 * Estimate the mesh size                                        *
 ****************************************************************/
std::vector<size_t> CircleMesh::estimateLogicalMeshSize( const MeshParameters::shared_ptr &params )
{
    auto db   = params->getDatabase();
    auto size = db->getIntegerArray( "Size" );
    AMP_ASSERT( size.size() == 1u );
    std::vector<size_t> size2( 1, 2 * size[0] );
    return size2;
}


/****************************************************************
 * Copy the mesh                                                 *
 ****************************************************************/
AMP::shared_ptr<Mesh> CircleMesh::clone() const { return AMP::make_shared<CircleMesh>( *this ); }


} // namespace Mesh
} // namespace AMP
