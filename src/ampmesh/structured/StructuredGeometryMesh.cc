#include "AMP/ampmesh/structured/StructuredGeometryMesh.h"


namespace AMP {
namespace Mesh {


/****************************************************************
 * Constructors                                                 *
 ****************************************************************/
StructuredGeometryMesh::StructuredGeometryMesh( MeshParameters::shared_ptr params )
    : BoxMesh( params )
{
    // Basic defaults
    d_globalSize.fill( 1 );
    d_isPeriodic.fill( false );
    d_numBlocks.fill( 1 );
    // Check for valid inputs
    AMP_INSIST( params.get(), "Params must not be null" );
    AMP_INSIST( d_comm != AMP_MPI( AMP_COMM_NULL ), "Communicator must be set" );
    AMP_INSIST( d_db.get(), "Database must exist" );
    // Construct the geometry
    d_geometry = AMP::Geometry::Geometry::buildGeometry( d_db );
    // Fill basic mesh information
    PhysicalDim = d_geometry->getDim();
    GeomDim     = static_cast<AMP::Mesh::GeomType>( d_geometry->getLogicalDim() );
    d_max_gcw   = d_db->getIntegerWithDefault( "GCW", 2 );
    AMP_ASSERT( PhysicalDim == d_db->getIntegerWithDefault( "dim", PhysicalDim ) );
    auto size = d_geometry->getLogicalGridSize( d_db->getIntegerArray( "Size" ) );
    for ( size_t d = 0; d < size.size(); d++ )
        d_globalSize[d] = size[d];
    auto isPeriodic = d_geometry->getPeriodicDim();
    for ( size_t d = 0; d < isPeriodic.size(); d++ )
        d_isPeriodic[d] = isPeriodic[d];
    auto surfaceIds = d_geometry->getLogicalSurfaceIds();
    for ( size_t d = 0; d < surfaceIds.size(); d++ )
        d_surfaceId[d] = surfaceIds[d];
    // Initialize the logical mesh
    BoxMesh::initialize();
    BoxMesh::finalize();
}
StructuredGeometryMesh::StructuredGeometryMesh( const StructuredGeometryMesh &mesh )
    : BoxMesh( mesh )
{
    PhysicalDim  = mesh.PhysicalDim;
    GeomDim      = mesh.GeomDim;
    d_max_gcw    = mesh.d_max_gcw;
    d_comm       = mesh.d_comm;
    d_name       = mesh.d_name;
    d_box        = mesh.d_box;
    d_box_local  = mesh.d_box_local;
    d_isPeriodic = mesh.d_isPeriodic;
    d_globalSize = mesh.d_globalSize;
    d_blockSize  = mesh.d_blockSize;
    d_numBlocks  = mesh.d_numBlocks;
    d_surfaceId  = mesh.d_surfaceId;
    for ( int d = 0; d < 4; d++ ) {
        for ( int i = 0; i < 6; i++ )
            d_globalSurfaceList[i][d] = mesh.d_globalSurfaceList[i][d];
    }
}

/****************************************************************
 * Basic functions                                               *
 ****************************************************************/
Mesh::Movable StructuredGeometryMesh::isMeshMovable() const { return Mesh::Movable::Displace; }
void StructuredGeometryMesh::displaceMesh( const std::vector<double> &x )
{
    for ( int i = 0; i < PhysicalDim; i++ ) {
        d_box[2 * i + 0] += x[i];
        d_box[2 * i + 1] += x[i];
        d_box_local[2 * i + 0] += x[i];
        d_box_local[2 * i + 1] += x[i];
    }
    d_geometry->displaceMesh( x.data() );
}
void StructuredGeometryMesh::displaceMesh( AMP::shared_ptr<const AMP::LinearAlgebra::Vector> )
{
    AMP_ERROR( "displaceMesh (vector) violates StructuredGeometryMesh properties" );
}
AMP::Geometry::Point
StructuredGeometryMesh::physicalToLogical( const AMP::Geometry::Point &x ) const
{
    return d_geometry->logical( x );
}
void StructuredGeometryMesh::coord( const MeshElementIndex &index, double *pos ) const
{
    AMP_ASSERT( index.type() == AMP::Mesh::GeomType::Vertex );
    double x = static_cast<double>( index.index( 0 ) ) / static_cast<double>( d_globalSize[0] );
    double y = static_cast<double>( index.index( 1 ) ) / static_cast<double>( d_globalSize[1] );
    double z = static_cast<double>( index.index( 2 ) ) / static_cast<double>( d_globalSize[2] );
    auto tmp = d_geometry->physical( AMP::Geometry::Point( x, y, z ) );
    for ( int d = 0; d < PhysicalDim; d++ )
        pos[d] = tmp[d];
}
AMP::shared_ptr<Mesh> StructuredGeometryMesh::clone() const
{
    return AMP::make_shared<StructuredGeometryMesh>( *this );
}


} // namespace Mesh
} // namespace AMP
