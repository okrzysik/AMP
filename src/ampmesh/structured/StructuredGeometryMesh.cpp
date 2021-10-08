#include "AMP/ampmesh/structured/StructuredGeometryMesh.h"
#include "AMP/ampmesh/MeshParameters.h"


namespace AMP {
namespace Mesh {


/****************************************************************
 * Constructors                                                 *
 ****************************************************************/
StructuredGeometryMesh::StructuredGeometryMesh( std::shared_ptr<const MeshParameters> params )
    : BoxMesh( params ), d_pos_hash( 0 )
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
    auto db2 = d_db->cloneDatabase();
    db2->erase( "x_offset", false );
    db2->erase( "y_offset", false );
    db2->erase( "z_offset", false );
    d_geometry  = AMP::Geometry::Geometry::buildGeometry( std::move( db2 ) );
    d_geometry2 = std::dynamic_pointer_cast<AMP::Geometry::LogicalGeometry>( d_geometry );
    AMP_ASSERT( d_geometry2 );
    // Fill basic mesh information
    PhysicalDim = d_geometry2->getDim();
    GeomDim     = static_cast<AMP::Mesh::GeomType>( d_geometry2->getLogicalDim() );
    d_max_gcw   = d_db->getWithDefault<int>( "GCW", 2 );
    AMP_ASSERT( PhysicalDim == d_db->getWithDefault<int>( "dim", PhysicalDim ) );
    auto size = d_geometry2->getLogicalGridSize( d_db->getVector<int>( "Size" ) );
    AMP_ASSERT( size.size() == static_cast<size_t>( GeomDim ) );
    for ( size_t d = 0; d < size.size(); d++ )
        d_globalSize[d] = size[d];
    auto isPeriodic = d_geometry2->getPeriodicDim();
    for ( size_t d = 0; d < isPeriodic.size(); d++ )
        d_isPeriodic[d] = isPeriodic[d];
    auto surfaceIds = d_geometry2->getLogicalSurfaceIds();
    for ( size_t d = 0; d < surfaceIds.size(); d++ )
        d_surfaceId[d] = surfaceIds[d];
    // Initialize the logical mesh
    BoxMesh::initialize();
    BoxMesh::finalize();
}
StructuredGeometryMesh::StructuredGeometryMesh( const StructuredGeometryMesh &mesh )
    : BoxMesh( mesh )
{
    d_geometry2 = std::dynamic_pointer_cast<AMP::Geometry::LogicalGeometry>( d_geometry );
    d_pos_hash  = mesh.d_pos_hash;
}


/********************************************************
 * Return the class name                                 *
 ********************************************************/
std::string StructuredGeometryMesh::meshClass() const { return "StructuredGeometryMesh"; }


/****************************************************************
 * Basic functions                                               *
 ****************************************************************/
Mesh::Movable StructuredGeometryMesh::isMeshMovable() const { return Mesh::Movable::Displace; }
uint64_t StructuredGeometryMesh::positionHash() const { return d_pos_hash; }
void StructuredGeometryMesh::displaceMesh( const std::vector<double> &x )
{
    for ( int i = 0; i < PhysicalDim; i++ ) {
        d_box[2 * i + 0] += x[i];
        d_box[2 * i + 1] += x[i];
        d_box_local[2 * i + 0] += x[i];
        d_box_local[2 * i + 1] += x[i];
    }
    d_geometry2->displace( x.data() );
    d_pos_hash++;
}
void StructuredGeometryMesh::displaceMesh( std::shared_ptr<const AMP::LinearAlgebra::Vector> )
{
    AMP_ERROR( "displaceMesh (vector) violates StructuredGeometryMesh properties" );
}
AMP::Geometry::Point
StructuredGeometryMesh::physicalToLogical( const AMP::Geometry::Point &x ) const
{
    return d_geometry2->logical( x );
}
void StructuredGeometryMesh::coord( const MeshElementIndex &index, double *pos ) const
{
    AMP_ASSERT( index.type() == AMP::Mesh::GeomType::Vertex );
    double x = static_cast<double>( index.index( 0 ) ) / static_cast<double>( d_globalSize[0] );
    double y = static_cast<double>( index.index( 1 ) ) / static_cast<double>( d_globalSize[1] );
    double z = static_cast<double>( index.index( 2 ) ) / static_cast<double>( d_globalSize[2] );
    auto tmp = d_geometry2->physical( AMP::Geometry::Point( x, y, z ) );
    for ( int d = 0; d < PhysicalDim; d++ )
        pos[d] = tmp[d];
}
std::unique_ptr<Mesh> StructuredGeometryMesh::clone() const
{
    return std::make_unique<StructuredGeometryMesh>( *this );
}


/****************************************************************
 * Check if two meshes are equal                                 *
 ****************************************************************/
bool StructuredGeometryMesh::operator==( const Mesh &rhs ) const
{
    // Check base class variables
    if ( !BoxMesh::operator==( rhs ) )
        return false;
    // Check if we can cast to a MovableBoxMesh
    auto mesh = dynamic_cast<const StructuredGeometryMesh *>( &rhs );
    if ( !mesh )
        return false;
    // Perform basic comparison
    return d_geometry2 == mesh->d_geometry2;
}


} // namespace Mesh
} // namespace AMP
