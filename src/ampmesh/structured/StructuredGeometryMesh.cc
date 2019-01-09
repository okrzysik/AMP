#include "AMP/ampmesh/structured/StructuredGeometryMesh.h"


namespace AMP {
namespace Mesh {


StructuredGeometryMesh::StructuredGeometryMesh( MeshParameters::shared_ptr params )
    : BoxMesh( params )
{
}
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


} // namespace Mesh
} // namespace AMP
