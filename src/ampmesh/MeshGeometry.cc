#include "AMP/ampmesh/MeshGeometry.h"
#include "AMP/ampmesh/Mesh.h"


namespace AMP::Geometry {


/********************************************************
 * Constructors                                          *
 ********************************************************/
MeshGeometry::MeshGeometry( std::unique_ptr<AMP::Mesh::Mesh> mesh ) : d_mesh( std::move( mesh ) )
{
    AMP_ASSERT( d_mesh );
    AMP_ASSERT( static_cast<int>( d_mesh->getGeomType() ) == getDim() - 1 );
    d_physicalDim = d_mesh->getDim();
}
std::unique_ptr<AMP::Geometry::Geometry> MeshGeometry::clone() const
{
    return std::unique_ptr<AMP::Geometry::Geometry>( new MeshGeometry( d_mesh->clone() ) );
}


/********************************************************
 * Get the distance to the surface                       *
 ********************************************************/
double MeshGeometry::distance( const Point &pos, const Point &dir ) const
{
    NULL_USE( pos );
    NULL_USE( dir );
    AMP_ERROR( "Not finished" );
    return 0;
}
bool MeshGeometry::inside( const Point &pos ) const
{
    double dist = distance( pos, Point( pos.ndim(), { 0, 0, 0 } ) );
    return dist <= 0;
}


/********************************************************
 * Get the surface                                       *
 ********************************************************/
int MeshGeometry::NSurface() const
{
    AMP_ERROR( "Not finished" );
    return 0;
}
int MeshGeometry::surface( const Point &x ) const
{
    NULL_USE( x );
    AMP_ERROR( "Not finished" );
    return 0;
}
Point MeshGeometry::surfaceNorm( const Point &x ) const
{
    NULL_USE( x );
    AMP_ERROR( "Not finished" );
    return { 0, 0, 0 };
}


/********************************************************
 * Get the centroid/box                                  *
 ********************************************************/
Point MeshGeometry::centroid() const
{
    AMP_ERROR( "Not finished" );
    return { 0, 0, 0 };
}
std::pair<Point, Point> MeshGeometry::box() const
{
    AMP_ERROR( "Not finished" );
    Point p0 = { 0, 0, 0 };
    return std::make_pair( p0, p0 );
}


/********************************************************
 * Get the centroid/box                                  *
 ********************************************************/
void MeshGeometry::displace( const double *x )
{
    std::vector<double> x2( d_mesh->getDim() );
    for ( size_t i = 0; i < x2.size(); i++ )
        x2[i] = x[i];
    d_mesh->displaceMesh( x2 );
}


} // namespace AMP::Geometry
