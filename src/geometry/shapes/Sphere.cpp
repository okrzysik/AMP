#include "AMP/geometry/shapes/Sphere.h"
#include "AMP/IO/HDF5.h"
#include "AMP/geometry/GeometryHelpers.h"
#include "AMP/utils/Database.h"
#include "AMP/utils/UtilityMacros.h"


namespace AMP::Geometry {


/********************************************************
 * Constructor                                           *
 ********************************************************/
Sphere::Sphere( std::shared_ptr<const AMP::Database> db )
{
    d_ids         = { 4, 4, 4, 4, 2, 1 };
    d_isPeriodic  = { false, false, false };
    d_physicalDim = 3;
    d_logicalDim  = 3;
    d_offset[0]   = 0;
    d_offset[1]   = 0;
    d_offset[2]   = 0;
    auto range    = db->getVector<double>( "Range" );
    AMP_INSIST( range.size() == 1u, "Range must be an array of length 1" );
    d_r = range[0];
}
Sphere::Sphere( double r ) : LogicalGeometry(), d_r( r )
{
    d_ids         = { 4, 4, 4, 4, 2, 1 };
    d_isPeriodic  = { false, false, false };
    d_physicalDim = 3;
    d_logicalDim  = 3;
    d_offset[0]   = 0;
    d_offset[1]   = 0;
    d_offset[2]   = 0;
}


/********************************************************
 * Compute the nearest point on the surface              *
 ********************************************************/
Point Sphere::nearest( const Point &pos ) const
{
    // Get the current point in the reference frame of the circle
    double x = pos.x() - d_offset[0];
    double y = pos.y() - d_offset[1];
    double z = pos.z() - d_offset[2];
    // Calculate the nearest point
    double r = sqrt( x * x + y * y + z * z );
    if ( r <= d_r ) {
        return pos;
    } else {
        x *= d_r / r;
        y *= d_r / r;
        z *= d_r / r;
        return { x + d_offset[0], y + d_offset[1], z + d_offset[2] };
    }
}


/********************************************************
 * Compute the distance to the object                    *
 ********************************************************/
double Sphere::distance( const Point &pos, const Point &ang ) const
{
    double x = pos.x() - d_offset[0];
    double y = pos.y() - d_offset[1];
    double z = pos.z() - d_offset[2];
    double d = GeometryHelpers::distanceToSphere( d_r, { x, y, z }, ang );
    return d;
}


/********************************************************
 * Check if the ray is inside the geometry               *
 ********************************************************/
bool Sphere::inside( const Point &pos ) const
{
    double x = pos.x() - d_offset[0];
    double y = pos.y() - d_offset[1];
    double z = pos.z() - d_offset[2];
    return x * x + y * y + z * z <= ( 1.0 + 1e-12 ) * d_r * d_r;
}


/********************************************************
 * Return the closest surface                            *
 ********************************************************/
Point Sphere::surfaceNorm( const Point &pos ) const
{
    double x = pos.x() - d_offset[0];
    double y = pos.y() - d_offset[1];
    double z = pos.z() - d_offset[2];
    double r = sqrt( x * x + y * y + z * z );
    return { x / r, y / r, z / r };
}


/********************************************************
 * Return the physical coordinates                       *
 ********************************************************/
Point Sphere::physical( const Point &pos ) const
{
    auto point = GeometryHelpers::map_logical_sphere( d_r, pos.x(), pos.y(), pos.z() );
    point[0] += d_offset[0];
    point[1] += d_offset[1];
    point[2] += d_offset[2];
    return point;
}


/********************************************************
 * Return the logical coordinates                        *
 ********************************************************/
Point Sphere::logical( const Point &pos ) const
{
    double x = pos.x() - d_offset[0];
    double y = pos.y() - d_offset[1];
    double z = pos.z() - d_offset[2];
    return GeometryHelpers::map_sphere_logical( d_r, x, y, z );
}


/********************************************************
 * Return the centroid and bounding box                  *
 ********************************************************/
Point Sphere::centroid() const { return { d_offset[0], d_offset[1], d_offset[2] }; }
std::pair<Point, Point> Sphere::box() const
{
    Point lb = { d_offset[0] - d_r, d_offset[1] - d_r, d_offset[2] - d_r };
    Point ub = { d_offset[0] + d_r, d_offset[1] + d_r, d_offset[2] + d_r };
    return { lb, ub };
}


/********************************************************
 * Return the volume                                     *
 ********************************************************/
double Sphere::volume() const
{
    constexpr double pi = 3.141592653589793;
    return 4.0 / 3.0 * pi * d_r * d_r * d_r;
}


/********************************************************
 * Return the logical grid                               *
 ********************************************************/
std::vector<int> Sphere::getLogicalGridSize( const std::vector<int> &x ) const
{
    AMP_INSIST( x.size() == 1u, "Size must be an array of length 1" );
    return { 2 * x[0], 2 * x[0], 2 * x[0] };
}
std::vector<int> Sphere::getLogicalGridSize( const std::vector<double> &res ) const
{
    AMP_INSIST( res.size() == 3u, "Resolution must be an array of length 3" );
    double res2 = std::min( { res[0], res[1], res[2] } );
    int N       = std::max<int>( 2 * d_r / res2, 1 );
    return { N, N, N };
}


/********************************************************
 * Displace the mesh                                     *
 ********************************************************/
void Sphere::displace( const double *x )
{
    d_offset[0] += x[0];
    d_offset[1] += x[1];
    d_offset[2] += x[2];
}


/********************************************************
 * Clone the object                                      *
 ********************************************************/
std::unique_ptr<AMP::Geometry::Geometry> Sphere::clone() const
{
    return std::make_unique<Sphere>( *this );
}


/********************************************************
 * Compare the geometry                                  *
 ********************************************************/
bool Sphere::operator==( const Geometry &rhs ) const
{
    auto geom = dynamic_cast<const Sphere *>( &rhs );
    if ( !geom )
        return false;
    return d_r == geom->d_r && d_offset == geom->d_offset;
}


/****************************************************************
 * Write/Read restart data                                       *
 ****************************************************************/
void Sphere::writeRestart( int64_t fid ) const
{
    AMP::writeHDF5( fid, "GeomType", std::string( "sphere" ) );
    AMP::writeHDF5( fid, "physical", d_physicalDim ); // Geometry
    AMP::writeHDF5( fid, "logical", d_logicalDim );   // LogicalGeometry
    AMP::writeHDF5( fid, "periodic", d_isPeriodic );  // LogicalGeometry
    AMP::writeHDF5( fid, "ids", d_ids );              // LogicalGeometry
    AMP::writeHDF5( fid, "offset", d_offset );
    AMP::writeHDF5( fid, "R", d_r );
}
Sphere::Sphere( int64_t fid )
{
    AMP::readHDF5( fid, "physical", d_physicalDim ); // Geometry
    AMP::readHDF5( fid, "logical", d_logicalDim );   // LogicalGeometry
    AMP::readHDF5( fid, "periodic", d_isPeriodic );  // LogicalGeometry
    AMP::readHDF5( fid, "ids", d_ids );              // LogicalGeometry
    AMP::readHDF5( fid, "offset", d_offset );
    AMP::readHDF5( fid, "R", d_r );
}


} // namespace AMP::Geometry
