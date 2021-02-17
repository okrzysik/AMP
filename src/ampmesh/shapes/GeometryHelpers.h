#ifndef included_AMP_GeometryHelpers
#define included_AMP_GeometryHelpers


#include <array>
#include <set>
#include <vector>


namespace AMP::Mesh {
template<class TYPE>
class MeshPoint;
using Point = class AMP::Mesh::MeshPoint<double>;
} // namespace AMP::Mesh

namespace AMP::Geometry::GeometryHelpers {


using Point2D = std::array<double, 2>;
using Point3D = std::array<double, 3>;


/**
 * \brief   Map logical coordinates to a circle
 * \details  This function will map logical coordinates in [0,1] to (x,y) coordinates
 *   in a circle.  It uses the mapping by:
 *   Dona Calhoun, Christiane Helzel, Randall LeVeque, "Logically Rectangular Grids
 *      and Finite GeomType::Volume Methods for PDEs in Circular and Spherical Domains",
 *      SIAM REVIEW, Vol. 50, No. 4, pp. 723–752 (2008)
 *   There are 3 methods to choose from:
 *      1 - D(d) = r*d/sqrt(2), R(d) = r*d
 *      2 - D(d) = r*d/sqrt(2), R(d) = r
 *      3 - D(d) = r*d*(2-d)/sqrt(2), R(d) = r1
 * \param[in] R         Radius of circle
 * \param[in] method    Method to map
 * \param[in] x         Logical x coordinate
 * \param[in] y         Logical y coordinate
 * @return              Returns a pair with the (x,y) value
 */
std::pair<double, double> map_logical_circle( double R, int method, double x, double y );


/**
 * \brief   Map physical coordinates to the logical coordinates
 * \details  This function will map physical coordinates (x,y) coordinates
 *   in a circle to [0,1] logical coordinates.  It uses the inverse of the mapping by:
 *   Dona Calhoun, Christiane Helzel, Randall LeVeque, "Logically Rectangular Grids
 *      and Finite GeomType::Volume Methods for PDEs in Circular and Spherical Domains",
 *      SIAM REVIEW, Vol. 50, No. 4, pp. 723–752 (2008)
 *   There are 3 methods to choose from:
 *      1 - D(d) = r*d/sqrt(2), R(d) = r*d
 *      2 - D(d) = r*d/sqrt(2), R(d) = r
 *      3 - D(d) = r*d*(2-d)/sqrt(2), R(d) = r1
 * \param[in] R         Radius of circle
 * \param[in] method    Method to map
 * \param[in] x         Physical x coordinate
 * \param[in] y         Physical y coordinate
 * @return              Returns a pair with the logical (x,y) value
 */
std::pair<double, double> map_circle_logical( double R, int method, double x, double y );


/**
 * \brief   Map logical coordinates to a regular polygon
 * \details  This function will map logical coordinates in [0,1] to (x,y) coordinates
 *   in a regular polygon.
 * \param[in] N         Number of faces
 * \param[in] R         Radius of circumcircle
 * \param[in] x         Logical x coordinate
 * \param[in] y         Logical y coordinate
 * @return              Returns a pair with the (x,y) value
 */
std::pair<double, double> map_logical_poly( int N, double R, double x, double y );


/**
 * \brief   Map physical coordinates to the logical coordinates
 * \details  This function will map physical coordinates (x,y) coordinates
 *   in a regular polygon to [0,1] logical coordinates.
 * \param[in] N         Number of faces
 * \param[in] R         Radius of circumcircle
 * \param[in] x         Physical x coordinate
 * \param[in] y         Physical y coordinate
 * @return              Returns a pair with the logical (x,y) value
 */
std::pair<double, double> map_poly_logical( int N, double R, double x, double y );


/**
 * \brief   Get the verticies of a regular polygon
 * \details  This function will return the verticies of a N-sided polygon
 *    that is compatible with the functions map_logical_poly and map_poly_logical.
 * \param[in] N         Number of faces
 * \param[in] R         Radius of circumcircle
 * @return              Returns the verticies clockwise (physical coordinates)
 */
std::vector<Point2D> get_poly_verticies( int N, double R );


/**
 * \brief   Map logical coordinates to a sphere
 * \details  This function will map logical coordinates in [0,1] to (x,y,z) coordinates
 *   in a sphere.  It uses the mapping by:
 *   Dona Calhoun, Christiane Helzel, Randall LeVeque, "Logically Rectangular Grids
 *      and Finite GeomType::Volume Methods for PDEs in Circular and Spherical Domains",
 *      SIAM REVIEW, Vol. 50, No. 4, pp. 723–752 (2008)
 * \param[in] R         Radius of sphere
 * \param[in] x         Logical x coordinate
 * \param[in] y         Logical y coordinate
 * \param[in] z         Logical z coordinate
 * @return              Returns the physical point
 */
Point3D map_logical_sphere( double R, double x, double y, double z );


/**
 * \brief   Map logical coordinates to a sphere
 * \details  This function will map physical coordinates in (x,y,z) to [0,1] locial
 *   coordinates in a sphere.  It uses the mapping by:
 *   Dona Calhoun, Christiane Helzel, Randall LeVeque, "Logically Rectangular Grids
 *      and Finite GeomType::Volume Methods for PDEs in Circular and Spherical Domains",
 *      SIAM REVIEW, Vol. 50, No. 4, pp. 723–752 (2008)
 * \param[in] R         Radius of sphere
 * \param[in] x         Physical x coordinate
 * \param[in] y         Physical y coordinate
 * \param[in] z         Physical z coordinate
 * @return              Returns the logical point
 */
Point3D map_sphere_logical( double R, double x, double y, double z );


/**
 * \brief   Map logical coordinates to the surface of a sphere
 * \details  This function will map logical coordinates in [0,1] to (x,y,z) coordinates
 *   on the surface of a sphere.  It uses the mapping by:
 *   Dona Calhoun, Christiane Helzel, Randall LeVeque, "Logically Rectangular Grids
 *      and Finite GeomType::Volume Methods for PDEs in Circular and Spherical Domains",
 *      SIAM REVIEW, Vol. 50, No. 4, pp. 723–752 (2008)
 * \param[in] R         Radius of sphere
 * \param[in] x         Logical x coordinate
 * \param[in] y         Logical y coordinate
 * @return              Returns the physical point
 */
Point3D map_logical_sphere_surface( double R, double x, double y );

/**
 * \brief   Map coordinates from the surface of a sphere to logical
 * \details  This function will map the (x,y,z) coordinates on the surface of a
 *   sphere to logical coordinates in [0,1].  It uses the mapping by:
 *   Dona Calhoun, Christiane Helzel, Randall LeVeque, "Logically Rectangular Grids
 *      and Finite GeomType::Volume Methods for PDEs in Circular and Spherical Domains",
 *      SIAM REVIEW, Vol. 50, No. 4, pp. 723–752 (2008)
 * \param[in] R         Radius of sphere
 * \param[in] x         Physical x coordinate
 * \param[in] y         Physical y coordinate
 * \param[in] z         Physical z coordinate
 * @return              Returns a pair with the logical (x,y) values
 */
std::pair<double, double> map_sphere_surface_logical( double R, double x, double y, double z );


/**
 * \brief   Map logical coordinates to a shell
 * \details  This function will map logical coordinates in [0,1] to (x,y,z) coordinates
 *   in a shell.  It uses the mapping by:
 *   Dona Calhoun, Christiane Helzel, Randall LeVeque, "Logically Rectangular Grids
 *      and Finite GeomType::Volume Methods for PDEs in Circular and Spherical Domains",
 *      SIAM REVIEW, Vol. 50, No. 4, pp. 723–752 (2008)
 * \param[in] r1        Inner radius of shell
 * \param[in] r2        Outer radius of shell
 * \param[in] x         Logical x coordinate
 * \param[in] y         Logical y coordinate
 * \param[in] z         Logical z coordinate
 * @return              Returns the physical point
 */
Point3D map_logical_shell( double r1, double r2, double x, double y, double z );


/**
 * \brief   Map a shell to logical coordinates
 * \details  This function will map (x,y,z) coordinates
 *   in a shell to logical coordinates in [0,1].  It uses the mapping by:
 *   Dona Calhoun, Christiane Helzel, Randall LeVeque, "Logically Rectangular Grids
 *      and Finite GeomType::Volume Methods for PDEs in Circular and Spherical Domains",
 *      SIAM REVIEW, Vol. 50, No. 4, pp. 723–752 (2008)
 * \param[in] r1        Inner radius of shell
 * \param[in] r2        Outer radius of shell
 * \param[in] x         Logical x coordinate
 * \param[in] y         Logical y coordinate
 * \param[in] z         Logical z coordinate
 * @return              Returns the physical point
 */
Point3D map_shell_logical( double r1, double r2, double x, double y, double z );


/**
 * \brief   Compute the intersection of a ray and a line segment
 * \details  This function will compute the intersection of a ray with a line segment.
 *    If the ray will never intersect the object, this distance is inf.
 * \param[in] pos       Starting point of ray
 * \param[in] ang       Direction of ray
 * \param[in] v1        First vertex
 * \param[in] v2        Second vertex
 * @return              Returns the distance
 */
double
distanceToLine( const Point2D &pos, const Point2D &ang, const Point2D &v1, const Point2D &v2 );


/**
 * \brief   Compute the intersection of a ray and a box
 * \details  This function will compute the intersection of a ray with a box.
 *    If the ray will never intersect the object, this distance is inf.
 * \param[in] pos       Starting point of ray
 * \param[in] ang       Direction of ray
 * \param[in] box       Box coordinates (x_min,x_max,y_min,y_max,z_min,z_max)
 * @return              Returns the distance
 */
double distanceToBox( const Point3D &pos, const Point3D &ang, const std::array<double, 6> &box );


/**
 * \brief   Compute the intersection of a ray and an infinite plane
 * \details  This function will compute the intersection of a ray with an infinite plane.
 *    The plane is described by the normal and a point on the plane.
 *    If the ray will never intersect the object, this distance is inf.
 * \param[in] n         Normal to the plane
 * \param[in] p0        A point on the plane
 * \param[in] pos       Starting point of ray
 * \param[in] ang       Direction of ray
 * @return              Returns the distance
 */
double
distanceToPlane( const Point3D &n, const Point3D &p0, const Point3D &pos, const Point3D &ang );


/**
 * \brief   Compute the intersection of a ray and circle (2D)
 * \details  This function will compute the intersection of a ray with a circle.
 *    It assumes a circle of radius r centered at the origin.
 *    If the ray is inside the cylinder the distance is negative.
 *    If the ray will never intersect the object, this distance is inf.
 * \param[in] r         Radius of circle
 * \param[in] pos       Starting point of ray
 * \param[in] ang       Direction of ray
 * @return              Returns the distance
 */
double distanceToCircle( double r, const Point2D &pos, const Point2D &ang );


/**
 * \brief   Compute the intersection of a ray and cylinder
 * \details  This function will compute the intersection of a ray with a cylinder.
 *    It assumes a cylinder of radius r and height h, centered at the origin
 *    along the z axis: z=[-h/2,h/2].  If the ray is inside the cylinder the distance
 *    is negative.  If the ray will never intersect the object, this distance is inf.
 * \param[in] r         Radius of cylinder
 * \param[in] h         Height of cylinder
 * \param[in] pos       Starting point of ray
 * \param[in] ang       Direction of ray
 * @return              Returns the distance
 */
double distanceToCylinder( double r, double h, const Point3D &pos, const Point3D &ang );


/**
 * \brief   Compute the intersection of a ray and tube
 * \details  This function will compute the intersection of a ray with a tube.
 *    It assumes a tube ith inner radius r1, outer radius r2, and height h
 *    centered at the origin along the z axis: z=[-h/2,h/2].
 *    If the ray is inside the cylinder the distance is negative.
 *    If the ray will never intersect the object, this distance is inf.
 * \param[in] r1        Inner radius of cylinder
 * \param[in] r2        Outer radius of cylinder
 * \param[in] h         Height of cylinder
 * \param[in] pos       Starting point of ray
 * \param[in] ang       Direction of ray
 * @return              Returns the distance
 */
double distanceToTube( double r1, double r2, double h, const Point3D &pos, const Point3D &ang );


/**
 * \brief   Compute the intersection of a ray and sphere
 * \details  This function will compute the intersection of a ray with a cylinder.
 *    It assumes a sphere of radius r and height h, centered at the origin.
 *    If the ray is inside the cylinder the distance is negative.
 *    If the ray will never intersect the object, this distance is inf.
 * \param[in] r         Radius of cylinder
 * \param[in] pos       Starting point of ray
 * \param[in] ang       Direction of ray
 * @return              Returns the distance
 */
double distanceToSphere( double r, const Point3D &pos, const Point3D &ang );


/**
 * \brief   Compute the intersection of a ray and cone
 * \details  This function will compute the intersection of a ray with a cone.
 *    It assumes a cone with the apex at the origin in the direction given by V.
 *    If the ray is inside the cone the distance is negative.
 *    If the ray will never intersect the object, this distance is inf.
 * \param[in] V         Direction of cone
 * \param[in] theta     Apex angle of cone
 * \param[in] pos       Starting point of ray
 * \param[in] ang       Direction of ray
 * @return              Returns the distance
 */
double distanceToCone( const Point3D &V, double theta, const Point3D &pos, const Point3D &ang );


/**
 * \brief   Compute the barycentric coordinates
 * \details  This function will compute the barycentric coordinates
 *    determine the normal.
 * \param[in] x         Verticies
 * \param[in] p         Point of interest
 * @return              Returns the barycentric coordinates
 */
template<int N1, int N2>
std::array<double, N1> barycentric( const std::array<double, N2> ( &x )[N1],
                                    const std::array<double, N2> &p );


/**
 * \brief   Compute the normal to a plane defined by 3 points
 * \details  This function will compute the normal to a plane
 *    defined by three points.  The order of the points will
 *    determine the normal.
 * \param[in] v1        First vertex
 * \param[in] v2        Second vertex
 * \param[in] v3        Third vertex
 * @return              Returns the normal
 */
Point3D normal( const Point3D &v1, const Point3D &v2, const Point3D &v3 );


/**
 * \brief   Find the nearest point to a line segment
 * \details  This function will compute the nearest point to a line segment in 2D.
 * \param[in] v1        First vertex
 * \param[in] v2        Second vertex
 * \param[in] p         Point of interest
 * @return              Returns the normal
 */
Point2D nearest( const Point2D &v1, const Point2D &v2, const Point2D &p );


/**
 * \brief   Find the nearest point to a line segment
 * \details  This function will compute the nearest point to a line segment in 3D.
 * \param[in] v1        First vertex
 * \param[in] v2        Second vertex
 * \param[in] p         Point of interest
 * @return              Returns the normal
 */
Point3D nearest( const Point3D &v1, const Point3D &v2, const Point3D &p );


/**
 * \brief   Find the nearest point to on a triangle
 * \details  This function will compute the nearest point to a triangle
 *    defined by three points in 3D.
 * \param[in] v         Verticies
 * \param[in] p         Point of interest
 * @return              Returns the normal
 */
Point3D nearest( const Point3D ( &v )[3], const Point3D &p );


/**
 * \brief   Subdivide a triangle
 * \details  Given the verticies of a triangle, sub-divide the triangle
 *     recursively until it is with the resolution, returning the new set of points
 * \param[in] v         Verticies
 * \param[in] res       Desired resolution
 * @return              Returns the new points (excluding the verticies)
 */
std::vector<AMP::Mesh::Point> subdivide( const std::array<AMP::Mesh::Point, 3> &v, double res );


//! Compute the normal to the plane formed by 3 points
Point3D normal( const Point3D &a, const Point3D &b, const Point3D &c );


} // namespace AMP::Geometry::GeometryHelpers

#endif
