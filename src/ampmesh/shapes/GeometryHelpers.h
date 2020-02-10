#ifndef included_AMP_GeometryHelpers
#define included_AMP_GeometryHelpers


#include "AMP/ampmesh/MeshPoint.h"
#include <set>


namespace AMP {
namespace Geometry {
namespace GeometryHelpers {


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
AMP::Mesh::Point map_logical_sphere( double R, double x, double y, double z );


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
AMP::Mesh::Point map_sphere_logical( double R, double x, double y, double z );


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
AMP::Mesh::Point map_logical_sphere_surface( double R, double x, double y );

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
AMP::Mesh::Point map_logical_shell( double r1, double r2, double x, double y, double z );


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
AMP::Mesh::Point map_shell_logical( double r1, double r2, double x, double y, double z );


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
double distanceToPlane( const AMP::Mesh::Point &n,
                        const AMP::Mesh::Point &p0,
                        const AMP::Mesh::Point &pos,
                        const AMP::Mesh::Point &ang );


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
double distanceToCircle( double r, const AMP::Mesh::Point &pos, const AMP::Mesh::Point &ang );


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
double
distanceToCylinder( double r, double h, const AMP::Mesh::Point &pos, const AMP::Mesh::Point &ang );


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
double distanceToSphere( double r, const AMP::Mesh::Point &pos, const AMP::Mesh::Point &ang );


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
double distanceToCone( const AMP::Mesh::Point &V,
                       double theta,
                       const AMP::Mesh::Point &pos,
                       const AMP::Mesh::Point &ang );


} // namespace GeometryHelpers
} // namespace Geometry
} // namespace AMP

#endif
