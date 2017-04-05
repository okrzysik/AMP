#ifndef included_AMP_BoxMeshHelpers
#define included_AMP_BoxMeshHelpers

#include "ampmesh/structured/BoxMeshHelpers.h"

#include <set>
#include <tuple>


namespace AMP {
namespace Mesh {
namespace BoxMeshHelpers {



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
std::pair<double,double> map_logical_circle( double R, int method, double x, double y );


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
std::pair<double,double> map_circle_logical( double R, int method, double x, double y );


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
 * @return              Returns a tuple with the (x,y,z) value
 */
std::tuple<double,double,double> map_logical_sphere( double R, double x, double y, double z );


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
 * @return              Returns a tuple with the (x,y,z) value
 */
std::tuple<double,double,double> map_logical_sphere_surface( double R, double x, double y );

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
std::pair<double,double> map_sphere_surface_logical( double R, double x, double y, double z );


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
 * @return              Returns a tuple with the (x,y,z) value
 */
std::tuple<double,double,double> map_logical_shell( double r1, double r2, double x, double y, double z );


} // BoxMeshHelpers namespace
} // Mesh namespace
} // AMP namespace

#endif
