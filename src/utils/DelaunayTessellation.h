#ifndef included_AMP_DelaunayTessellation
#define included_AMP_DelaunayTessellation

#include <array>
#include <stdint.h>
#include <stdlib.h>
#include <tuple>
#include <vector>

#include "AMP/utils/Array.h"


namespace AMP {


/** \namespace DelaunayTessellation
 *
 * This namespace provides Delaunay Tessellation.
 * Note:  The functions are not threaded, but are thread safe.
 * Note:  Currently these functions have not been tested (and likely do not work) with more than
 *    2^31 elements per array.  This limits the number of triangles to 2^31/(ndim+1).
 *    In 3d, this limit is 536,870,912 triangles and ~65,000,000 nodes.
 *    In 2d, this limit is 715,827,882 triangles and 357,913,941 nodes.
 *    It should be a relatively simple matter to increase the number of triangle to 2^31
 *        regardless of dimension.
 *    Moving beyond 2^31 triangles will require a different storage for the triangle ids (tri_nab).
 *    Moving beyond 2^31 nodes will require all data structures to be 64-bit.
 *    Moving to 64-bit storage for the arrays will have a significant impact on memory (2x).
 *    Assuming ~500M triangles and 65M nodes in 3d, the current memory requirements is ~24 GB for
 *        the internal data, with peak memory usage as high as 32GB.
 * Note:  MATLAB's delaunay3 command with 1M nodes in R3 in [0 1] requires ~200 MB to store tri,
 *    peak memory usage of ~3.5 GB (excluding MATLAB and x), and ~200 s to process (on euv).
 *    Using 10M nodes in R2 in [0 1] requires a peak memory usage of ~10GB and ~600 s
 *       to process (oneuv).
 *    Currently with 1M nodes in R3 in [0 1], the program requires a peak memory usage of ~350 MB
 *        and ~ 135 s (on laptop).
 *    In 2D with 10M nodes in [0,1], it requires ~750 MB and ~110 s to process (on laptop).
 *    Running with 1M nodes in R3 in a structured grid takes ~756s (on laptop).
 *    Most of this time is spent in add_node.
 */
namespace DelaunayTessellation {


//! Function that creates the Delaunay Tessellation
/*!
 * This function will create a valid Delaunay Tessellation in multiple dimensions.
 * Currently only 2D and 3D are supported.  If successful, it will return the number of
 * triangles, if unsuccessful it will throw a std::exception.  Additionally, there are
 * several optional stuctures.
 * @param x         The coordinates of the vertices (ndim x N)
 * @return          Returns the triangles and triangle neighbors <tri,tri_nab>
 *                  tri - The returned pointer where the triangles are stored (ndim+1,N)
 *                  tri_nab - The returned pointer where the triangle neighbors are stored
 *                      (ndim+1,N)
 */
template<class TYPE>
std::tuple<AMP::Array<int>, AMP::Array<int>> create_tessellation( const Array<TYPE> &x );


//! Function to calculate the volume of a simplex
/*!
 * This function calculates the volume of a N-dimensional simplex
 * Note: the sign of the volume depends on the order of the points.
 *   It will be positive for points stored in a clockwise mannor.
 * Note:  If the volume is zero, then the simplex is invalid.
 *   Eg. a line in 2D or a plane in 3D.
 * @param ndim      The number of dimensions (currently only 2D and 3D are supported)
 * @param x         The coordinates of the vertices of the simplex ( NDIM x NDIM+1 )
 */
double calc_volume( int ndim, const double x[] );


//! Function to check if a point is inside the circumsphere of a simplex.
/*!
 * This function checks if a point is inside the circumsphere of a simplex.
 * It returns -1 if the point is outside the circumsphere, 1 if it is inside the sphere,
 * and 0 if it is within the tolerance of the sphere.
 * Note:  For this function to work properly, the volume of the simplex
 *    (as computed by calc_volume) must be positive.
 * Note:  If we are checking the surface between 2 simplicies and they are both valid
 *    (have a positive, non-zero volume), it is suffcient to check the vertix of 1 volume
 *    against the circumcircle of the other.  We do not need to perform both checks.
 * @param ndim      The number of dimensions
 * @param x         The coordinates of the vertices of the simplex
 * @param xi        The coordinates of the vertex to check
 * @param TOL_VOL   A tolerance on the volume to use
 */
int test_in_circumsphere( const int ndim,
                          const double x[],
                          const double xi[],
                          const double TOL_VOL );
int test_in_circumsphere( const int ndim, const int x[], const int xi[], const double TOL_VOL );


//! Function to return the circumsphere containing a simplex
/*!
 * This function computes the circumsphere that contains a simplex
 * @param[in]  ndim     The number of dimensions
 * @param[in]  x        The coordinates of the vertices of the simplex
 * @param[out] R        The radius of the circumsphere
 * @param[out] c        The center of the circumsphere
 */
void get_circumsphere( const int ndim, const double x[], double &R, double *c );
void get_circumsphere( const int ndim, const int x[], double &R, double *c );


//! Subroutine to compute the Barycentric coordinates
/**
 * This function computes the Barycentric coordinates.
 * @param[in]  ndim     The number of dimensions
 * @param x     Coordinates of the triangle vertices ( NDIM x NDIM+1 )
 * @param xi    Coordinates of the desired point ( NDIM )
 * @param L     (output) The Barycentric coordinates of the point ( NDIM+1 )
 */
void compute_Barycentric( const int ndim, const double *x, const double *xi, double *L );


} // namespace DelaunayTessellation
} // namespace AMP

#endif
