#ifndef included_AMP_DelaunayInterpolation
#define included_AMP_DelaunayInterpolation

#include <stdlib.h>

#include "AMP/utils/kdtree.h"


/** \class DelaunayInterpolation
 *
 * This class provides Delaunay based N-dimensional simplex interpolation.
 */
template<class TYPE>
class DelaunayInterpolation
{
public:
    //! Constructor that specifies the number of dimensions
    /*!
     * This is the intended constructor that initializes the class.  It has a single argument
     * consisting of the number of dimensions.
     * @param ndim              The number of dimensions
     */
    explicit DelaunayInterpolation( int ndim );


    // Deleted constructors
    DelaunayInterpolation()                                = delete;
    DelaunayInterpolation( const DelaunayInterpolation & ) = delete;
    DelaunayInterpolation &operator=( const DelaunayInterpolation & ) = delete;


    //! Empty destructor.
    ~DelaunayInterpolation();


    //! Function to set the desired level of data to be kept in the class
    /*!
     * This function allows the user to control the number of internal structures that are
     * saved between function calls.  Many of the computational routines need intermediate
     * data structures that can represent a significant fraction of the computational time.
     * By saving these structures between calls, this time can be saved at the cost of
     * additional memory usage.  Note: the level of storage can be changed at any time.
     * If it is increased, the internal structures will not be created until needed, if it
     * is decreased, the internal structures will be deleted immediately.  There is an exception
     * for the coordinates which will copied immediately.
     * @param level     Ammount of internal data that will be saved.  (default = 4)
     *                  1 - Save only the minimum amount of information to store the tesselation
     *                      This requires ~6*N ints in 2D and ~24*N ints in 3D
     *                  2 - Save the tesselation and the coordinates
     *                      This requires the storage of (1) + dim*N values for the coordinates
     *                      It is strongly recommended that the coordinates are stored internally.
     *                      If not, then the user is responsible for maintaining the coordinate
     * storage.
     *                      See update_coordinates for more information.
     *                  3 - Save the tesselation and common data for searching
     *                      This requires ~ 5x the storage (1).
     *                      This does NOT store the coordinates, see notes in (2).
     *                  4 - Save all data
     *                      This requires the most storage, but will give the best performance.
     *                      This works out to ~150 bytes/node and ~500 bytes/node in 3D.
     *                      It is recommended unless memory is an issue.
     */
    void set_storage_level( const int level = 4 );


    //! Function to get the current memory usage
    /*!
     * This function returns the current number of bytes in use by the structure, or estimates the
     * memory usage if the storage level is changed.
     * @param level     Optional argument that specifies a new storage level to use to estimate
     *                  the ammount of memory needed.  If it is 0, the current storage in use is
     * returned.
     *                  See set_storage_level for more information on the storage options
     */
    size_t memory_usage( const int level = 0 ) const;


    //! Function to return the number of triangles in the tessellation
    /*!
     * This function returns the number of triangles in the tessellation.
     */
    size_t get_N_tri() const;


    //! Function to return the triangles in the tessellation
    /*!
     * This function returns the triangles in the tessellation.
     * Note: this function will allocate memory the the user will have to delete
     * @param base      The base for the indexing of tri( 0: use 0-based indexing,  1: use 1-based
     * indexing)
     */
    int *get_tri( int base ) const;


    //! Function to return the triangles neignbors
    /*!
     * This function returns the triangles neighbors
     * Note: this function will allocate memory the the user will have to delete
     * @param base      The base for the indexing of tri( 0: use 0-based indexing,  1: use 1-based
     * indexing)
     */
    int *get_tri_nab( int base ) const;


    //! Function to construct the tessellation
    /*!
     * This function creates the tessellation using the given points.
     * If sucessful, this routine returns 0.
     * @param N         The number of verticies
     * @param x         The coordinates of the verticies( ndim x N )
     *                  Note: The coordinates are either copied or the pointer is copied depending
     * on
     *                  the storage level.  If the pointer is used, it is the user's responsibility
     * to
     *                  make sure the memory is not deleted, or to update the coordinate pointers
     * before
     *                  each call.  See update_coordinates for more information.
     */
    int create_tessellation( const int N, const TYPE x[] );


    //! Function to construct the tessellation
    /*!
     * This function creates the tessellation using the given points.
     * If sucessful, this routine returns 0.
     * Note: this version always copies the points, changing the storage level if necessary
     * @param N         The number of verticies
     * @param x         The coordinates of the x verticies
     * @param y         The coordinates of the y verticies (may be NULL for 1D)
     * @param z         The coordinates of the z verticies (may be NULL for 1D/2D)
     */
    int create_tessellation( const int N, const TYPE *x, const TYPE *y, const TYPE *z );


    //! Function to construct the tessellation using a given tessellation
    /*!
     * This function sets the internal tessellation to match a provided tessellation.
     * It does not check if the provided tessellation is valid.  It allows the user to
     * provide their own tesselation if desired.
     * If sucessful, this routine returns 0.
     * @param N         The number of verticies
     * @param x         The coordinates of the verticies( ndim x N )
     *                  Note: The coordinates are either copied or the pointer is copied depending
     * on
     *                  the storage level.  If the pointer is used, it is the user's responsibility
     * to
     *                  make sure the memory is not deleted, or to update the coordinate pointers
     * before
     *                  each call.  See update_coordinates for more information.
     * @param base      The base for the indexing of tri( 0: use 0-based indexing,  1: use 1-based
     * indexing)
     * @param N_tri     The number of simplexes (triangles in 2D)
     * @param tri       The tesselation( ndim+1 x N_tri )
     */
    void create_tessellation(
        const int N, const TYPE x[], const int base, const int N_tri, const int tri[] );


    //! Function to update the coordinates
    /*!
     * This function updates the coordinates pointer.  This function is NOT meant to allow for
     * moving
     * verticies.  If the storage scheme does not include storing the coordinates internally, then
     * it
     * is the user's responsibility to ensure that the original coordinates are not deleted or
     * moved.
     * This function provides an interface to allow the user to change the memory address of the
     * coordinates.
     * This should only be done by experienced users.  It is strongly recommended for the user to
     * use a
     * storage scheme that includes the coordinates.
     * @param x         The coordinates of the verticies( ndim x N )
     *                  Note: The coordinates are either copied or the pointer is copied depending
     * on
     *                  the storage level.  If the pointer is used, it is the user's responsibility
     * to
     *                  make sure the memory is not deleted, or to update the coordinate pointers
     * before
     *                  each call.
     */
    void update_coordinates( const TYPE x[] );


    //! Function to copy the tessellation
    /*!
     * This function copies the internal tessellation to a user-provided array.
     * @param N         Output: The number of verticies
     * @param x         Output: The coordinates of the verticies, should be a preallocated of size(
     * ndim x N ).
     *                  If x==NULL then the coordinates will not be returned
     * @param base      Input:  The base for the indexing of tri( 0: use 0-based indexing,  1: use
     * 1-based indexing)
     * @param N_tri     Output: The number of simplexes (triangles in 2D)
     * @param tri       Output: The tesselation on input, should be a preallocated of size( ndim+1 x
     * N_tri ).
     *                  If tri==NULL then the triangle data will not be returned
     */
    void copy_tessellation( int *N, TYPE *x, const int base, int *N_tri, int *tri ) const;


    //! Subroutine to find the nearest neighbor to a point
    /*!
     * This function finds the nearest neighbor to each point.
     * It is able to do perform the search in O(N^(1/ndim)) on average.
     * Note: this function requires the calculate of the node lists if they are not stored (see
     * set_storage_level)
     * @param Ni        The number of points in xi and yi
     * @param xi        Coordinates of the query points( ndim x Ni )
     * @param base      Input:  The base for the indexing of tri( 0: use 0-based indexing,  1: use
     * 1-based indexing)
     * @param nearest   Ouput index of the nearest neighbor
     */
    template<class TYPE2>
    void
    find_nearest( const unsigned int Ni, const TYPE2 xi[], const int base, unsigned int *nearest );


    //! Subroutine to find the triangle that contains the point
    /*!
     * This function finds the triangle that contains the given points.  This uses a simple search
     * method,
     * which may not be the most efficient.  It is O(N*Ni).
     * Note: this function requires the calculate of the node lists and triangle neighbor lists if
     * they are not stored (see set_storage_level)
     * @param Ni        The number of points in xi
     * @param xi        Coordinates of the query points( ndim x Ni )
     * @param base      The starting index value( 0: use 0-based indexing,  1: use 1-based indexing)
     * @param index     Ouput index of triangle containing the point (-1: Point is outside convex
     * hull, -2: Search failed)
     * @param extrap    If the point is outside the convex hull, return the nearest triangle instead
     * of -1
     */
    template<class TYPE2>
    void find_tri(
        const unsigned int Ni, const TYPE2 xi[], const int base, int *index, bool extrap = false );


    //! Subroutine to calculate the gradient at each node
    /*!
     * This function gets a list of the nodes that connect to each node
     * @param f         Function values at the verticies( ndim )
     * @param method    Gradient method to use
     *                  1 - Use a simple least squares method using only the local nodes.  This
     * method is
     *                      relatively fast, but only first order in the gradient, causing the
     * truncation
     *                      error of the interpolate to be O(x^3).  Note that it will still give a
     * better cubic
     *                      interpolation than MATLAB's cubic in griddata.
     *                  2 - Least squares method that solves a sparse system.  This method is second
     * order in the
     *                      gradient yielding an interpolant that has truncation error O(x^4), but
     * requires soving
     *                      a sparse 2nx2n system.  This can be reasonably fast for most systems,
     * but is memory
     *                      limited for large systems (it can grow as O(n^2*ndim^2)).
     *                      Currently this method is NOT implimented in the C++ version.
     *                 3 - Least squares method that uses the matrix computed by method 2 and block
     * Gauss-Seidel
     *                      iteration to improve the gradient calculated by method 1.  Usually only
     * a finite number of
     *                      iterations are needed to significantly reduce the error.  Technically,
     * this method has a
     *                      trunction error that is still O(x^3), but for most purposes, this term
     * is reduced so it is
     *                      less than the O(x^4) term.  For most systems, this method is not
     * significantly faster than
     *                      method 2, but it does not have the same memory limitations and is
     * O(n*ndim^2) in memory.
     *                      Typically 10-20 iterations are all that is necessary.
     *                 4 - This is the same as method 3, but does not store any internal data.  This
     * saves us from creating
     *                      a large temporary structure ~10*ndim^2*N at a cost of ~2x in
     * performance.
     * @param grad      (Output) Calculated gradient at the nodes( ndim x N )
     * @param n_it      Optional argument specifying the number of Gauss-Seidel iterations.  Only
     * used if method = 3 or 4.
     */
    void calc_node_gradient( const double *f, const int method, double *grad, const int n_it = 20 );


    //! Subroutine to perform nearest-neighbor interpoaltion
    /*!
     * This function performs nearest-neighbor interpoaltion.
     * @param f         Function values at the triangle verticies( 1 x N )
     * @param Ni        Number of points to perform the interpolation
     * @param xi        Coordinates of the query points( ndim x Ni )
     * @param nearest   The nearest-neighbor points (see find_nearest)( 1 x Ni)
     * @param fi        (Output) The interpolated function values at xi( 1 x Ni)
     */
    template<class TYPE2>
    void interp_nearest( const double f[],
                         const unsigned int Ni,
                         const TYPE2 xi[],
                         const unsigned int nearest[],
                         double *fi );


    //! Subroutine to perform linear interpoaltion
    /*!
     * This function performs linear interpoaltion.
     * If a valid triangle index is not given, NaN will be returned.
     * If extrap is false and the point is not within the triangle, NaN will be returned.
     * @param f         Function values at the triangle verticies( 1 x N )
     * @param Ni        Number of points to perform the interpolation
     * @param xi        Coordinates of the query points( ndim x Ni )
     * @param index     The index of the triangle containing the point (see find_tri)
     * @param fi        (Output) The interpolated function values at xi( 1 x Ni)
     * @param gi        (Optional Output) The interpolated gradient at xi( ndim x Ni )
     * @param extrap    Do we want to extrapolate from the current triangle
     *                  Note: extrapolating can incure large error if sliver triangles on the
     * boundary are present
     */
    template<class TYPE2>
    void interp_linear( const double f[],
                        const unsigned int Ni,
                        const TYPE2 xi[],
                        const int index[],
                        double *fi,
                        double *gi  = nullptr,
                        bool extrap = false );


    //! Subroutine to perform cubic interpoaltion
    /*!
     * This function performs cubic interpoaltion.
     * Note: If the point is not contained within a triangle NaN will be returned.
     * @param f         Function values at the triangle verticies( 1 x N )
     * @param g         Gradient of f(x) at the triangle verticies( ndim x N ) (see
     * calc_node_gradient if unknown)
     * @param Ni        Number of points to perform the interpolation
     * @param xi        Coordinates of the query points( ndim x Ni )
     * @param index     The index of the triangle containing the point (see find_tri)
     * @param fi        (Output) The interpolated function values at xi( 1 x Ni )
     * @param gi        (Optional Output) The interpolated gradient at xi( ndim x Ni )
     * @param extrap    Do we want to extrapolate from the current triangle
     *                  0: Do not extrapolate (NaNs will be used for points outside the domain)
     *                  1: Extrapolate using linear interpolation (using the nearest point and it's
     * gradient)
     *                  2: Extrapolate using quadratic interpolation (using linear extrapolation for
     * the gradient)
     */
    template<class TYPE2>
    void interp_cubic( const double f[],
                       const double g[],
                       const unsigned int Ni,
                       const TYPE2 xi[],
                       const int index[],
                       double *fi,
                       double *gi = nullptr,
                       int extrap = 0 );


    //! Subroutine to compute the Barycentric coordinates
    /**
     * This function computes the Barycentric coordinates and the matrix T to convert
     * from Barycentric coordinates to cartesian (x=T*L).
     * @param ndim  Number of dimensions
     * @param x     Coordinates of the triangle verticies( ndim x ndim+1 )
     * @param xi    Coordinates of the desired point( ndim x 1 )
     * @param L     (output) The Barycentric coordinates of the point( ndim+1 x 1 )
     */
    static void compute_Barycentric( const int ndim, const double *x, const double *xi, double *L );


private:
    //! Subroutine to get the list of nodes that are connected to each node
    void create_node_neighbors() const;

    //! Subroutine to get the list of triangle that are connected to each triangle
    void create_tri_neighbors() const;

    //! Subroutine to get the starting triangle for each node
    void create_node_tri() const;

    //! Subroutine to get the list of triangle that are connected to each triangle
    void create_kdtree() const;

    // Subroutine to perform cubic interpolation for a single point
    void interp_cubic_single( const double f[],
                              const double g[],
                              const double xi[],
                              const int index,
                              double &fi,
                              double *gi,
                              int extrap );


    // Internal Data
    unsigned char d_ndim; // The number of dimensions
    int d_level;          // The storage level
    unsigned int d_N;     // The number of verticies
    unsigned int d_N_tri; // The number of triangles
    int *d_tri;           // The triangles (ndim+1 x N_tri)
    TYPE *d_x; // Pointer to the coordinates (they may be stored inside or outside the class)
               // (ndim+1 x N)
    mutable size_t d_N_node_sum;        // The sum of the number of node neighbors
    mutable unsigned int *d_N_node;     // The number of neighbor nodes for each node (1xN)
    mutable unsigned int **d_node_list; // The list of neighbor nodes for each node (1xN)
    // Note: The first element points to an array of size N_node_sum
    mutable int *d_tri_nab;  // List of neighbor triangles for each triangle( ndim+1 x N_tri )
    mutable int *d_node_tri; // For each node, a triangle that contains that node (will be used as
                             // the starting triangle)
    mutable kdtree *d_tree;  // Nearest neighbor search tree
};

#endif
