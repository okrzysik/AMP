#ifndef included_AMP_DelaunayTessellation
#define included_AMP_DelaunayTessellation

#include <stdint.h>
#include <stdlib.h>
#include <vector>


/** \class DelaunayInterpolation
 *
 * This class provides Delaunay Tessellation.
 * Note:  This class is not threaded, but is thread safe.
 * Note:  Currently these functions have not been tested (and likely do not work) with more than
 *    2^31 elements per array.  This limits the number of triangles to 2^31/(ndim+1).
 *    In 3d, this limit is 536,870,912 triangles and ~65,000,000 nodes.
 *    In 2d, this limit is 715,827,882 triangles and 357,913,941 nodes.
 *    It should be a relativly simple matter to increase the number of triangle to 2^31 regradless
 * of dimension.
 *    Moving beyond 2^31 triangles will require a different storage for the triangle ids (tri_nab).
 *    Moving beyond 2^31 nodes will require all data structures to be 64-bit.
 *    Moving to 64-bit storage for the arrays will have a significant impact on memory (2x).
 *    Assumming ~500M triangles and 65M nodes in 3d, the current memory requirements is ~24 GB for
 * the internal
 *    data, with peak memory usage as high as 32GB.
 * Note:  MATLAB's delaunay3 command with 1M nodes in R3 in [0 1] requires ~200 MB to store tri,
 *    peak memory usage of ~3.5 GB (excluding MATLAB and x), and ~200 s to process (on euv).
 *    Using 10M nodes in R2 in [0 1] requires a peak memory usage of ~10GB and ~600 s to process
 * (on euv).
 *    Currently with 1M nodes in R3 in [0 1], the program requires a peak memory usage of ~350 MB
 * and ~ 135 s (on laptop).
 *    In 2D with 10M nodes in [0,1], it requires ~750 MB and ~110 s to process (on laptop).
 *    Running with 1M nodes in R3 in a structured grid takes ~756s (on laptop).  Most of this time
 * is spent in add_node.
 */
class DelaunayTessellation
{
public:
    //! Function that creates the Delaunay Tessellation
    /*!
     * This function will create a valid Delaunay Tessellation in multiple dimensions.
     * Currently only 2D and 3D are supported.  If successful, it will return the number of
     * triangles, if unsuccessful it will return an error code.  Additionally, there are
     * several optional stuctures.  The error codes are:
     *    -1:  The number of dimensions is not currently supported
     *    -2:  All the points are collinear (2D), or coplanar (3D)
     *    -3:  Error creating a new triangle on the convex hull
     *    -4:  Internal error with structures
     *    -5:  Unable to converge using the edge flips
     *    -6:  Unable to find a valid flip to fix the tessellation
     *    -7:  Duplicate (or nearly duclicate points detected)
     * @param ndim      The number of dimensions (currently only 2D and 3D are supported)
     * @param N         The number of verticies
     * @param x         The coordinates of the verticies (ndim x N)
     * @param tri       (output) The returned pointer where the triangles are stored.
     *                  Note: since it is not possible to know how many triangles will be created,
     *                  this function will allocate the memory for the triangles.  The user is
     *                  responsible for deallocating this memory using delete.
     * @param tri_nab  (output) The returned pointer where the triangle neighbors are stored.
     *                  Note: since it is not possible to know how many triangles will be created,
     *                  this function will allocate the memory for the triangles.
     *                  Note: if the triangle neighbors are not needed, this may be NULL.
     */
    static int create_tessellation(
        const int ndim, const int N, const double x[], int *tri[], int *tri_nab[] );
    static int create_tessellation(
        const int ndim, const int N, const short int x[], int *tri[], int *tri_nab[] );
    static int
    create_tessellation( const int ndim, const int N, const int x[], int *tri[], int *tri_nab[] );
    static int create_tessellation(
        const int ndim, const int N, const int64_t x[], int *tri[], int *tri_nab[] );


    //! Function to calculate the volume of a simplex
    /*!
     * This function calculates the volume of a N-dimensional simplex
     * Note: the sign of the volume depends on the order of the points.
     *   It will be positive for points stored in a clockwise mannor.
     * Note:  If the volume is zero, then the simplex is invalid.
     *   Eg. a line in 2D or a plane in 3D.
     * @param ndim      The number of dimensions (currently only 2D and 3D are supported)
     * @param x         The coordinates of the verticies of the simplex ( NDIM x NDIM+1 )
     */
    static double calc_volume( int ndim, const double x[] );


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
     * @param x         The coordinates of the verticies of the simplex
     * @param xi        The coordinates of the vertex to check
     * @param TOL_VOL   A tolerance on the volume to use
     */
    static int test_in_circumsphere( const int ndim,
                                     const double x[],
                                     const double xi[],
                                     const double TOL_VOL );
    static int
    test_in_circumsphere( const int ndim, const int x[], const int xi[], const double TOL_VOL );


    //! Function to return the circumsphere containing a simplex
    /*!
     * This function computes the circumsphere that contains a simplex
     * @param[in]  ndim     The number of dimensions
     * @param[in]  x        The coordinates of the verticies of the simplex
     * @param[out] R        The radius of the circumsphere
     * @param[out] c        The center of the circumsphere
     */
    static void get_circumsphere( const int ndim, const double x[], double &R, double *c );
    static void get_circumsphere( const int ndim, const int x[], double &R, double *c );


    //! Subroutine to compute the Barycentric coordinates
    /**
     * This function computes the Barycentric coordinates.
     * @param[in]  ndim     The number of dimensions
     * @param x     Coordinates of the triangle verticies ( NDIM x NDIM+1 )
     * @param xi    Coordinates of the desired point ( NDIM )
     * @param L     (output) The Barycentric coordinates of the point ( NDIM+1 )
     */
    static void compute_Barycentric( const int ndim, const double *x, const double *xi, double *L );


private:
    //! Empty constructor.
    DelaunayTessellation() {}

    //! Empty destructor.
    ~DelaunayTessellation() {}

    // Templated form of create_tessellation
    template<int NDIM, class TYPE, class ETYPE>
    static int create_tessellation( const int N, const TYPE x[], int *tri[], int *tri_nab[] );

    // Templated form of test_in_circumsphere
    template<int NDIM, class TYPE, class ETYPE>
    static int test_in_circumsphere( const TYPE x[], const TYPE xi[], const double TOL_VOL );

    // Templated form of get_circumsphere
    template<int NDIM, class TYPE>
    static void get_circumsphere( const TYPE x[], double &R, double *c );

    // Templated form of calc_volume
    template<int NDIM, class TYPE, class ETYPE>
    static double calc_volume( const TYPE x[] );

    // Templated form of compute_Barycentric
    template<int NDIM, class TYPE, class ETYPE>
    static void compute_Barycentric( const TYPE *x, const TYPE *xi, double *L );

    // Function to remove sliver triangles on the surface
    template<int NDIM, class TYPE, class ETYPE>
    static void
    clean_triangles( const int N, const TYPE *x, size_t &N_tri, int *tri, int *tri_nab );

    // Function to swap the indicies of two triangles
    template<int NDIM>
    static void swap_triangles( size_t N_tri, int i1, int i2, int *tri, int *tri_nab );


    // Structure to hold surfaces that need to be tested
    struct check_surface_struct {
        unsigned char test;
        unsigned char f1;
        unsigned char f2;
        int t1;
        int t2;
    };


    // Function to find a valid flip
    template<int NDIM, class TYPE, class ETYPE>
    static bool find_flip( const TYPE *x,
                           const int *tri,
                           const int *tri_nab,
                           const double TOL_VOL,
                           std::vector<check_surface_struct> &check_surface,
                           int &N_tri_old,
                           int &N_tri_new,
                           int *index_old,
                           int *new_tri,
                           int *new_tri_nab );
    template<class TYPE, class ETYPE>
    static inline bool find_flip_2D( const TYPE *x,
                                     const int *tri,
                                     const int *tri_nab,
                                     const double TOL_VOL,
                                     std::vector<check_surface_struct> &check_surface,
                                     int &N_tri_old,
                                     int &N_tri_new,
                                     int *index_old,
                                     int *new_tri,
                                     int *new_tri_nab );
    template<class TYPE, class ETYPE>
    static inline bool find_flip_3D( const TYPE *x,
                                     const int *tri,
                                     const int *tri_nab,
                                     const double TOL_VOL,
                                     std::vector<check_surface_struct> &check_surface,
                                     int &N_tri_old,
                                     int &N_tri_new,
                                     int *index_old,
                                     int *new_tri,
                                     int *new_tri_nab );


    //! Function to check if a simple flip is valid
    /*!
     * This function checks if a simple flip is valid.  To be valid, the line
     * between the two points that are not on the surface and the intersection
     * of the face must lie within the face.
     * @param x         The coordinates of the verticies of the simplex
     * @param i         The face to check
     * @param xi        The coordinates of the vertex to check
     */
    template<int NDIM, class TYPE, class ETYPE>
    static bool test_flip_valid( const TYPE x[], const int i, const TYPE xi[] );


    //! Function to perform a flip in 2D
    /*!
     * This function perform a 2-2 flip in 2D.  If sucessful, the function
     * will return true.  This function also needs to calculate the triangle neighbors.
     * Note: in 2D there is only 1 type of flip (the 2-2 flip).
     * @param N             The number of verticies
     * @param x             The coordinates of the verticies (2xN)
     * @param N_tri         The number of triangles
     * @param tri           The current list of triangles (3xN_tri)
     * @param tri_nab       The current list of triangle neighbors (3xN_tri)
     * @param t1            The index of the first triangle
     * @param s1            The surface index on the first triangle
     * @param t2            The index of the second triangle
     * @param s2            The surface index on the second triangle
     * @param index_old     (output) The list of triangles that have been replaced (1x2)
     * @param new_tri       (output) The new triangles that were created (3x2)
     * @param new_tri_nab   (output) The new triangles neighbors (3x2)
     *                          new_tri_nab >= 0   - Triangle neighbor is an existing triangle (with
     * the given index)
     *                          new_tri_nab -(i+1) - Triangle neighbor is the ith new triangle
     *                          new_tri_nab ==-1   - Triangle face is on the convex hull
     * @param TOL_VOL       The tolerance to use for determining if a simplex is valid
     */
    template<class TYPE, class ETYPE>
    static bool flip_2D( const TYPE x[],
                         const int tri[],
                         const int tri_nab[],
                         const int t1,
                         const int s1,
                         const int t2,
                         const int s2,
                         int *index_old,
                         int *new_tri,
                         int *new_tri_nab,
                         const double TOL_VOL );


    //! Function to perform a 2-2 flip in 3D
    /*!
     * This function performs a 2-2 flip in 3D.  If the flip can be applied, the function
     * will return true.  Note: a valid 2-2 flip in 3D requires that 4 of the
     * verticies are coplanar and on the convex hull.
     * In 3D the possible edge flips are:
     *   2-3 flip in which we break the 2 triangles into 3
     *   3-2 flip which is the reverse of the 2-3 flip
     *   4-4 flip in which we change 4 triangles into 4 new triangles and is
     *       necessary when 4 verticies are coplanar
     *   2-2 flip which is the 4-4 flip, but the coplanar verticies lie on the convex hull
     * @param N             The number of verticies
     * @param x             The coordinates of the verticies (2xN)
     * @param N_tri         The number of triangles
     * @param tri           The current list of triangles (3xN_tri)
     * @param tri_nab       The current list of triangle neighbors (3xN_tri)
     * @param t1            The index of the first triangle
     * @param s1            The surface index on the first triangle
     * @param t2            The index of the second triangle
     * @param s2            The surface index on the second triangle
     * @param index_old     (output) The list of triangles that have been replaced (1x2)
     * @param new_tri       (output) The new triangles that were created (3x2)
     * @param new_tri_nab   (output) The new triangles neighbors (3x2)
     *                          new_tri_nab >= 0   - Triangle neighbor is an existing triangle (with
     * the given index)
     *                          new_tri_nab -(i+1) - Triangle neighbor is the ith new triangle
     *                          new_tri_nab ==-1   - Triangle face is on the convex hull
     * @param TOL_VOL       The tolerance to use for determining if a simplex is valid
     */
    template<class TYPE, class ETYPE>
    static bool flip_3D_22( const TYPE x[],
                            const int tri[],
                            const int tri_nab[],
                            const int t1,
                            const int s1,
                            const int t2,
                            const int s2,
                            int *index_old,
                            int *new_tri,
                            int *new_tri_nab,
                            const double TOL_VOL );


    //! Function to perform a 2-3 flip in 3D
    /*!
     * This function performs a 2-3 flip in 3D.  If the flip can be applied, the function
     * will return true.
     * In 3D the possible edge flips are:
     *   2-3 flip in which we break the 2 triangles into 3
     *   3-2 flip which is the reverse of the 2-3 flip
     *   4-4 flip in which we change 4 triangles into 4 new triangles and is
     *       necessary when 4 verticies are coplanar
     *   2-2 flip which is the 4-4 flip, but the coplanar verticies lie on the convex hull
     * @param N             The number of verticies
     * @param x             The coordinates of the verticies (2xN)
     * @param N_tri         The number of triangles
     * @param tri           The current list of triangles (3xN_tri)
     * @param tri_nab       The current list of triangle neighbors (3xN_tri)
     * @param t1            The index of the first triangle
     * @param s1            The surface index on the first triangle
     * @param t2            The index of the second triangle
     * @param s2            The surface index on the second triangle
     * @param index_old     (output) The list of triangles that have been replaced (1x2)
     * @param new_tri       (output) The new triangles that were created (3x2)
     * @param new_tri_nab   (output) The new triangles neighbors (3x2)
     *                          new_tri_nab >= 0   - Triangle neighbor is an existing triangle (with
     * the given index)
     *                          new_tri_nab -(i+1) - Triangle neighbor is the ith new triangle
     *                          new_tri_nab ==-1   - Triangle face is on the convex hull
     * @param TOL_VOL       The tolerance to use for determining if a simplex is valid
     */
    template<class TYPE, class ETYPE>
    static bool flip_3D_23( const TYPE x[],
                            const int tri[],
                            const int tri_nab[],
                            const int t1,
                            const int s1,
                            const int t2,
                            const int s2,
                            int *index_old,
                            int *new_tri,
                            int *new_tri_nab,
                            const double TOL_VOL );


    //! Function to perform a 3-2 flip in 3D
    /*!
     * This function performs a 3-2 flip in 3D.  If the flip can be applied, the function
     * will return true.
     * In 3D the possible edge flips are:
     *   2-3 flip in which we break the 2 triangles into 3
     *   3-2 flip which is the reverse of the 2-3 flip
     *   4-4 flip in which we change 4 triangles into 4 new triangles and is
     *       necessary when 4 verticies are coplanar
     *   2-2 flip which is the 4-4 flip, but the coplanar verticies lie on the convex hull
     * @param N             The number of verticies
     * @param x             The coordinates of the verticies (2xN)
     * @param N_tri         The number of triangles
     * @param tri           The current list of triangles (3xN_tri)
     * @param tri_nab       The current list of triangle neighbors (3xN_tri)
     * @param t1            The index of the first triangle
     * @param s1            The surface index on the first triangle
     * @param t2            The index of the second triangle
     * @param s2            The surface index on the second triangle
     * @param index_old     (output) The list of triangles that have been replaced (1x2)
     * @param new_tri       (output) The new triangles that were created (3x2)
     * @param new_tri_nab   (output) The new triangles neighbors (3x2)
     *                          new_tri_nab >= 0   - Triangle neighbor is an existing triangle (with
     * the given index)
     *                          new_tri_nab -(i+1) - Triangle neighbor is the ith new triangle
     *                          new_tri_nab ==-1   - Triangle face is on the convex hull
     * @param TOL_VOL       The tolerance to use for determining if a simplex is valid
     */
    template<class TYPE, class ETYPE>
    static bool flip_3D_32( const TYPE x[],
                            const int tri[],
                            const int tri_nab[],
                            const int t1,
                            const int s1,
                            const int t2,
                            const int s2,
                            int *index_old,
                            int *new_tri,
                            int *new_tri_nab,
                            const double TOL_VOL );


    //! Function to perform a 4-4 flip in 3D
    /*!
     * This function performs a 4-4 flip in 3D.  If the flip can be applied, the function
     * will return true.  Note: a valid 4-4 flip in 3D requires 4 triangles to form
     * an octahedron with 4 verticies sharing a plane.
     * Note: there may be many valid flips (especially for a structured grid), but this function
     * will only return 1.
     * In 3D the possible edge flips are:
     *   2-3 flip in which we break the 2 triangles into 3
     *   3-2 flip which is the reverse of the 2-3 flip
     *   4-4 flip in which we change 4 triangles into 4 new triangles and is
     *       necessary when 4 verticies are coplanar
     *   2-2 flip which is the 4-4 flip, but the coplanar verticies lie on the convex hull
     * @param N             The number of verticies
     * @param x             The coordinates of the verticies (2xN)
     * @param N_tri         The number of triangles
     * @param tri           The current list of triangles (3xN_tri)
     * @param tri_nab       The current list of triangle neighbors (3xN_tri)
     * @param t1            The index of the first triangle
     * @param s1            The surface index on the first triangle
     * @param t2            The index of the second triangle
     * @param s2            The surface index on the second triangle
     * @param index_old     (output) The list of triangles that have been replaced (1x2)
     * @param new_tri       (output) The new triangles that were created (3x2)
     * @param new_tri_nab   (output) The new triangles neighbors (3x2)
     *                          new_tri_nab >= 0   - Triangle neighbor is an existing triangle (with
     * the given index)
     *                          new_tri_nab -(i+1) - Triangle neighbor is the ith new triangle
     *                          new_tri_nab ==-1   - Triangle face is on the convex hull
     * @param TOL_VOL       The tolerance to use for determining if a simplex is valid
     */
    template<class TYPE, class ETYPE>
    static bool flip_3D_44( const TYPE x[],
                            const int tri[],
                            const int tri_nab[],
                            const int t1,
                            const int s1,
                            const int t2,
                            const int s2,
                            int *index_old,
                            int *new_tri,
                            int *new_tri_nab,
                            const double TOL_VOL );

    //! Function to check the triangle and triangle neighbors
    /*!
     * This function performs a set of test taht checks that each triangle
     * is valid and that the triangle neighbors are valid.
     * Note: This can be an expensive test and should only be used when debugging
     * @param N             The number of verticies
     * @param x             The coordinates of the verticies (NDIMxN)
     * @param N_tri         The number of triangles
     * @param tri           The current list of triangles (NDIM+1xN_tri)
     * @param tri_nab       The current list of triangle neighbors (NDIM+1xN_tri)
     * @param N_unused      The number of unused coordinates
     * @param TOL_VOL       The tolerance to use for determining if a simplex is valid
     */
    template<int NDIM, class TYPE, class ETYPE>
    static bool check_current_triangles( int N,
                                         const TYPE x[],
                                         size_t N_tri,
                                         const int tri[],
                                         const int tri_nab[],
                                         const std::vector<size_t> &unused,
                                         double TOL_VOL );


    /*! Class for storing the faces on the convex hull
     *  Storing the list of faces on the convex hull is necessary, but requires maintaining
     *  a list of triangles and their faces.  This class simplifies this storage.
     *  In addition, this class will and new nodes to the convex hull, returning the triangles
     *  that were created.
     *  Note:  Seperate instantiations of this class are thread safe, but a single instance is not.
     */
    template<int NDIM, class TYPE, class ETYPE>
    class FaceList
    {
    public:
        /*! @brief  Standard constructor
         *  @detailed  Default constructor to be used
         * @param N         The number of verticies
         * @param x         The coordinates of the nodes (NDIM x N)
         *                  Note: coordinates must not change or be deleted during lifetime of
         * FaceList.
         * @param tri_id    The initial triangle id
         * @param new_tri   The triangle list (NDIM+1)
         */
        FaceList(
            const int N, const TYPE *x, const int tri_id, const int tri[], const TYPE TOL_VOL );

        //! Empty destructor.
        ~FaceList();

        //! Function to get the number of faces on the convex hull
        int get_N_face() { return N_face; }

        /*! @brief  Function to add a node to the convex hull
         *  @detailed  This function will add a new node to the convex hull.  In the process of
         *     of adding the node, some faces will be removed, new faces will be added,
         *     and new triangles will be generated.  This function will return the new triangles
         *     which must be added.  Note: the new triangles are not necessarilly Delaunay.
         *     This function returns the number of new triangles if successful, 0 if the node
         *     is inside the convex hull, and an error code for all other errors.
         *  Error codes:
         *     -1 - Invalid triangle created (eg. the triangle is flat in 3d or a line in 2D)
         *     -2 - All other errors
         * @param[in] node_id       The vertex index to add
         * @param[in/out] unused    A list of unused triangle ids
         * @param[inout]  id0       The first index to use for new triangle ids
         * @param[out] new_tri_id   A list of valid ids to use for new triangles
         * @param[out] new_tri      The list of new triangles that were created
         * @param[out] new_tri_nab  The list of triangle neighbors for the new triangles
         * @param[out] neighbor     The list of existing triangles that are neighbors to the new
         * triangles
         * @param[out] neighbor     The list of existing triangle faces that are neighbors to the
         * new triangles
         */
        int add_node( const int node_id,
                      std::vector<size_t> &unused,
                      size_t &N_tri,
                      unsigned int *new_tri_id,
                      int *new_tri,
                      int *new_tri_nab,
                      int *neighbor,
                      int *face_id );


        //! Function to update faces on the convex hull
        /*! This function will update faces on the convex hull.  This is necessary if the flips
         * changed the faces that lie on the convex hull.  There are two possibilites, the triangles
         * could have been changed but the faces are the same, in which case we need to update the
         * triangle numbers, face ids, and internal strucutures.  The second possiblity is the
         * entire
         * face configuration could change (eg a 2-2 flip in 3d).  This requires updating the faces
         * on
         * the convex hull.
         * Note: the number of faces on the convex hull should never change due to flips.
         * @param N         The number of faces that have changed
         * @param old_tid   The old triangle numbers (N)
         * @param old_fid   The old face ids (N)
         * @param new_tid   The new triangle numbers (N)
         * @param new_fid   The new face ids (N)
         * @param tri       The complete triangle list ( N_tri x NDIM+1 )
         */
        void update_face( const int N,
                          const int old_tid[],
                          const int old_fid[],
                          const int new_tid[],
                          const int new_fid[],
                          const int tri[] );

    private:
        // Private constructors
        FaceList();                              // Empty constructor.
        FaceList( const FaceList & );            // no implementation for copy
        FaceList &operator=( const FaceList & ); // no implementation for copy

        // Structure to store face information
        struct face_data_struct {
            int prev;
            int next;
            int tri_id;         // Triangle id
            int face_id;        // Face id
            int index[NDIM];    // Indicies of the face verticies
            TYPE x[NDIM][NDIM]; // Coordinates of the face verticies
            // Function to reset data to a NULL state
            void reset();
        };

        // Data members
        const int Nx; // The number of verticies
        int N_face;   // The number of faces on the convex hull
        int size;
        int hash_table[1024]; // Internal hash table to improve performance when search for a given
                              // face
        const TYPE *x0;       // The vertex coordinates
        double xc[NDIM];      // A point within the centroid
        double TOL_vol;       // Tolerance to use for calculation
        face_data_struct *data; // The stored data

        // Function that calculates the distance between a plane and a point
        double calc_surface_distance( const TYPE x[NDIM][NDIM], const TYPE xi[] ) const;
        bool outside_triangle( const TYPE x[NDIM][NDIM], const TYPE xi[] ) const;

        // Function to get a unique index for each face
        inline size_t get_face_index( int face, int tri ) { return face + tri * ( NDIM + 1 ); }

        // Function to delete a set of faces
        void delete_faces( int N_delete, int *ids );

        // Function to check that the internal data is valid
        void check_data();
    };
};

#endif
