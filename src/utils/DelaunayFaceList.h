#ifndef included_AMP_DelaunayFaceList
#define included_AMP_DelaunayFaceList

#include <stdint.h>
#include <stdlib.h>
#include <vector>


namespace AMP::DelaunayTessellation {


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
     *                  FaceList.
     * @param tri_id    The initial triangle id
     * @param new_tri   The triangle list (NDIM+1)
     */
    FaceList( const int N, const TYPE *x, const int tri_id, const int tri[], const TYPE TOL_VOL );

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
    int hash_table[1024]; // Internal hash table to improve performance when search for a given face
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

} // namespace AMP::DelaunayTessellation

#endif
