#ifndef included_AMP_TriangleMesh
#define included_AMP_TriangleMesh

#include "AMP/ampmesh/Geometry.h"
#include "AMP/ampmesh/Mesh.h"
#include "AMP/ampmesh/MeshID.h"
#include "AMP/ampmesh/MeshIterator.h"

#ifdef USE_AMP_VECTORS
namespace AMP {
namespace LinearAlgebra {
class Vector;
}
} // namespace AMP
#endif

#include "AMP/utils/shared_ptr.h"
#include <array>
#include <map>
#include <vector>


namespace AMP {
namespace Mesh {


template<size_t NG, size_t NP>
class TriangleMeshIterator;
template<size_t NG, size_t NP>
class TriangleMeshElement;


/**
 * \class TriangleMesh
 * \brief A class used to represent an unstructured mesh of Triangles/Tetrahedrals
 */
template<size_t NG, size_t NP>
class TriangleMesh : public AMP::Mesh::Mesh
{
public: // Convenience typedefs
    typedef std::array<double, NP> Point;
    typedef std::array<ElementID, 2> Edge;
    typedef std::array<ElementID, 3> Triangle;
    typedef std::array<ElementID, 4> Tetrahedron;

public:
    /**
     * \brief Read in mesh files, partition domain, and prepare environment for simulation
     * \details Create triangle mesh data from the given parameters
     * \param params  Parameters for constructing a mesh from an input database
     */
    static AMP::shared_ptr<TriangleMesh<NG, NP>> generate( MeshParameters::shared_ptr params );

    /**
     * \brief Generate a triangle mesh from local triangle coordinates
     * \details  Create a triangle mesh from the local triangle coordinates.
     *    Note: Triangle list should be unique for each rank, load balance will be automatically
     * adjusted. \param triangles  List of triangles (each rank may contribute a unique list) \param
     * comm       Communicator to use (load balance wil be automatically generated on this comm)
     * \param tol        Relative tolerance (based on range of points) to use to determine if two
     * points are the same
     */
    static AMP::shared_ptr<TriangleMesh<NG, NP>>
    generate( const std::vector<std::array<Point, NG + 1>> &triangles,
              const AMP_MPI &comm,
              double tol = 1e-12 );


    //! Virtual function to copy the mesh (allows use to proply copy the derived class)
    virtual AMP::shared_ptr<Mesh> clone() const override final;


    /**
     * \brief   Estimate the number of elements in the mesh
     * \details  This function will estimate the number of elements in the mesh.
     *   This is used so that we can properly balance the meshes across multiple processors.
     *   Ideally this should be both an accurate estimate and very fast.  It should not require
     *   any communication and should not have to actually load a mesh.
     * \param params Parameters for constructing a mesh from an input database
     */
    static size_t estimateMeshSize( const MeshParameters::shared_ptr &params );

    /**
     * \brief   Estimate the number of elements in the mesh
     * \details  This function will estimate the number of elements in the mesh.
     *   This is used so that we can properly balance the meshes across multiple processors.
     *   Ideally this should be both an accurate estimate and very fast.  It should not require
     *   any communication and should not have to actually load a mesh.
     * \param params Parameters for constructing a mesh from an input database
     */
    static std::vector<size_t> estimateLogicalMeshSize( const MeshParameters::shared_ptr &params );


    /**
     * \brief   Return the maximum number of processors that can be used with the mesh
     * \details  This function will return the maximum number of processors that can
     *   be used with the mesh.
     * \param params Parameters for constructing a mesh from an input database
     */
    static size_t maxProcs( const MeshParameters::shared_ptr &params );


    // Copy/move constructors
    TriangleMesh( const TriangleMesh & );
    TriangleMesh( TriangleMesh && ) = default;
    TriangleMesh &operator=( const TriangleMesh & ) = delete;
    TriangleMesh &operator=( TriangleMesh && ) = default;

    //! Deconstructor
    virtual ~TriangleMesh();


    /* Return the number of local element of the given type
     * \param type   Geometric type
     */
    virtual size_t numLocalElements( const GeomType type ) const override final;


    /* Return the global number of elements of the given type
     * Note: depending on the mesh this routine may require global communication across the mesh.
     * \param type   Geometric type
     */
    virtual size_t numGlobalElements( const GeomType type ) const override final;


    /* Return the number of ghost elements of the given type on the current processor
     * \param type   Geometric type
     */
    virtual size_t numGhostElements( const GeomType type, const int gcw ) const override final;


    /**
     * \brief    Return an MeshIterator over the given geometric objects
     * \details  Return an MeshIterator over the given geometric objects
     * \param type   Geometric type to iterate over
     * \param gcw    Desired ghost cell width
     */
    virtual MeshIterator getIterator( const GeomType type, const int gcw = 0 ) const override final;


    /**
     * \brief    Return an MeshIterator over the given geometric objects on the surface
     * \details  Return an MeshIterator over the given geometric objects on the surface
     * \param type   Geometric type to iterate over
     * \param gcw    Desired ghost cell width
     */
    virtual MeshIterator getSurfaceIterator( const GeomType type,
                                             const int gcw = 0 ) const override final;


    /**
     * \brief    Return the list of all boundary ID sets in the mesh
     * \details  Return the list of all boundary ID sets in the mesh
     * Note: depending on the mesh this routine may require global communication across the mesh.
     */
    virtual std::vector<int> getBoundaryIDs() const override final;


    /**
     * \brief    Return an MeshIterator over the given geometric objects on the given boundary ID
     * set
     * \details  Return an MeshIterator over the given geometric objects on the given boundary ID
     * set
     * \param type   Geometric type to iterate over
     * \param id     Boundary id for the elements (example: sideset id)
     * \param gcw    Desired ghost cell width
     */
    virtual MeshIterator getBoundaryIDIterator( const GeomType type,
                                                const int id,
                                                const int gcw = 0 ) const override final;

    /**
     * \brief    Return the list of all boundary ID sets in the mesh
     * \details  Return the list of all boundary ID sets in the mesh
     * Note: depending on the mesh this routine may require global communication across the mesh.
     */
    virtual std::vector<int> getBlockIDs() const override final;


    /**
     * \brief    Return an MeshIterator over the given geometric objects on the given block ID set
     * \details  Return an MeshIterator over the given geometric objects on the given block ID set
     * \param type   Geometric type to iterate over
     * \param id     Block id for the elements (example: block id in cubit, subdomain in libmesh)
     * \param gcw    Desired ghost cell width
     */
    virtual MeshIterator
    getBlockIDIterator( const GeomType type, const int id, const int gcw = 0 ) const override final;


    /**
     * \brief    Return an MeshIterator constructed through a set operation of two other
     * MeshIterators.
     * \details  Return an MeshIterator constructed through a set operation of two other
     * MeshIterators.
     * \param OP Set operation to perform.
     *           SetOP::Union - Perform a union of the iterators ( A U B )
     *           SetOP::Intersection - Perform an intersection of the iterators ( A n B )
     *           SetOP::Complement - Perform a compliment of the iterators ( A - B )
     * \param A  Pointer to MeshIterator A
     * \param B  Pointer to MeshIterator B
     */
    static MeshIterator getIterator( SetOP OP, const MeshIterator &A, const MeshIterator &B );


    /**
     * \brief    Return a mesh element given it's id.
     * \details  This function queries the mesh to get an element given the mesh id.
     *    This function is only required to return an element if the id is local.
     *    Ideally, this should be done in O(1) time, but the implimentation is up to
     *    the underlying mesh.  The base class provides a basic implimentation, but
     *    uses mesh iterators and requires O(N) time on the number of elements in the mesh.
     * \param id    Mesh element id we are requesting.
     */
    virtual MeshElement getElement( const MeshElementID &id ) const override final;


    /**
     * \brief    Return the parent elements of the given mesh element
     * \details  This function queries the mesh to get an element given the mesh id,
     *    then returns the parent elements that have the element as a child
     * \param elem  Mesh element of interest
     * \param type  Element type of the parents requested
     */
    virtual std::vector<MeshElement> getElementParents( const MeshElement &elem,
                                                        const GeomType type ) const override final;

    /**
     * \brief    Is the mesh movable
     * \details  This function will check if the mesh can be displaced.
     * @return   enum indicating the extent the mesh can be moved
     */
    virtual Movable isMeshMovable() const override { return Movable::Deform; };


    /**
     * \brief    Displace the entire mesh
     * \details  This function will displace the entire mesh by a scalar value.
     *   This function is a blocking call for the mesh communicator, and requires
     *   the same value on all processors.  The displacement vector should be the
     *   size of the physical dimension.
     * \param x  Displacement vector
     */
    virtual void displaceMesh( const std::vector<double> &x ) override;


#ifdef USE_AMP_VECTORS
    /**
     * \brief    Displace the entire mesh
     * \details  This function will displace the entire mesh by displacing
     *   each node by the values provided in the vector.  This function is
     *   a blocking call for the mesh communicator
     * \param x  Displacement vector.  Must have N DOFs per node where N
     *           is the physical dimension of the mesh.
     */
    virtual void displaceMesh( AMP::shared_ptr<const AMP::LinearAlgebra::Vector> x ) override;
#endif


protected:
    // Constructors
    TriangleMesh();
    explicit TriangleMesh( MeshParameters::shared_ptr );
    explicit TriangleMesh( const std::vector<std::array<double, NP>> &verticies,
                           const std::vector<std::array<int64_t, NG + 1>> &triangles,
                           const std::vector<std::array<int64_t, NG + 1>> &tri_nab,
                           const AMP_MPI &comm );
    void initialize();

protected:
    // Return the IDs of the elements composing the current element
    void getElementsIDs( const ElementID &id, const GeomType type, ElementID *IDs ) const;
    inline void getVerticies( const ElementID &id, int &N, ElementID *IDs ) const;

    // Return the IDs of the neighboring elements
    void getNeighborIDs( const ElementID &id, std::vector<ElementID> &IDs ) const;

    // Return the IDs of the parent elements
    std::vector<ElementID> getElementParents( const ElementID &id, const GeomType type ) const;

    // Return the coordinated of the given vertex
    // Note: no error checking is done to make sure it is a valid vertex
    const Point &getPos( const ElementID &id ) const;

    // Check if the element is on the given boundry, block, etc
    bool isOnSurface( const ElementID &elemID ) const;
    bool isOnBoundary( const ElementID &elemID, int id ) const;
    bool isInBlock( const ElementID &elemID, int id ) const;
    static bool inIterator( const ElementID &id, const TriangleMeshIterator<NG, NP> &it );

    // Friends
    friend TriangleMeshIterator<NG, NP>;
    friend TriangleMeshElement<NG, NP>;

private: // Internal data
    // Store the locat start indicies
    std::array<size_t, 4> d_N_global;

    // Store the local triangle data
    std::vector<Point> d_vert;
    std::vector<Edge> d_edge;
    std::vector<Triangle> d_tri;
    std::vector<Tetrahedron> d_tet;
    std::vector<std::array<ElementID, NG + 1>> d_neighbors;

    // Store the ghost data
    std::map<ElementID, Point> d_remote_vert;
    std::map<ElementID, Edge> d_remote_edge;
    std::map<ElementID, Triangle> d_remote_tri;
    std::map<ElementID, Tetrahedron> d_remote_tet;

    // Store the parent data
    std::vector<size_t> d_parent_size[NG][NG];
    std::vector<size_t> d_parent_offset[NG][NG];
    std::vector<ElementID> d_parent_ids[NG][NG];

    // Store children data
    std::vector<std::array<ElementID, 3>> d_tri_edge;
    std::vector<std::array<ElementID, 4>> d_tet_tri;
    std::vector<std::array<ElementID, 6>> d_tet_edge;

    // Store common iterators
    std::vector<TriangleMeshIterator<NG, NP>> d_iterators;
    std::vector<TriangleMeshIterator<NG, NP>> d_surface_iterators;
    std::vector<std::vector<TriangleMeshIterator<NG, NP>>> d_boundary_iterators;
    std::vector<std::vector<TriangleMeshIterator<NG, NP>>> d_block_iterators;
};


} // namespace Mesh
} // namespace AMP


#endif
