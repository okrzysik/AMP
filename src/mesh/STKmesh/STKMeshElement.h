#ifndef included_AMP_STKMeshElement
#define included_AMP_STKMeshElement

#include "AMP/mesh/MeshElement.h"
#include "AMP/mesh/STKmesh/STKMesh.h"
#include "AMP/mesh/STKmesh/STKMeshIterator.h"
#include <memory>
#include <vector>

namespace AMP::Mesh {


/**
 * \class STKMeshElement
 * \brief A derived class used to define a mesh element
 * \details  This class provides routines for accessing and using a mesh element.
 * A mesh element can be thought of as the smallest unit of a mesh.  It is of a type
 * of GeomType.  This class is derived to store a STKMesh element.
 */
class STKMeshElement : public MeshElement
{
public:
    //! Empty constructor for a MeshElement
    STKMeshElement();

    //! Copy constructor
    STKMeshElement( const STKMeshElement & );

    //! Assignment operator
    STKMeshElement &operator=( const STKMeshElement & );

    //! De-constructor for a MeshElement
    virtual ~STKMeshElement();

    //! Return the element class
    virtual inline std::string elementClass() const { return "STKMeshElement"; }

    //! Return the elements composing the current element
    std::vector<MeshElement> getElements( const GeomType type ) const override;

    //! Return the elements neighboring the current element
    std::vector<std::shared_ptr<MeshElement>> getNeighbors() const override;

    //! Return the volume of the current element (does not apply to vertices)
    double volume() const override;

    //! Return the normal to the current element (does not apply to all elements)
    Point norm() const override;

    //! Return the coordinates of all vertices composing the element
    std::vector<double> coord() const override;

    /**
     * \brief     Return the centroid of the element
     * \details   This function returns the centroid of the element.  The
     *   centroid is defined as the average of the coordinates of the vertices.
     *   The centroid of a vertex is the vertex and will return the same result as coord().
     */
    std::vector<double> centroid() const override;

    /**
     * \brief     Return true if the element contains the point
     * \details   This function checks if the given point is inside or
     *   within TOL of the given element.  If the current element is a vertex,
     *   this function checks if the point is with TOL of the vertex.
     * \param pos   The coordinates of the point to check.
     * \param TOL   The tolerance to use for the computation.
     */
    bool containsPoint( const Point &pos, double TOL = 1e-12 ) const override;

    //! Check if the element is on the surface
    bool isOnSurface() const override;

    /**
     * \brief     Check if the current element is on the given boundary
     * \details   Check if the current element is on the boundary specified by the given id
     * \param id  The boundary id to check
     */
    bool isOnBoundary( int id ) const override;

    /**
     * \brief     Check if the current element is in the given block
     * \details   Check if the current element is in the block specified by the given id
     * \param id  The block id to check
     */
    bool isInBlock( int id ) const override;


protected:
    /** Default constructors
     * \param dim       Spatial dimension
     * \param element   Underlying STKmesh element
     * \param mesh      Underlying mesh
     * \param rank      Rank of the current processor (must agree with STKmesh->processor_id())
     * \param meshID    ID of the current mesh
     *        type      Element type
     */
    STKMeshElement( int dim,
                    stk::mesh::Entity *element,
                    unsigned int rank,
                    MeshID meshID,
                    const STKMesh *mesh );
    STKMeshElement( int dim,
                    std::shared_ptr<stk::mesh::Entity> element,
                    unsigned int rank,
                    MeshID meshID,
                    const STKMesh *mesh );

    //! Clone the iterator
    MeshElement *clone() const override;

    // Internal data
    int d_dim;                      // The dimension of the mesh
    unsigned int d_rank;            // The rank of the current processor
    stk::mesh::Entity *ptr_element; // The underlying STKmesh element properties (raw pointer)
    const STKMesh *d_mesh;          // The pointer to the current mesh
    MeshID d_meshID;                // The ID of the current mesh
    bool d_delete_elem;             // Do we need to delete the STKMesh element

    friend class AMP::Mesh::STKMesh;
    friend class AMP::Mesh::STKMeshIterator;

private:
    static constexpr uint32_t getTypeID() { return AMP::Utilities::hash_char( "STKMeshElement" ); }
};
} // namespace AMP::Mesh

#endif
