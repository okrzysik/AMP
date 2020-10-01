#ifndef included_AMP_libmeshMeshElement
#define included_AMP_libmeshMeshElement


#include "AMP/ampmesh/MeshElement.h"
#include "AMP/ampmesh/libmesh/libmeshElemIterator.h"
#include "AMP/ampmesh/libmesh/libmeshMesh.h"
#include "AMP/ampmesh/libmesh/libmeshNodeIterator.h"
#include <memory>
#include <vector>


namespace AMP {
namespace Mesh {


/**
 * \class libmeshMeshElement
 * \brief A derived class used to define a mesh element
 * \details  This class provides routines for accessing and using a mesh element.
 * A mesh element can be thought of as the smallest unit of a mesh.  It is of a type
 * of GeomType.  This class is derived to store a libMesh element.
 */
class libmeshMeshElement : public MeshElement
{
public:
    //! Empty constructor for a MeshElement
    libmeshMeshElement();

    //! Copy constructor
    libmeshMeshElement( const libmeshMeshElement & );

    //! Assignment operator
    libmeshMeshElement &operator=( const libmeshMeshElement & );

    //! De-constructor for a MeshElement
    virtual ~libmeshMeshElement();

    //! Return the unique global ID of the element
    MeshElementID globalID() const override { return d_globalID; }

    //! Return the element class
    inline std::string elementClass() const override { return "libmeshMeshElement"; }

    //! Return the elements composing the current element
    virtual void getElements( const GeomType type,
                              std::vector<MeshElement> &elements ) const override;

    //! Return the elements neighboring the current element
    void getNeighbors( std::vector<MeshElement::shared_ptr> &neighbors ) const override;

    //! Return the volume of the current element (does not apply to verticies)
    double volume() const override;

    //! Return the normal to the current element (does not apply to all elements)
    Point norm() const override;

    //! Return the coordinates of all verticies composing the element
    Point coord() const override;

    /**
     * \brief     Return the centroid of the element
     * \details   This function returns the centroid of the element.  The
     *   centroid is defined as the average of the coordinates of the verticies.
     *   The centroid of a vertex is the vertex and will return the same result as coord().
     */
    Point centroid() const override;

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

    //! Return the owner rank according to AMP_COMM_WORLD
    unsigned int globalOwnerRank() const override;


protected:
    /** Default constructors
     * \param dim       Spatial dimension
     * \param type      Element type
     * \param element   Underlying libmesh element
     * \param mesh      Underlying mesh
     * \param rank      Rank of the current processor (must agree with libmesh->processor_id())
     * \param meshID    ID of the current mesh
     */
    libmeshMeshElement( int dim,
                        GeomType type,
                        void *element,
                        unsigned int rank,
                        MeshID meshID,
                        const libmeshMesh *mesh );
    libmeshMeshElement( int dim,
                        GeomType type,
                        std::shared_ptr<libMesh::Elem> element,
                        unsigned int rank,
                        MeshID meshID,
                        const libmeshMesh *mesh );

    //! Clone the iterator
    MeshElement *clone() const override;

    // Internal data
    int d_dim;                           // The dimension of the mesh
    unsigned int d_rank;                 // The rank of the current processor
    void *ptr_element;                   // The underlying libmesh element properties (raw pointer)
    std::shared_ptr<libMesh::Elem> ptr2; // Optional smart pointer to the element (to hold a copy)
    const libmeshMesh *d_mesh;           // The pointer to the current mesh
    MeshID d_meshID;                     // The ID of the current mesh
    bool d_delete_elem;                  // Do we need to delete the libMesh element

    friend class AMP::Mesh::libmeshMesh;
    friend class AMP::Mesh::libmeshNodeIterator;
    friend class AMP::Mesh::libmeshElemIterator;

private:
    static constexpr uint32_t getTypeID()
    {
        return AMP::Utilities::hash_char( "libmeshMeshElement" );
    }
    MeshElementID d_globalID;
};


} // namespace Mesh
} // namespace AMP

#endif
