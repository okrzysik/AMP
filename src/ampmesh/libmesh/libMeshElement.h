#ifndef included_AMP_libMeshElement
#define included_AMP_libMeshElement


#include "ampmesh/MeshElement.h"
#include "ampmesh/libmesh/libMesh.h"
#include "ampmesh/libmesh/libMeshIterator.h"
#include "utils/shared_ptr.h"
#include <vector>


namespace AMP {
namespace Mesh {


/**
 * \class libMeshElement
 * \brief A derived class used to define a mesh element
 * \details  This class provides routines for accessing and using a mesh element.
 * A mesh element can be thought of as the smallest unit of a mesh.  It is of a type
 * of GeomType.  This class is derived to store a libMesh element.
 */
class libMeshElement : public MeshElement
{
public:
    //! Empty constructor for a MeshElement
    libMeshElement();

    //! Copy constructor
    libMeshElement( const libMeshElement & );

    //! Assignment operator
    libMeshElement &operator=( const libMeshElement & );

    //! De-constructor for a MeshElement
    virtual ~libMeshElement();

    //! Return the element class
    virtual inline std::string elementClass() const override{ return "libMeshElement"; }

    //! Return the elements composing the current element
    virtual std::vector<MeshElement> getElements( const GeomType type ) const override;

    //! Return the elements neighboring the current element
    virtual std::vector<MeshElement::shared_ptr> getNeighbors() const override;

    //! Return the volume of the current element (does not apply to verticies)
    virtual double volume() const override;

    //! Return the coordinates of all verticies composing the element
    virtual std::vector<double> coord() const override;

    /**
     * \brief     Return the coordinate of the vertex
     * \details   This function returns the coordinates of the vertex
     *   in the given direction (only applies to verticies).
     * \param i     The direction requested.  Equivalent to coord()[i]
     */
    virtual double coord( int i ) const override;

    /**
     * \brief     Return the centroid of the element
     * \details   This function returns the centroid of the element.  The
     *   centroid is defined as the average of the coordinates of the verticies.
     *   The centroid of a vertex is the vertex and will return the same result as coord().
     */
    virtual std::vector<double> centroid() const override;

    /**
     * \brief     Return true if the element contains the point
     * \details   This function checks if the given point is inside or
     *   within TOL of the given element.  If the current element is a vertex,
     *   this function checks if the point is with TOL of the vertex.
     * \param pos   The coordinates of the point to check.
     * \param TOL   The tolerance to use for the computation.
     */
    virtual bool containsPoint( const std::vector<double> &pos, double TOL = 1e-12 ) const override;

    //! Check if the element is on the surface
    virtual bool isOnSurface() const override;

    /**
     * \brief     Check if the current element is on the given boundary
     * \details   Check if the current element is on the boundary specified by the given id
     * \param id  The boundary id to check
     */
    virtual bool isOnBoundary( int id ) const override;

    /**
     * \brief     Check if the current element is in the given block
     * \details   Check if the current element is in the block specified by the given id
     * \param id  The block id to check
     */
    virtual bool isInBlock( int id ) const override;

    //! Return the owner rank according to AMP_COMM_WORLD
    virtual unsigned int globalOwnerRank() const override;


protected:
    /** Default constructors
     * \param dim       Spatial dimension
     * \param type      Element type
     * \param element   Underlying libmesh element
     * \param mesh      Underlying mesh
     * \param rank      Rank of the current processor (must agree with libmesh->processor_id())
     * \param meshID    ID of the current mesh
     */
    libMeshElement( int dim,
                    GeomType type,
                    void *element,
                    unsigned int rank,
                    MeshID meshID,
                    const libMesh *mesh );
    libMeshElement( int dim,
                    GeomType type,
                    AMP::shared_ptr<::Elem>
                        element,
                    unsigned int rank,
                    MeshID meshID,
                    const libMesh *mesh );

    //! Clone the iterator
    virtual MeshElement *clone() const override;

    // Internal data
    int d_dim;                    // The dimension of the mesh
    unsigned int d_rank;          // The rank of the current processor
    void *ptr_element;            // The underlying libmesh element properties (raw pointer)
    AMP::shared_ptr<::Elem> ptr2; // Optional smart pointer to the element (to hold a copy)
    const libMesh *d_mesh;        // The pointer to the current mesh
    MeshID d_meshID;              // The ID of the current mesh
    bool d_delete_elem;           // Do we need to delete the libMesh element

    friend class AMP::Mesh::libMesh;
    friend class AMP::Mesh::libMeshIterator;
};
}
}

#endif
