#ifndef included_AMP_libMeshElement
#define included_AMP_libMeshElement

#include <vector>
#include <boost/shared_ptr.hpp>
#include "ampmesh/Mesh.h"
#include "ampmesh/MeshElement.h"
#include "ampmesh/libmesh/libMesh.h"
#include "ampmesh/libmesh/libMeshIterator.h"

// libMesh includes
#include "libmesh/elem.h"

namespace AMP {
namespace Mesh {


/**
 * \class libMeshElement
 * \brief A derived class used to define a mesh element
 * \details  This class provides routines for accessing and using a mesh element.
 * A mesh element can be thought of as the smallest unit of a mesh.  It is of a type
 * of GeomType.  This class is derived to store a libMesh element.
 */
class libMeshElement: public MeshElement
{
public:

    //! Empty constructor for a MeshElement
    libMeshElement ( );

    //! Copy constructor
    libMeshElement(const libMeshElement&);

    //! Assignment operator
    libMeshElement& operator=(const libMeshElement&);

    //! De-constructor for a MeshElement
    virtual ~libMeshElement ( );

    //! Return the elements composing the current element
    virtual std::vector<MeshElement> getElements(const GeomType type) const;

    //! Return the elements neighboring the current element
    virtual std::vector< MeshElement::shared_ptr > getNeighbors() const;

    //! Return the volume of the current element (does not apply to verticies)
    virtual double volume() const;

    //! Return the coordinates of all verticies composing the element
    virtual std::vector<double> coord() const;

    /**
     * \brief     Return the centroid of the element
     * \details   This function returns the centroid of the element.  The
     *   centroid is defined as the average of the coordinates of the verticies.
     *   The centroid of a vertex is the vertex and will return the same result as coord().
     */
    virtual std::vector<double> centroid() const;

    /**
     * \brief     Return true if the element contains the point
     * \details   This function checks if the given point is inside or 
     *   within TOL of the given element.  If the current element is a vertex,
     *   this function checks if the point is with TOL of the vertex.
     * \param pos   The coordinates of the point to check.
     * \param TOL   The tolerance to use for the computation.
     */
    virtual bool containsPoint( const std::vector<double> &pos, double TOL=1e-12 ) const;

    //! Check if the element is on the surface
    virtual bool isOnSurface() const;

    /**
     * \brief     Check if the current element is on the given boundary
     * \details   Check if the current element is on the boundary specified by the given id
     * \param id  The boundary id to check
     */
    virtual bool isOnBoundary(int id) const;

    /**
     * \brief     Check if the current element is in the given block
     * \details   Check if the current element is in the block specified by the given id
     * \param id  The block id to check
     */
    virtual bool isInBlock(int id) const;


protected:

    /** Default constructors
     * \param dim       Spatial dimension
     * \param type      Element type
     * \param element   Underlying libmesh element
     * \param mesh      Underlying mesh
     * \param rank      Rank of the current processor (must agree with libmesh->processor_id())
     * \param meshID    ID of the current mesh
     */
    libMeshElement(int dim, GeomType type, void* element, unsigned int rank, MeshID meshID, const libMesh* mesh );
    libMeshElement(int dim, GeomType type, boost::shared_ptr< ::Elem > element, unsigned int rank, MeshID meshID, const libMesh* mesh );

    //! Clone the iterator
    virtual MeshElement* clone() const;

    // Internal data
    int d_dim;                  // The dimension of the mesh
    unsigned int d_rank;        // The rank of the current processor
    void* ptr_element;          // The underlying libmesh element properties (raw pointer)
    boost::shared_ptr< ::Elem> ptr2; // Optional smart pointer to the element (to hold a copy)
    const libMesh* d_mesh;      // The pointer to the current mesh
    MeshID d_meshID;            // The ID of the current mesh
    bool d_delete_elem;         // Do we need to delete the libMesh element

    friend class AMP::Mesh::libMesh;
    friend class AMP::Mesh::libMeshIterator;

};



}
}

#endif

