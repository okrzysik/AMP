#ifndef included_AMP_MeshElement
#define included_AMP_MeshElement

#include <vector>
#include <boost/shared_ptr.hpp>
#include "ampmesh/MeshID.h"


namespace AMP {
namespace Mesh {


/**
 * \class Mesh
 * \brief A class used to define a mesh element
 * \details  This class provides routines for accessing and using a mesh element.
 * A mesh element can be thought of as the smallest unit of a mesh.  It is of a type
 * of GeomType.  It contains the composing pieces of the element
 */
class MeshElement
{
public:
    /**
     *\typedef shared_ptr
     *\brief  Name for the shared pointer.
     *\details  Use this typedef for a reference counted pointer to a mesh manager object.
     */
    typedef boost::shared_ptr<MeshElement>  shared_ptr;

    //! Empty constructor for a MeshElement
    MeshElement ( );

    //! Copy constructor
    MeshElement(const MeshElement&);

    //! Assignment operator
    MeshElement& operator=(const MeshElement&);

    //! De-constructor for a MeshElement
    virtual ~MeshElement ( );

    // Overload operators
    inline bool operator== (const MeshElement& rhs ) const { return d_globalID == rhs.d_globalID; }
    inline bool operator!= (const MeshElement& rhs ) const { return d_globalID != rhs.d_globalID; }
    inline bool operator<  (const MeshElement& rhs ) const { return d_globalID <  rhs.d_globalID; }
    inline bool operator>  (const MeshElement& rhs ) const { return d_globalID >  rhs.d_globalID; }
    inline bool operator<= (const MeshElement& rhs ) const { return d_globalID <= rhs.d_globalID; }
    inline bool operator>= (const MeshElement& rhs ) const { return d_globalID >= rhs.d_globalID; }
    inline bool operator== (const MeshElementID& rhs ) const { return d_globalID == rhs; }
    inline bool operator!= (const MeshElementID& rhs ) const { return d_globalID != rhs; }
    inline bool operator<  (const MeshElementID& rhs ) const { return d_globalID <  rhs; }
    inline bool operator>  (const MeshElementID& rhs ) const { return d_globalID >  rhs; }
    inline bool operator<= (const MeshElementID& rhs ) const { return d_globalID <= rhs; }
    inline bool operator>= (const MeshElementID& rhs ) const { return d_globalID >= rhs; }

    //! Return the element type
    virtual GeomType elementType() const { return d_elementType; }

    //! Return the unique global ID of the element
    virtual MeshElementID globalID() const { return d_globalID; }

    //! Return the elements composing the current element
    virtual std::vector<MeshElement> getElements(const GeomType type) const;

    /**
     *  Return the elements neighboring the current element.
     *  One neighbor is returned for each side of the element.
     *  If the side is on the surface, then it's neighbor is null.
     *  For Verticies, a list of all verticies that share an element is returned.
     *  This list is in unsorted order.
     */  
    virtual std::vector< MeshElement::shared_ptr >  getNeighbors() const;

    //! Return the volume of the current element (does not apply to verticies)
    virtual double volume() const;

    //! Return the coordinates of the vertex (only applies to verticies)
    virtual std::vector<double> coord() const;

    /**
     * \brief     Return the centroid of the element
     * \details   This function returns the centroid of the element.  The
     *   centroid is defined as the average of the coordinates of the verticies.
     *   The centroid of a vertex is the vertex and will return the same result as coord().
     */
    virtual std::vector<double> centroid() const;

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

    //! The MeshElement type
    GeomType d_elementType;

    //! The unique global id of the element
    MeshElementID d_globalID;

    // A pointer to the derived class
    MeshElement *element;

    // Clone the iterator
    virtual MeshElement* clone() const;

    // Unique (per class) ID for identifing the underlying iterator
    unsigned int typeID;
};



}
}

#endif

