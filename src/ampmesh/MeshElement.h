#ifndef included_AMP_MeshElement
#define included_AMP_MeshElement

#include <vector>
#include <boost/shared_ptr.hpp>

namespace AMP {
namespace Mesh {


//! Enumeration for basic mesh-based quantities
enum GeomType { Vertex=0, Edge=1, Face=2, Volume=3, null=-1 };

typedef size_t MeshElementID;


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

    //! Return the element type
    virtual GeomType elementType() { return d_elementType; }

    //! Return the unique global ID of the element
    virtual MeshElementID globalID() { return d_globalID; }

    //! Return the elements composing the current element
    virtual std::vector<MeshElement> getElements(GeomType &type);

    //! Return the elements neighboring the current element
    virtual std::vector<MeshElement> getNeighbors();

    //! Return the volume of the current element (does not apply to verticies)
    virtual double volume();

    //! Return the coordinates of all verticies composing the element
    virtual std::vector<double> coord();

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

