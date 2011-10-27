#ifndef included_AMP_MeshElement
#define included_AMP_MeshElement

#include <vector>

namespace AMP {
namespace Mesh {


//! Enumeration for basic mesh-based quantities
enum GeomType { Vertex=0, Edge=1, Face=2, Volume=3 };


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
    //! Return the element type
    virtual GeomType getElementType() { return elementType; }

    //! Return the unique global ID of the element
    virtual ID getGlobalId() { return globalID; }

    //! Return the elements composing the current element
    virtual vector<MeshElement> getElements(GeomType &type);

    //! Return the elements neighboring the current element
    virtual vector<MeshElement> getNeighbors();

    //! Return the volume of the current element (does not apply to verticies)
    virtual double getVolume();

    //! Return the coordinates of all verticies composing the element
    virtual void getCoord(std::vector<double> &coords);

protected:

    //!  Empty constructor for a MeshElement
    MeshElement ( );

    //! The MeshElement type
    const GeomType elementType;

    //! The unique global id of the element
    const MeshID globalID;
};



}
}

#endif

