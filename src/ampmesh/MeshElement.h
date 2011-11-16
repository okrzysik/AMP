#ifndef included_AMP_MeshElement
#define included_AMP_MeshElement

#include <vector>
#include <boost/shared_ptr.hpp>

namespace AMP {
namespace Mesh {


//! Enumeration for basic mesh-based quantities
enum GeomType { Vertex=0, Edge=1, Face=2, Volume=3, null=-1 };


/**
 * \struct MeshElementID
 * \brief  A structure used to identify the mesh element
 * \details  This structure provides a unique id that can be used to identify a mesh element.
 */
struct MeshElementID{
    bool is_local;              //!<  Is the current element local
    GeomType type;              //!<  The geometric type of the element
    unsigned int local_id;      //!<  The local id of the element on the owner processor and mesh
    unsigned int owner_rank;    //!<  The rank of the owner proccessor (as defined from rank of the comm on the owner mesh)
    size_t meshID;              //!<  The ID of the mesh that owns the rank (not implimented yet)
    // Constructors used to initialize key values
	MeshElementID() {
        is_local = false;
        type = null;
        local_id = static_cast<unsigned int>(-1);
        owner_rank = static_cast<unsigned int>(-1);
        meshID = static_cast<size_t>(-1);
    }
	MeshElementID(bool isLocal, GeomType type_id, unsigned int local_ID, unsigned int owner_rank_id, size_t mesh_ID) {
        is_local = isLocal;
        type = type_id;
        local_id = local_ID;
        owner_rank = owner_rank_id;
        meshID = mesh_ID;
    }
    // Overload key operators
    inline bool operator== (const MeshElementID& rhs ) const {
        return type==rhs.type && local_id==rhs.local_id && owner_rank==rhs.owner_rank && meshID==rhs.meshID;
    }
    inline bool operator!= (const MeshElementID& rhs ) const {
        return type!=rhs.type || local_id!=rhs.local_id || owner_rank!=rhs.owner_rank || meshID!=rhs.meshID;
    }
    inline bool operator>= (const MeshElementID& rhs ) const {
        // Sort by meshID first
        if ( meshID < rhs.meshID )
            return false;
        else if ( meshID > rhs.meshID )
            return true;
        // Sort by processor id next
        if ( owner_rank < rhs.owner_rank )
            return false;
        else if ( owner_rank > rhs.owner_rank )
            return true;
        // Sort by type next
        if ( type < rhs.type )
            return false;
        else if ( type > rhs.type )
            return true;
        // Finally check the local id
        if ( local_id < rhs.local_id )
            return false;
        return true;
    }
    inline bool operator> (const MeshElementID& rhs ) const {
        if ( this->operator>=(rhs) && this->operator!=(rhs) )
            return true;
        return false;
    }
    inline bool operator< (const MeshElementID& rhs ) const {
        return !(this->operator>=(rhs));
    }
    inline bool operator<= (const MeshElementID& rhs ) const {
        return !(this->operator>=(rhs)) || this->operator==(rhs);
    }
};


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
    virtual GeomType elementType() const { return d_elementType; }

    //! Return the unique global ID of the element
    virtual MeshElementID globalID() const { return d_globalID; }

    //! Return the elements composing the current element
    virtual std::vector<MeshElement> getElements(const GeomType type) const;

    //! Return the elements neighboring the current element
    virtual std::vector<MeshElement> getNeighbors() const;

    //! Return the volume of the current element (does not apply to verticies)
    virtual double volume() const;

    //! Return the coordinates of all verticies composing the element
    virtual std::vector<double> coord() const;

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

