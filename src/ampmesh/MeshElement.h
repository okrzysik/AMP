#ifndef included_AMP_MeshElement
#define included_AMP_MeshElement

#include "AMP/ampmesh/MeshID.h"
#include "AMP/utils/shared_ptr.h"
#include <vector>


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
    typedef AMP::shared_ptr<MeshElement> shared_ptr;

    //! Empty constructor for a MeshElement
    inline MeshElement();

    //! Copy constructor
    inline MeshElement( const MeshElement & );

    //! Move constructor
    inline MeshElement( MeshElement && );

    //! Assignment operator
    inline MeshElement &operator=( const MeshElement & );

    //! Move assignment operator
    inline MeshElement &operator=( MeshElement && );

    //! Destructor for a MeshElement
    virtual inline ~MeshElement();

    //! Is the mesh element null
    inline bool isNull() const;


public: // non-virtual functions
    //! Function to get a pointer to the raw mesh element (structuredMeshElement, etc.)
    inline MeshElement *getRawElement();

    //! Function to get a pointer to the raw mesh element (structuredMeshElement, etc.)
    inline const MeshElement *getRawElement() const;

    //! Return the element type
    inline GeomType elementType() const;

    //! Return the unique global ID of the element
    inline const MeshElementID &globalID() const;

    // Overload operators
    inline bool operator==( const MeshElement &rhs ) const { return d_globalID == rhs.d_globalID; }
    inline bool operator!=( const MeshElement &rhs ) const { return d_globalID != rhs.d_globalID; }
    inline bool operator<( const MeshElement &rhs ) const { return d_globalID < rhs.d_globalID; }
    inline bool operator>( const MeshElement &rhs ) const { return d_globalID > rhs.d_globalID; }
    inline bool operator<=( const MeshElement &rhs ) const { return d_globalID <= rhs.d_globalID; }
    inline bool operator>=( const MeshElement &rhs ) const { return d_globalID >= rhs.d_globalID; }
    inline bool operator==( const MeshElementID &rhs ) const { return d_globalID == rhs; }
    inline bool operator!=( const MeshElementID &rhs ) const { return d_globalID != rhs; }
    inline bool operator<( const MeshElementID &rhs ) const { return d_globalID < rhs; }
    inline bool operator>( const MeshElementID &rhs ) const { return d_globalID > rhs; }
    inline bool operator<=( const MeshElementID &rhs ) const { return d_globalID <= rhs; }
    inline bool operator>=( const MeshElementID &rhs ) const { return d_globalID >= rhs; }

    //! Return the coordinates of the vertex (only applies to verticies)
    inline std::vector<double> coord() const;


    //! Return the elements composing the current element
    inline std::vector<MeshElement> getElements( const GeomType type ) const;

    /**
     *  Return the elements neighboring the current element.
     *  One neighbor is returned for each side of the element.
     *  If the side is on the surface, then it's neighbor is null.
     *  For Verticies, a list of all verticies that share an element is returned.
     *  This list is in unsorted order.
     */
    inline std::vector<MeshElement::shared_ptr> getNeighbors() const;

    /**
     * \brief     Return the coordinate of the vertex
     * \details   This function returns the coordinates of the vertex
     *   in the given direction (only applies to verticies).
     * \param i     The direction requested.  Equivalent to coord()[i]
     */
    inline double coord( int i ) const;

    /**
     * \brief     Return the centroid of the element
     * \details   This function returns the centroid of the element.  The
     *   centroid is defined as the average of the coordinates of the verticies.
     *   The centroid of a vertex is the vertex and will return the same result as coord().
     */
    inline std::vector<double> centroid() const;


public: // Virtual functions
    //! Return the owner rank according to AMP_COMM_WORLD
    virtual inline unsigned int globalOwnerRank() const;

    //! Return the element class
    virtual inline std::string elementClass() const;

    //! Return the volume of the current element (does not apply to verticies)
    virtual inline double volume() const;

    /**
     * \brief     Return true if the element contains the point
     * \details   This function checks if the given point is inside or
     *   within TOL of the given element.  If the current element is a vertex,
     *   this function checks if the point is with TOL of the vertex.
     * \param pos   The coordinates of the point to check.
     * \param TOL   The tolerance to use for the computation.
     */
    virtual inline bool containsPoint( const std::vector<double> &pos, double TOL = 1e-12 ) const;

    //! Check if the element is on the surface
    virtual inline bool isOnSurface() const;

    /**
     * \brief     Check if the current element is on the given boundary
     * \details   Check if the current element is on the boundary specified by the given id
     * \param id  The boundary id to check
     */
    virtual inline bool isOnBoundary( int id ) const;

    /**
     * \brief     Check if the current element is in the given block
     * \details   Check if the current element is in the block specified by the given id
     * \param id  The block id to check
     */
    virtual inline bool isInBlock( int id ) const;


public: // Advanced functions
    //! Return the elements composing the current element
    virtual void getElements( const GeomType type, std::vector<MeshElement> &elements ) const;

    //! Return the IDs of the elements composing the current element
    virtual void getElementsID( const GeomType type, std::vector<MeshElementID> &ID ) const;

    /**
     *  Return the elements neighboring the current element.
     *  One neighbor is returned for each side of the element.
     *  If the side is on the surface, then it's neighbor is null.
     *  For Verticies, a list of all verticies that share an element is returned.
     *  This list is in unsorted order.
     */
    virtual void getNeighbors( std::vector<MeshElement::shared_ptr> &neighbors ) const;

    //! Return the coordinates of the vertex (only applies to verticies)
    virtual void coord( size_t &N, double *pos ) const;

    /**
     * \brief     Return the centroid of the element
     * \details   This function returns the centroid of the element.  The
     *   centroid is defined as the average of the coordinates of the verticies.
     *   The centroid of a vertex is the vertex and will return the same result as coord().
     */
    virtual void centroid( size_t &N, double *pos ) const;


protected:
    // Unique (per class) ID for identifing the underlying iterator
    unsigned int typeID;

    // A pointer to the derived class
    MeshElement *element;

    //! The unique global id of the element
    MeshElementID d_globalID;

    // Clone the element
    virtual inline MeshElement *clone() const;
};


} // namespace Mesh
} // namespace AMP

#include "AMP/ampmesh/MeshElement.inline.h"

#endif
