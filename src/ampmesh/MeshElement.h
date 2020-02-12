#ifndef included_AMP_MeshElement
#define included_AMP_MeshElement

#include "AMP/ampmesh/MeshID.h"
#include "AMP/ampmesh/MeshPoint.h"
#include "AMP/utils/Utilities.h"
#include <memory>
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
    typedef std::shared_ptr<MeshElement> shared_ptr;

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
    inline GeomType elementType() const { return globalID().type(); }


    // Overload operators
    inline bool operator==( const MeshElement &rhs ) const { return globalID() == rhs.globalID(); }
    inline bool operator!=( const MeshElement &rhs ) const { return globalID() != rhs.globalID(); }
    inline bool operator<( const MeshElement &rhs ) const { return globalID() < rhs.globalID(); }
    inline bool operator>( const MeshElement &rhs ) const { return globalID() > rhs.globalID(); }
    inline bool operator<=( const MeshElement &rhs ) const { return globalID() <= rhs.globalID(); }
    inline bool operator>=( const MeshElement &rhs ) const { return globalID() >= rhs.globalID(); }
    inline bool operator==( const MeshElementID &rhs ) const { return globalID() == rhs; }
    inline bool operator!=( const MeshElementID &rhs ) const { return globalID() != rhs; }
    inline bool operator<( const MeshElementID &rhs ) const { return globalID() < rhs; }
    inline bool operator>( const MeshElementID &rhs ) const { return globalID() > rhs; }
    inline bool operator<=( const MeshElementID &rhs ) const { return globalID() <= rhs; }
    inline bool operator>=( const MeshElementID &rhs ) const { return globalID() >= rhs; }

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
    inline double coord( int i ) const { return coord()[i]; }


    /**
     * \brief     Return true if the element contains the point
     * \details   This function checks if the given point is inside or
     *   within TOL of the given element.  If the current element is a vertex,
     *   this function checks if the point is with TOL of the vertex.
     * \param pos   The coordinates of the point to check.
     * \param TOL   The tolerance to use for the computation.
     */
    inline bool containsPoint( const std::vector<double> &pos, double TOL = 1e-12 ) const;


public: // Virtual functions
    //! Return the unique global ID of the element
    virtual MeshElementID globalID() const;

    //! Return the owner rank according to AMP_COMM_WORLD
    virtual unsigned int globalOwnerRank() const;

    //! Return the element class
    virtual std::string elementClass() const;

    //! Return the coordinates of the vertex (only applies to verticies)
    virtual Point coord() const;

    /**
     * \brief     Return the centroid of the element
     * \details   This function returns the centroid of the element.  The
     *   centroid is defined as the average of the coordinates of the verticies.
     *   The centroid of a vertex is the vertex and will return the same result as coord().
     */
    virtual Point centroid() const;

    //! Return the volume of the current element (does not apply to verticies)
    virtual double volume() const;

    //! Return the normal to the current element (does not apply to all elements)
    virtual Point norm() const;

    /**
     * \brief     Return true if the element contains the point
     * \details   This function checks if the given point is inside or
     *   within TOL of the given element.  If the current element is a vertex,
     *   this function checks if the point is with TOL of the vertex.
     * \param pos   The coordinates of the point to check.
     * \param TOL   The tolerance to use for the computation.
     */
    virtual bool containsPoint( const Point &pos, double TOL = 1e-12 ) const;

    /**
     * \brief    Calculate the nearest point on the element
     * \details  This function computes nearest point on/in the element to the given point
     * \param[in] pos   Current position of the point
     */
    virtual Point nearest( const Point &pos ) const;

    /**
     * \brief    Calculate the distance to the element given a ray
     * \details  This function computes the distance to the element given a ray.
     *     If the ray will never intersect the element, this distance is inf.
     * \param[in] pos   Current position of ray
     * \param[in] dir   Direction of ray (should be normalized for most uses)
     * @return          Returns the distance to the element surface
     */
    virtual double distance( const Point &pos, const Point &dir ) const;

    //! Check if the element is on the surface
    virtual bool isOnSurface() const;

    /**
     * \brief     Check if the current element is on the given boundary
     * \details   Check if the current element is on the boundary specified by the given id
     * \param id  The boundary id to check
     */
    virtual bool isOnBoundary( int id ) const;

    /**
     * \brief     Check if the current element is in the given block
     * \details   Check if the current element is in the block specified by the given id
     * \param id  The block id to check
     */
    virtual bool isInBlock( int id ) const;


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


protected:
    // Unique (per class) ID for identifing the underlying iterator
    unsigned int typeID;

    // A pointer to the derived class
    MeshElement *element;

    // Clone the element
    virtual inline MeshElement *clone() const;

private:
    static constexpr uint32_t getTypeID() { return AMP::Utilities::hash_char( "MeshElement" ); }
};


} // namespace Mesh
} // namespace AMP

#include "AMP/ampmesh/MeshElement.inline.h"

#endif
