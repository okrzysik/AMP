#ifndef included_AMP_TriangleMeshElement
#define included_AMP_TriangleMeshElement


#include "AMP/ampmesh/MeshElement.h"
#include "AMP/utils/shared_ptr.h"
#include <vector>


namespace AMP {
namespace Mesh {


template<size_t NG, size_t NP>
class TriangleMesh;
template<size_t NG, size_t NP>
class TriangleMeshIterator;


/**
 * \class TriangleMeshElement
 * \brief A derived class used to define a mesh element
 * \details  This class provides routines for accessing and using a mesh element.
 * A mesh element can be thought of as the smallest unit of a mesh.  It is of a type
 * of GeomType.  This class is derived to store a TriangleMesh element.
 */
template<size_t NG, size_t NP>
class TriangleMeshElement final : public MeshElement
{
public:
    //! Empty constructor for a MeshElement
    TriangleMeshElement();

    //! Copy constructor
    TriangleMeshElement( const TriangleMeshElement & );

    //! Assignment operator
    TriangleMeshElement &operator=( const TriangleMeshElement & );

    //! Move operator
    TriangleMeshElement( TriangleMeshElement && );

    //! Move assignment operator
    TriangleMeshElement &operator=( TriangleMeshElement && );

    //! De-constructor for a MeshElement
    virtual ~TriangleMeshElement() = default;

    //! Return the element class
    virtual inline std::string elementClass() const override { return "TriangleMeshElement"; }

    //! Return the elements composing the current element
    virtual void getElements( const GeomType type,
                              std::vector<MeshElement> &elements ) const override;

    //! Return the IDs of the elements composing the current element
    virtual void getElementsID( const GeomType type,
                                std::vector<MeshElementID> &ID ) const override;

    //! Return the elements neighboring the current element
    virtual void getNeighbors( std::vector<MeshElement::shared_ptr> &neighbors ) const override;

    //! Return the volume of the current element (does not apply to verticies)
    virtual double volume() const override;

    //! Return the coordinates of all verticies composing the element
    virtual Point coord() const override;

    /**
     * \brief     Return the centroid of the element
     * \details   This function returns the centroid of the element.  The
     *   centroid is defined as the average of the coordinates of the verticies.
     *   The centroid of a vertex is the vertex and will return the same result as coord().
     */
    virtual Point centroid() const override;

    /**
     * \brief     Return true if the element contains the point
     * \details   This function checks if the given point is inside or
     *   within TOL of the given element.  If the current element is a vertex,
     *   this function checks if the point is with TOL of the vertex.
     * \param pos   The coordinates of the point to check.
     * \param TOL   The tolerance to use for the computation.
     */
    virtual bool containsPoint( const Point &pos, double TOL = 1e-12 ) const override;

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
    // Default constructors
    TriangleMeshElement( const MeshElementID &id, const TriangleMesh<NG, NP> *mesh );

    // Reset the element data
    inline void resetElemId( const ElementID &id ) { d_globalID.resetElemID( id ); }

    //! Clone the iterator
    virtual MeshElement *clone() const override;

    // The pointer to the current mesh
    const TriangleMesh<NG, NP> *d_mesh;

    friend class AMP::Mesh::TriangleMesh<NG, NP>;
    friend class AMP::Mesh::TriangleMeshIterator<NG, NP>;

private:
    static constexpr uint32_t getTypeID();
};
} // namespace Mesh
} // namespace AMP

#endif
