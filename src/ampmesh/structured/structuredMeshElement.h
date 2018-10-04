#ifndef included_AMP_structuredMeshElement
#define included_AMP_structuredMeshElement

#include "AMP/ampmesh/structured/BoxMesh.h"
#include <vector>

namespace AMP {
namespace Mesh {


/**
 * \class structuredMeshElement
 * \brief A derived class used to define a mesh element
 * \details  This class provides routines for accessing and using a mesh element.
 * A mesh element can be thought of as the smallest unit of a mesh.  It is of a type
 * of GeomType.  This class is derived to store a libMesh element.
 */
class structuredMeshElement : public MeshElement
{
public:
    //! Empty constructor for a MeshElement
    structuredMeshElement();

    //! Copy constructor
    structuredMeshElement( const structuredMeshElement & );

    //! Assignment operator
    structuredMeshElement &operator=( const structuredMeshElement & );

    //! De-constructor for a MeshElement
    virtual ~structuredMeshElement();

    //! Return the element class
    virtual inline std::string elementClass() const override { return "structuredMeshElement"; }

    //! Return the elements composing the current element
    virtual void getElements( const GeomType type,
                              std::vector<MeshElement> &elements ) const override;

    //! Return the IDs of the elements composing the current element
    virtual void getElementsID( const GeomType type,
                                std::vector<MeshElementID> &ID ) const override;

    //! Return the elements neighboring the current element
    virtual void getNeighbors( std::vector<MeshElement::shared_ptr> &neighbors ) const override;

    /**
     * \brief     Return the centroid of the element
     * \details   This function returns the centroid of the element.  The
     *   centroid is defined as the average of the coordinates of the verticies.
     *   The centroid of a vertex is the vertex and will return the same result as coord().
     */
    virtual Point centroid() const override;

    //! Return the volume of the current element (does not apply to verticies)
    virtual double volume() const override;

    //! Return the coordinates of the vertex (only applies to verticies)
    virtual Point coord() const override final
    {
        Point x;
        x.setNdim( d_physicalDim );
        d_mesh->coord( d_index, x.data() );
        return x;
    }

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


    /**
     * \brief     Get the parents of the given element
     * \details   This function will get the parent elements of the current element
     * \param type  The desired type of the parents to get
     */
    virtual std::vector<MeshElement> getParents( GeomType type ) const;


    //! Return the owner rank according to AMP_COMM_WORLD
    virtual unsigned int globalOwnerRank() const override;


    //! Return the index of the element
    BoxMesh::MeshElementIndex getIndex() const { return d_index; }

public:
    /** Default constructor
     * \param index     Index for the current elements
     * \param mesh      Underlying mesh
     */
    structuredMeshElement( const BoxMesh::MeshElementIndex &index, const AMP::Mesh::BoxMesh *mesh );

    //! Reset the internal data to an empty element
    void reset();

    //! Reset the internal data to the given element
    void reset( const BoxMesh::MeshElementIndex &index, const AMP::Mesh::BoxMesh *mesh );


protected:
    // Clone the iterator
    virtual MeshElement *clone() const override;

    // Internal data
    GeomType d_meshType;               // Mesh logical dimension
    unsigned char d_physicalDim;       // Mesh physical dimension
    BoxMesh::MeshElementIndex d_index; // Index of element
    const AMP::Mesh::BoxMesh *d_mesh;  // Pointer to mesh

    // Helper functions
    void getNeighborIndex( int &N, BoxMesh::MeshElementIndex *index ) const;
    void getElementIndex( const GeomType type, int &N, BoxMesh::MeshElementIndex *index ) const;

    // Reset just the mesh element
    void reset( const BoxMesh::MeshElementIndex &index )
    {
        d_index    = index;
        d_globalID = d_mesh->convert( index );
    }

    friend class AMP::Mesh::BoxMesh;
    friend class AMP::Mesh::structuredMeshIterator;

private:
    static constexpr uint32_t getTypeID()
    {
        return AMP::Utilities::hash_char( "structuredMeshElement" );
    }
};


} // namespace Mesh
} // namespace AMP


#endif
