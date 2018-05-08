#ifndef included_AMP_libMeshIterators
#define included_AMP_libMeshIterators

#include "AMP/ampmesh/MeshIterator.h"
#include "AMP/ampmesh/libmesh/libMesh.h"


namespace AMP {
namespace Mesh {


class libMeshIterator : public MeshIterator
{
public:
    //! Empty MeshIterator constructor
    libMeshIterator();

    //! Deconstructor
    virtual ~libMeshIterator();

    //! Copy constructor
    libMeshIterator( const libMeshIterator & );

    //! Assignment operator
    libMeshIterator &operator=( const libMeshIterator & );

    // Increment
    virtual MeshIterator &operator++() override;

    // Increment
    virtual MeshIterator operator++( int ) override;

    // Decrement
    virtual MeshIterator &operator--() override;

    // Decrement
    virtual MeshIterator operator--( int ) override;

    // Arithmetic operator+
    virtual MeshIterator operator+( int ) const override;

    // Arithmetic operator+=
    virtual MeshIterator &operator+=( int N ) override;

    // Check if two iterators are equal
    virtual bool operator==( const MeshIterator &rhs ) const override;

    // Check if two iterators are not equal
    virtual bool operator!=( const MeshIterator &rhs ) const override;

    // Return an iterator to the begining
    virtual MeshIterator begin() const override;

    // Return an iterator to the begining
    virtual MeshIterator end() const override;

    using MeshIterator::operator+;
    using MeshIterator::operator+=;

protected:
    /** Default constructor
     * \param type      Entity type:  0: node, 1: element
     * \param mesh      Pointer to the libMesh mesh
     * \param gcw       gcw to use
     * \param begin     Pointer to iterator with the begining position
     * \param end       Pointer to iterator with the end position
     * \param pos       Pointer to iterator with the current position
     * \param size      Number of elements in the iterator (-1: unknown)
     * \param pos2      Index of the current position in the iterator (-1: unknown)
     */
    libMeshIterator( int type,
                     const AMP::Mesh::libMesh *mesh,
                     int gcw,
                     void *begin,
                     void *end,
                     void *pos,
                     int size = -1,
                     int pos2 = -1 );

    //! Clone the iterator
    virtual MeshIterator *clone() const override;

    friend class AMP::Mesh::libMesh;

private:
    // Data members
    int d_gcw;
    int d_dim;
    int d_type;
    int d_rank;
    void *d_begin2;
    void *d_end2;
    void *d_pos2;
    MeshID d_meshID;
    const AMP::Mesh::libMesh *d_mesh;
    MeshElement d_cur_element;

    void setCurrentElement();
};
} // namespace Mesh
} // namespace AMP

#endif
