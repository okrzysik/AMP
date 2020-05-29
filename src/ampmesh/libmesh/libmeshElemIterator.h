#ifndef included_AMP_libmeshElemIterator
#define included_AMP_libmeshElemIterator

#include "AMP/ampmesh/MeshIterator.h"
#include "AMP/ampmesh/libmesh/libmeshMesh.h"

// libMesh includes
#include "libmesh/elem.h"


namespace AMP {
namespace Mesh {


class libmeshElemIterator : public MeshIterator
{
public:
    //! Empty MeshIterator constructor
    libmeshElemIterator() = delete;

    //! Deconstructor
    virtual ~libmeshElemIterator() = default;

    //! Copy constructor
    libmeshElemIterator( const libmeshElemIterator & );

    //! Assignment operator
    libmeshElemIterator &operator=( const libmeshElemIterator & );

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
     * \param mesh      Pointer to the libMesh mesh
     * \param gcw       gcw to use
     * \param begin     Pointer to iterator with the begining position
     * \param end       Pointer to iterator with the end position
     * \param pos       Pointer to iterator with the current position
     * \param size      Number of elements in the iterator (-1: unknown)
     * \param pos2      Index of the current position in the iterator (-1: unknown)
     */
    libmeshElemIterator( const AMP::Mesh::libmeshMesh *mesh,
                         int gcw,
                         const libMesh::Mesh::element_iterator &begin,
                         const libMesh::Mesh::element_iterator &end,
                         const libMesh::Mesh::element_iterator &pos,
                         int size = -1,
                         int pos2 = -1 );

    //! Clone the iterator
    virtual MeshIterator *clone() const override;

    friend class AMP::Mesh::libmeshMesh;

private:
    // Data members
    int d_gcw;
    int d_dim;
    int d_rank;
    libMesh::Mesh::element_iterator d_begin2;
    libMesh::Mesh::element_iterator d_end2;
    libMesh::Mesh::element_iterator d_pos2;
    MeshID d_meshID;
    const AMP::Mesh::libmeshMesh *d_mesh;
    MeshElement d_cur_element;

    void setCurrentElement();

private:
    static constexpr uint32_t getTypeID()
    {
        return AMP::Utilities::hash_char( "libmeshElemIterator" );
    }
};
} // namespace Mesh
} // namespace AMP

#endif
