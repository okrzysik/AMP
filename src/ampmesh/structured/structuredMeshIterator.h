#ifndef included_AMP_structuredMeshIterators
#define included_AMP_structuredMeshIterators

#include "ampmesh/MeshIterator.h"
#include "ampmesh/structured/BoxMesh.h"
#include "ampmesh/structured/structuredMeshElement.h"
#include "utils/shared_ptr.h"

namespace AMP {
namespace Mesh {


class structuredMeshIterator : public MeshIterator {
public:
    //! Empty MultiVectorIterator constructor
    structuredMeshIterator();

    //! Default MultiVectorIterator constructor
    structuredMeshIterator( AMP::shared_ptr<std::vector<BoxMesh::MeshElementIndex>> elements,
                            const AMP::Mesh::BoxMesh *mesh,
                            size_t pos = 0 );

    //! Deconstructor
    virtual ~structuredMeshIterator();

    //! Copy constructor
    structuredMeshIterator( const structuredMeshIterator & );

    //! Assignment operator
    structuredMeshIterator &operator=( const structuredMeshIterator & );

    //! Increment
    MeshIterator &operator++();

    //! Increment
    MeshIterator operator++( int );

    //! Decrement
    MeshIterator &operator--();

    //! Decrement
    MeshIterator operator--( int );

    // Arithmetic operator+
    virtual MeshIterator operator+( int ) const;

    // Arithmetic operator+=
    virtual MeshIterator &operator+=( int N );

    //! Check if two iterators are equal
    bool operator==( const MeshIterator &rhs ) const;

    //! Check if two iterators are not equal
    bool operator!=( const MeshIterator &rhs ) const;

    //! Dereference the iterator
    MeshElement &operator*( void );

    //! Dereference the iterator
    MeshElement *operator->( void );

    //! Return an iterator to the begining
    MeshIterator begin() const;

    //! Return an iterator to the begining
    MeshIterator end() const;

    //! Return the number of elements in the iterator
    virtual size_t size() const;

    //! Return the current position (from the beginning) in the iterator
    virtual size_t position() const;

    using MeshIterator::operator+;
    using MeshIterator::operator+=;

protected:
    //! Clone the iterator
    virtual MeshIterator *clone() const;

    friend class AMP::Mesh::BoxMesh;

private:
    // Data members
    size_t d_pos;
    AMP::shared_ptr<std::vector<BoxMesh::MeshElementIndex>> d_elements;
    const AMP::Mesh::BoxMesh *d_mesh;
    structuredMeshElement d_cur_element;
};
}
}

#endif
