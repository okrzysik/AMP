#ifndef included_AMP_structuredMeshIterators
#define included_AMP_structuredMeshIterators

#include "AMP/mesh/MeshIterator.h"
#include "AMP/mesh/structured/BoxMesh.h"
#include "AMP/mesh/structured/structuredMeshElement.h"
#include <memory>

#include <array>


namespace AMP::Mesh {


class structuredMeshIterator final : public MeshIterator
{
public:
    //! Empty MultiVectorIterator constructor
    structuredMeshIterator();

    //! Range base constructor
    structuredMeshIterator( const BoxMesh::MeshElementIndex &first,
                            const BoxMesh::MeshElementIndex &last,
                            const AMP::Mesh::BoxMesh *mesh,
                            size_t pos = 0 );

    //! Element list constructor
    structuredMeshIterator( std::shared_ptr<const std::vector<BoxMesh::MeshElementIndex>> elements,
                            const AMP::Mesh::BoxMesh *mesh,
                            size_t pos = 0 );

    //! Deconstructor
    virtual ~structuredMeshIterator();

    //! Move constructor
    structuredMeshIterator( structuredMeshIterator && ) = default;

    //! Copy constructor
    structuredMeshIterator( const structuredMeshIterator & );

    //! Move operator
    structuredMeshIterator &operator=( structuredMeshIterator && ) = default;

    //! Assignment operator
    structuredMeshIterator &operator=( const structuredMeshIterator & );

    //! Increment
    MeshIterator &operator++() override;

    //! Increment
    MeshIterator operator++( int ) override;

    //! Decrement
    MeshIterator &operator--() override;

    //! Decrement
    MeshIterator operator--( int ) override;

    // Arithmetic operator+
    MeshIterator operator+( int ) const override;

    // Arithmetic operator+=
    MeshIterator &operator+=( int N ) override;

    //! Check if two iterators are equal
    bool operator==( const MeshIterator &rhs ) const override;

    //! Check if two iterators are not equal
    bool operator!=( const MeshIterator &rhs ) const override;

    //! Return an iterator to the begining
    MeshIterator begin() const override;

    //! Return an iterator to the begining
    MeshIterator end() const override;

    using MeshIterator::operator+;
    using MeshIterator::operator+=;

protected:
    //! Clone the iterator
    MeshIterator *clone() const override;

    // Get the index given the position
    inline BoxMesh::MeshElementIndex getIndex( int pos ) const;

    // Get the elements in the iterator
    std::shared_ptr<const std::vector<BoxMesh::MeshElementIndex>> getElements() const;

    friend class AMP::Mesh::BoxMesh;

private:
    // Data members
    bool d_checkBoundary;
    std::array<bool, 3> d_isPeriodic;
    std::array<int, 3> d_globalSize;
    BoxMesh::MeshElementIndex d_first;
    BoxMesh::MeshElementIndex d_last;
    std::shared_ptr<const std::vector<BoxMesh::MeshElementIndex>> d_elements;
    const AMP::Mesh::BoxMesh *d_mesh;
    mutable structuredMeshElement d_cur_element;

private:
    static constexpr uint32_t getTypeID()
    {
        return AMP::Utilities::hash_char( "structuredMeshIterator" );
    }
};


} // namespace AMP::Mesh

#endif
