#ifndef included_AMP_TriangleMeshIterators
#define included_AMP_TriangleMeshIterators


#include "AMP/ampmesh/MeshIterator.h"
#include "AMP/ampmesh/triangle/TriangleMeshElement.h"


namespace AMP {
namespace Mesh {


template<size_t NG, size_t NP>
class TriangleMesh;


template<uint8_t NG, uint8_t NP, uint8_t TYPE>
class TriangleMeshIterator final : public MeshIterator
{
public:
    //! Empty MeshIterator constructor
    TriangleMeshIterator();

    //! Deconstructor
    virtual ~TriangleMeshIterator() = default;

    //! Copy constructor
    TriangleMeshIterator( const TriangleMeshIterator & );

    //! Assignment operator
    TriangleMeshIterator &operator=( const TriangleMeshIterator & );

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
     * \param list      List of elements
     * \param pos       Pointer to iterator with the current position
     */
    explicit TriangleMeshIterator( const AMP::Mesh::TriangleMesh<NG, NP> *mesh,
                                   std::shared_ptr<const std::vector<ElementID>> list,
                                   size_t pos = 0 );

    //! Clone the iterator
    virtual MeshIterator *clone() const override;

    friend class AMP::Mesh::TriangleMesh<NG, NP>;

private:
    // Data members
    const AMP::Mesh::TriangleMesh<NG, NP> *d_mesh;
    std::shared_ptr<const std::vector<ElementID>> d_list;
    TriangleMeshElement<NG, NP, TYPE> d_cur_element;

private:
    static constexpr uint32_t getTypeID();
};


} // namespace Mesh
} // namespace AMP

#endif
