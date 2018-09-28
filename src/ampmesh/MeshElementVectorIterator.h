#ifndef included_AMP_MultiVectorIterator
#define included_AMP_MultiVectorIterator

#include "AMP/ampmesh/MeshIterator.h"
#include "AMP/utils/shared_ptr.h"
#include <iterator>

namespace AMP {
namespace Mesh {


/**
 * \class MultiVectorIterator
 * \brief A class used to iterate over a set of mesh elements.
 * \details  This class provides routines for iterating over a set
 * of mesh elments that are in a std::vector.
 */
class MultiVectorIterator : public MeshIterator
{
public:
    //! Empty MultiVectorIterator constructor
    MultiVectorIterator();

    //! Default MultiVectorIterator constructor
    MultiVectorIterator( AMP::shared_ptr<std::vector<MeshElement>> elements, size_t pos = 0 );

    /** MultiVectorIterator constructor
     *  Note that this version of the constructor will create a copy of the elements
     */
    MultiVectorIterator( const std::vector<MeshElement> &elements, size_t pos = 0 );

    //! Deconstructor
    virtual ~MultiVectorIterator();

    //! Copy constructor
    MultiVectorIterator( const MultiVectorIterator & );

    //! Assignment operator
    MultiVectorIterator &operator=( const MultiVectorIterator & );

    //! Increment
    MeshIterator &operator++() override;

    //! Increment
    MeshIterator operator++( int ) override;

    //! Decrement
    MeshIterator &operator--() override;

    //! Decrement
    MeshIterator operator--( int ) override;

    // Arithmetic operator+
    virtual MeshIterator operator+( int ) const override;

    // Arithmetic operator+=
    virtual MeshIterator &operator+=( int N ) override;

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
    virtual MeshIterator *clone() const override;

    // A pointer to a std::vector containing the desired mesh elements
    AMP::shared_ptr<std::vector<MeshElement>> d_elements;

private:
    static constexpr uint32_t getTypeID()
    {
        return AMP::Utilities::hash_char( "MultiVectorIterator" );
    }
};
} // namespace Mesh
} // namespace AMP

#endif
