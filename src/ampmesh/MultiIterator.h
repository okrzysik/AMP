#ifndef included_AMP_MultiIterator
#define included_AMP_MultiIterator

#include "ampmesh/MeshIterator.h"
#include "utils/shared_ptr.h"

namespace AMP {
namespace Mesh {


/**
 * \class MultiIterator
 * \brief A class used to combine multiple iterators
 * \details  This class provides routines for iterating over a mesh.  More
 *  specifically, this class combines multiple iterators into one.  This
 *  is primarily needed for MultiMesh, but may be used for other applicaitons.
 */
class MultiIterator : public MeshIterator
{
public:
    //! Empty MultiIterator constructor
    MultiIterator();

    //! Default MultiIterator constructor
    MultiIterator( std::vector<AMP::shared_ptr<MeshIterator>> iterators, size_t global_pos = 0 );

    //! Deconstructor
    virtual ~MultiIterator();

    //! Copy constructor
    MultiIterator( const MultiIterator & );

    //! Assignment operator
    MultiIterator &operator=( const MultiIterator & );

    // Increment
    MeshIterator &operator++();

    // Increment
    MeshIterator operator++( int );

    // Decrement
    MeshIterator &operator--();

    // Decrement
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

private:
    // Data members
    size_t d_localPos, d_globalPos, d_iteratorNum, d_globalSize;
    std::vector<size_t> d_iteratorSize;
    std::vector<AMP::shared_ptr<MeshIterator>> d_iterators;
    MeshIterator cur_iterator;
};
}
}

#endif
