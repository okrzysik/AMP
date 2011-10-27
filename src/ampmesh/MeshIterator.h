#ifndef included_AMP_MeshIterators
#define included_AMP_MeshIterators

#include <iterator>
#include "MeshElement.h"

namespace AMP { 
namespace Mesh {


class MeshIterator: public iterator
{

public:

    /**
     * \brief MeshIterator constructor
     * \details  This constructor will construct a new MeshIterator from an existing mesh iterator.
     * \param mit Existing MeshIterator
     */
    MeshIterator(const MeshIterator& mit);

    //! Increment
    virtual MeshIterator& operator++();
    
    //! Increment
    virtual MeshIterator operator++(int);

    //! Decrement
    virtual MeshIterator& operator--();
    
    //! Decrement
    virtual MeshIterator operator--(int);

    //! Check if two iterators are equal
    virtual bool operator==(const MeshIterator& rhs);

    //! Check if two iterators are not equal
    virtual bool operator!=(const MeshIterator& rhs);
    
    //! Dereference the iterator
    virtual MeshElement& operator*();
};

protected:
    //! Empty MeshIterator constructor
    MeshIterator();


}
}

#endif

