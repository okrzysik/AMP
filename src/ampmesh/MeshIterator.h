#ifndef included_AMP_MeshIterators
#define included_AMP_MeshIterators

#include <iterator>
#include <boost/shared_ptr.hpp>
#include "MeshElement.h"

namespace AMP { 
namespace Mesh {


class MeshIterator
{

public:
    /**
     *\typedef shared_ptr
     *\brief  Name for the shared pointer.
     *\details  Use this typedef for a reference counted pointer to a mesh manager object.
     */
    typedef boost::shared_ptr<MeshIterator>  shared_ptr;

    //! Empty MeshIterator constructor
    MeshIterator();

    //! Copy constructor
    MeshIterator(const MeshIterator&);

    //! Assignment operator
    MeshIterator& operator=(const MeshIterator&);

    //! Deconstructor
    virtual ~MeshIterator ();

    //! Pre-Increment
    virtual MeshIterator& operator++();
    
    //! Post-Increment (creates a temporary object)
    virtual MeshIterator operator++(int);

    //! Pre-Decrement
    virtual MeshIterator& operator--();
    
    //! Post-Decrement (creates a temporary object)
    virtual MeshIterator operator--(int);

    //! Check if two iterators are equal
    virtual bool operator==(const MeshIterator& rhs);

    //! Check if two iterators are not equal
    virtual bool operator!=(const MeshIterator& rhs);
    
    //! Dereference the iterator
    virtual MeshElement &operator*(void);

    //! Dereference the iterator
    virtual MeshElement *operator->(void);

    //! Return an iterator to the begining
    virtual MeshIterator begin();

    //! Return an iterator to the begining
    virtual MeshIterator end();

protected:
    // A pointer to the derived class
    MeshIterator *iterator;
    // Clone the iterator
    virtual MeshIterator* clone() const;
    // Unique (per class) ID for identifing the underlying iterator
    unsigned int typeID;
};


}
}

#endif

