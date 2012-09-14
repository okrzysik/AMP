#ifndef included_AMP_MeshIterators
#define included_AMP_MeshIterators

#include <iterator>
#include <boost/shared_ptr.hpp>
#include "ampmesh/MeshElement.h"

namespace AMP { 
namespace Mesh {


/**
 * \class MeshIterator
 * \brief A class used to iterate over elements in a Mesh
 *
 * \details  This class provides routines for iterating over a set of elements.
 *   It is inherited from std::iterator.  The base class provides routines for
 *   the random access iterators, but does so using the increment/decrement routines.
 *   Derived classes may (or may not) override these routines for performance optimizations.
 */
class MeshIterator : public std::iterator<std::random_access_iterator_tag,AMP::Mesh::MeshElement>
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
    virtual bool operator==(const MeshIterator& rhs) const;

    //! Check if two iterators are not equal
    virtual bool operator!=(const MeshIterator& rhs) const;
    
    //! Dereference the iterator
    virtual MeshElement &operator*(void);

    //! Dereference the iterator
    virtual MeshElement *operator->(void);

    //! Return an iterator to the begining
    virtual MeshIterator begin() const;

    //! Return an iterator to the begining
    virtual MeshIterator end() const;

    //! Return the number of elements in the iterator
    virtual size_t size() const;

    //! Return the current position (from the beginning) in the iterator
    virtual size_t position() const;

    //! Return the Unique (per class) ID for identifing the underlying iterator
    unsigned int type_id() const {return typeID;}

    //! Arithmetic operator+
    virtual MeshIterator operator+(int) const;

    //! Arithmetic operator+
    virtual MeshIterator operator+(const MeshIterator&) const;

    //! Arithmetic operator-
    virtual MeshIterator operator-(int) const;

    //! Arithmetic operator-
    virtual MeshIterator operator-(const MeshIterator&) const;

    //! Arithmetic operator+=
    virtual MeshIterator& operator+=(int);

    //! Arithmetic operator+=
    virtual MeshIterator& operator+=(const MeshIterator&);

    //! Arithmetic operator-=
    virtual MeshIterator& operator-=(int);

    //! Arithmetic operator-=
    virtual MeshIterator& operator-=(const MeshIterator&);

    //! Operator <
    virtual bool operator<(const MeshIterator&);

    //! Operator <=
    virtual bool operator<=(const MeshIterator&);

    //! Operator >
    virtual bool operator>(const MeshIterator&);

    //! Operator >=
    virtual bool operator>=(const MeshIterator&);

    //! Dereference iterator with offset
    virtual MeshElement& operator[](int);

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

