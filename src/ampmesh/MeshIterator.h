#ifndef included_AMP_MeshIterators
#define included_AMP_MeshIterators

#include <iterator>
#include "utils/shared_ptr.h"
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
    typedef AMP::shared_ptr<MeshIterator>  shared_ptr;

    //! Empty MeshIterator constructor
    MeshIterator();

    //! Copy constructor
    MeshIterator(const MeshIterator&);

    //! Assignment operator
    MeshIterator& operator=(const MeshIterator&);

    //! Deconstructor
    virtual ~MeshIterator ();

    /** 
     * \brief Pre-Increment
     * \details  Pre-Increment the mesh iterator and return the reference to the iterator.
     *   This should be the fastest way to increment the iterator.
     */
    virtual MeshIterator& operator++();
    
    /** 
     * \brief Post-Increment
     * \details  Post-Increment the mesh iterator and return a reference to a temporary iterator.
     *   This should be avoided and pre-increment used whenever possible.
     */
    virtual MeshIterator operator++(int);

    /** 
     * \brief Pre-Decrement
     * \details  Pre-Decrement the mesh iterator and return the reference to the iterator.
     *   This should be the fastest way to decrement the iterator.
     *   Note: not all iterators support decrementing the iterator (libmesh).
     */
    virtual MeshIterator& operator--();
    
    /** 
     * \brief Post-Decrement
     * \details  Post-Decrement the mesh iterator and return a reference to a temporary iterator.
     *   This should be avoided and pre-decrement used whenever possible.
     *   Note: not all iterators support decrementing the iterator (libmesh).
     */
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

    /** 
     * \brief Arithmetic operator+
     * \details  Random access increment to advance the iterator by N.
     *   Note: not all iterators support random access (libmesh).
     *   In this case, the pre-increment will be used instead and may reduce performance.
     *   Note: the default behavior of all random access iterators will be to call this function
     *     so derived classes only need to impliment this function for improved performance.
     * \param N  Number to increment by (may be negitive)
     */
    virtual MeshIterator operator+(int N) const;

    /** 
     * \brief Arithmetic operator+
     * \details  Random access increment to advance the iterator by the given iterator.
     *   Note: not all iterators support random access (libmesh).
     *   In this case, the pre-increment will be used instead and may reduce performance.
     * \param it  Iterator to add
     */
    virtual MeshIterator operator+(const MeshIterator& it) const;

    /** 
     * \brief Arithmetic operator-
     * \details  Random access decrement to reverse the iterator by N.
     *   Note: not all iterators support random access (libmesh).
     *   In this case, the pre-decrement will be used instead and may reduce performance.
     * \param N  Number to decrement by (may be negitive)
     */
    virtual MeshIterator operator-(int N) const;

    /** 
     * \brief Arithmetic operator-
     * \details  Random access decrement to reverse the iterator by the given iterator.
     *   Note: not all iterators support random access (libmesh).
     *   In this case, the pre-decrement will be used instead and may reduce performance.
     * \param it  Iterator to subtract
     */
    virtual MeshIterator operator-(const MeshIterator& it) const;

    /** 
     * \brief Arithmetic operator+=
     * \details  Random access increment to advance the iterator by N.
     *   Note: not all iterators support random access (libmesh).
     *   In this case, the pre-increment will be used instead and may reduce performance.
     *   Note: the default behavior of all random access assignement iterators will be to call this function
     *     so derived classes only need to impliment this function for improved performance.
     * \param N  Number to increment by (may be negitive)
     */
    virtual MeshIterator& operator+=(int N);

    /** 
     * \brief Arithmetic operator+=
     * \details  Random access increment to advance the iterator by the given iterator.
     *   Note: not all iterators support random access (libmesh).
     *   In this case, the pre-increment will be used instead and may reduce performance.
     * \param it  Iterator to add
     */
    virtual MeshIterator& operator+=(const MeshIterator& it);

    /** 
     * \brief Arithmetic operator-=
     * \details  Random access decrement to reverse the iterator by N.
     *   Note: not all iterators support random access (libmesh).
     *   In this case, the pre-decrement will be used instead and may reduce performance.
     * \param N  Number to decrement by (may be negitive)
     */
    virtual MeshIterator& operator-=(int N);

    /** 
     * \brief Arithmetic operator-=
     * \details  Random access decrement to reverse the iterator by the given iterator.
     *   Note: not all iterators support random access (libmesh).
     *   In this case, the pre-decrement will be used instead and may reduce performance.
     * \param it  Iterator to subtract
     */
    virtual MeshIterator& operator-=(const MeshIterator& it);

    //! Operator <
    virtual bool operator<(const MeshIterator&);

    //! Operator <=
    virtual bool operator<=(const MeshIterator& it);

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

