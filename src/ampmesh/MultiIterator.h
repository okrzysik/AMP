#ifndef included_AMP_MultiIterator
#define included_AMP_MultiIterator

#include <boost/shared_ptr.hpp>
#include "ampmesh/Mesh.h"
#include "ampmesh/MeshIterator.h"

namespace AMP { 
namespace Mesh {


/**
 * \class MultiIterator
 * \brief A class used to combine multiple iterators
 * \details  This class provides routines for iterating over a mesh.  More
 *  specifically, this class combines multiple iterators into one.  This
 *  is primarily needed for MultiMesh, but may be used for other applicaitons.
 */
class MultiIterator: public MeshIterator {
public:

    //! Empty MultiIterator constructor
    MultiIterator();

    //! Default MultiIterator constructor
    MultiIterator( std::vector<boost::shared_ptr<MeshIterator> > iterators, size_t global_pos=0 );

    //! Deconstructor
    ~MultiIterator ();

    //! Copy constructor
    MultiIterator(const MultiIterator&);

    //! Assignment operator
    MultiIterator& operator=(const MultiIterator&);

    //! Increment
    MeshIterator& operator++();
    
    //! Increment
    MeshIterator operator++(int);

    //! Decrement
    MeshIterator& operator--();
    
    //! Decrement
    MeshIterator operator--(int);

    //! Check if two iterators are equal
    bool operator==(const MeshIterator& rhs);

    //! Check if two iterators are not equal
    bool operator!=(const MeshIterator& rhs);
    
    //! Dereference the iterator
    MeshElement &operator*(void);

    //! Dereference the iterator
    MeshElement *operator->(void);

    //! Return an iterator to the begining
    MeshIterator begin();

    //! Return an iterator to the begining
    MeshIterator end();

    //! Return the number of elements in the iterator
    virtual size_t size() const;

protected:

    //! Clone the iterator
    virtual MeshIterator* clone() const;

private:
    // Data members
    size_t d_localPos, d_globalPos, d_iteratorNum, d_globalSize;
    std::vector<size_t> d_iteratorSize;
    std::vector<boost::shared_ptr<MeshIterator> > d_iterators;
    MeshIterator cur_iterator;
};


}
}

#endif

