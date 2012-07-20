#ifndef included_AMP_STKMeshIterators
#define included_AMP_STKMeshIterators

#include <boost/shared_ptr.hpp>
#include "ampmesh/STKmesh/STKMesh.h"
#include "ampmesh/MeshIterator.h"

namespace AMP {
namespace Mesh {


class STKMeshIterator: public MeshIterator {
public:

    //! Empty MeshIterator constructor
    STKMeshIterator();

    //! Deconstructor
    ~STKMeshIterator ();

    //! Copy constructor
    STKMeshIterator(const STKMeshIterator&);

    //! Assignment operator
    STKMeshIterator& operator=(const STKMeshIterator&);

    //! Increment
    MeshIterator& operator++();

    //! Increment
    MeshIterator operator++(int);

    //! Decrement
    MeshIterator& operator--();

    //! Decrement
    MeshIterator operator--(int);

    //! Check if two iterators are equal
    bool operator==(const MeshIterator& rhs) const;

    //! Check if two iterators are not equal
    bool operator!=(const MeshIterator& rhs) const;

    //! Dereference the iterator
    MeshElement &operator*(void);

    //! Dereference the iterator
    MeshElement *operator->(void);

    //! Return an iterator to the begining
    MeshIterator begin() const; 

    //! Return an iterator to the begining
    MeshIterator end() const ;

    //! Return the number of elements in the iterator
    virtual size_t size() const;

    //! Return the current position (from the beginning) in the iterator
    virtual size_t position() const;

protected:
    /** Default constructor
     * \param type      Entity type:  0: node, 1: element
     * \param mesh      Pointer to the STKMesh mesh
     * \param gcw       gcw to use
     * \param begin     Pointer to iterator with the begining position
     * \param end       Pointer to iterator with the end position
     * \param pos       Pointer to iterator with the current position
     * \param size      Number of elements in the iterator (-1: unknown)
     * \param pos2      Index of the current position in the iterator (-1: unknown)
     */
    STKMeshIterator(const AMP::Mesh::STKMesh *mesh, int gcw, std::vector< stk::mesh::Entity*> &entities );

    //! Clone the iterator
    virtual MeshIterator* clone() const;

friend class AMP::Mesh::STKMesh;

private:
    // Data members
    int    d_gcw;
    int    d_dim;
    int    d_rank;
    MeshID d_meshID;
    const AMP::Mesh::STKMesh *d_mesh;
    std::vector<stk::mesh::Entity*> d_entries;
    std::vector<stk::mesh::Entity*>::iterator d_pos;
    MeshElement  d_cur_element;
};


}
}

#endif
