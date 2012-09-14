#ifndef included_AMP_STKMeshIterators
#define included_AMP_STKMeshIterators

#include <boost/shared_ptr.hpp>
#include "ampmesh/STKmesh/STKMesh.h"
#include "ampmesh/MeshIterator.h"

namespace AMP {
namespace Mesh {


class STKMeshIterator: public MeshIterator {
public:

    typedef boost::shared_ptr<std::vector<stk::mesh::Entity*> > MeshPtr ;

    //! Empty MeshIterator constructor
    STKMeshIterator();

    //! Deconstructor
    virtual ~STKMeshIterator ();

    //! Copy constructor
    STKMeshIterator(const STKMeshIterator&);

    //! Assignment operator
    STKMeshIterator& operator=(const STKMeshIterator&);

    //! Increment
    virtual MeshIterator& operator++();

    //! Increment
    virtual MeshIterator operator++(int);

    //! Decrement
    virtual MeshIterator& operator--();

    //! Decrement
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
    virtual MeshIterator end() const ;

    //! Return the number of elements in the iterator
    virtual size_t size() const;

    //! Return the current position (from the beginning) in the iterator
    virtual size_t position() const;

protected:
    /** Default constructor
     * \param entities  Entity type:  0: node, 1: element
     * \param mesh      Pointer to the STKMesh mesh
     * \param gcw       gcw to use
     *        begin     Pointer to iterator with the begining position
     *        end       Pointer to iterator with the end position
     *        pos       Pointer to iterator with the current position
     *        size      Number of elements in the iterator (-1: unknown)
     *        pos2      Index of the current position in the iterator (-1: unknown)
     */
    STKMeshIterator(const AMP::Mesh::STKMesh *mesh, int gcw, std::vector< stk::mesh::Entity*> &entities );
    STKMeshIterator(const AMP::Mesh::STKMesh *mesh, int gcw, MeshPtr entities );

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
    MeshPtr d_entries;
    std::vector<stk::mesh::Entity*>::iterator d_pos;
    MeshElement  d_cur_element;
};


}
}

#endif
