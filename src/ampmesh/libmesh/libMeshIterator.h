#ifndef included_AMP_libMeshIterators
#define included_AMP_libMeshIterators

#include "utils/shared_ptr.h"
#include "ampmesh/libmesh/libMesh.h"
#include "ampmesh/MeshIterator.h"

#include "libmesh/mesh.h"

namespace AMP { 
namespace Mesh {


class libMeshIterator: public MeshIterator {
public:

    //! Empty MeshIterator constructor
    libMeshIterator();

    //! Deconstructor
    virtual ~libMeshIterator ();

    //! Copy constructor
    libMeshIterator(const libMeshIterator&);

    //! Assignment operator
    libMeshIterator& operator=(const libMeshIterator&);

    // Increment
    MeshIterator& operator++();
    
    // Increment
    MeshIterator operator++(int);

    // Decrement
    MeshIterator& operator--();
    
    // Decrement
    MeshIterator operator--(int);

    // Arithmetic operator+
    virtual MeshIterator operator+(int) const;

    // Arithmetic operator+=
    virtual MeshIterator& operator+=(int N);

    // Check if two iterators are equal
    bool operator==(const MeshIterator& rhs) const;

    // Check if two iterators are not equal
    bool operator!=(const MeshIterator& rhs) const;
    
    // Dereference the iterator
    MeshElement &operator*(void);

    // Dereference the iterator
    MeshElement *operator->(void);

    // Return an iterator to the begining
    MeshIterator begin() const;

    // Return an iterator to the begining
    MeshIterator end() const;

    // Return the number of elements in the iterator
    virtual size_t size() const;

    // Return the current position (from the beginning) in the iterator
    virtual size_t position() const;

    using MeshIterator::operator+;
    using MeshIterator::operator+=;

protected:
    /** Default constructor
     * \param type      Entity type:  0: node, 1: element
     * \param mesh      Pointer to the libMesh mesh
     * \param gcw       gcw to use
     * \param begin     Pointer to iterator with the begining position
     * \param end       Pointer to iterator with the end position
     * \param pos       Pointer to iterator with the current position
     * \param size      Number of elements in the iterator (-1: unknown)
     * \param pos2      Index of the current position in the iterator (-1: unknown)
     */
    libMeshIterator(int type, const AMP::Mesh::libMesh *mesh, int gcw, void *begin, void *end, void *pos, int size=-1, int pos2=-1 );

    //! Clone the iterator
    virtual MeshIterator* clone() const;

friend class AMP::Mesh::libMesh;

private:
    // Data members
    int d_gcw;
    int d_dim;
    int d_type;
    int d_size;
    int d_pos2;
    int d_rank;
    void *d_begin;
    void *d_end;
    void *d_pos;
    MeshID d_meshID;
    const AMP::Mesh::libMesh *d_mesh;
    MeshElement d_cur_element;
};


}
}

#endif

