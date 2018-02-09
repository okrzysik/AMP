#ifndef included_AMP_STKMeshIterators
#define included_AMP_STKMeshIterators

#include "AMP/ampmesh/MeshIterator.h"
#include "AMP/ampmesh/STKmesh/STKMesh.h"
#include "AMP/utils/shared_ptr.h"

namespace AMP {
namespace Mesh {


class STKMeshIterator : public MeshIterator
{
public:
    typedef AMP::shared_ptr<std::vector<stk::mesh::Entity *>> MeshPtr;

    //! Empty MeshIterator constructor
    STKMeshIterator();

    //! Deconstructor
    virtual ~STKMeshIterator();

    //! Copy constructor
    STKMeshIterator( const STKMeshIterator & );

    //! Assignment operator
    STKMeshIterator &operator=( const STKMeshIterator & );

    //! Increment
    virtual MeshIterator &operator++() override;

    //! Increment
    virtual MeshIterator operator++(int) override;

    //! Decrement
    virtual MeshIterator &operator--() override;

    //! Decrement
    virtual MeshIterator operator--(int) override;

    //! Check if two iterators are equal
    virtual bool operator==( const MeshIterator &rhs ) const override;

    //! Check if two iterators are not equal
    virtual bool operator!=( const MeshIterator &rhs ) const override;

    //! Dereference the iterator
    virtual MeshElement &operator*(void) override;

    //! Dereference the iterator
    virtual MeshElement *operator->(void) override;

    //! Dereference the iterator
    virtual const MeshElement &operator*(void) const override;

    //! Dereference the iterator
    virtual const MeshElement *operator->(void) const override;

    //! Return an iterator to the begining
    virtual MeshIterator begin() const override;

    //! Return an iterator to the begining
    virtual MeshIterator end() const override;

    //! Return the number of elements in the iterator
    virtual size_t size() const override;

    //! Return the current position (from the beginning) in the iterator
    virtual size_t position() const override;

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
    STKMeshIterator( const AMP::Mesh::STKMesh *mesh,
                     int gcw,
                     std::vector<stk::mesh::Entity *> &entities );
    STKMeshIterator( const AMP::Mesh::STKMesh *mesh, int gcw, MeshPtr entities );

    //! Clone the iterator
    virtual MeshIterator *clone() const override;

    friend class AMP::Mesh::STKMesh;

private:
    // Data members
    int d_gcw;
    int d_dim;
    int d_rank;
    MeshID d_meshID;
    const AMP::Mesh::STKMesh *d_mesh;
    MeshPtr d_entries;
    std::vector<stk::mesh::Entity *>::iterator d_pos;
    MeshElement d_cur_element;
};
} // namespace Mesh
} // namespace AMP

#endif
