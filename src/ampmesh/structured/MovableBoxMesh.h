#ifndef included_AMP_MovableBoxMesh
#define included_AMP_MovableBoxMesh

#include "ampmesh/structured/BoxMesh.h"


namespace AMP {
namespace Mesh {



/**
 * \class MovableBoxMesh
 * \brief A general BoxMesh
 * \details A concrete implimentation of BoxMesh in which the
 *    coordinates of the nodes may be moved.  The base class
 *    BoxMesh is an arbitrary logically rectangular mesh class.
 *    To allow for derived implimentations which can leverage specific
 *    properties, a derived class is required.  
 */
class MovableBoxMesh : public AMP::Mesh::BoxMesh
{
public:
    //! Construct a movable box mesh from any existing box mesh
    MovableBoxMesh( const AMP::Mesh::BoxMesh& mesh );
    
    /**
     * \brief    Is the mesh movable
     * \details  This function will check if the mesh can be displaced.
     *    It will return 0 if the mesh cannont be moved, 1 if it can be displaced,
     *    and 2 if the individual nodes can be moved.
     * @return  The if
     */
    virtual int isMeshMovable( ) const override;

    /**
     * \brief    Displace the entire mesh
     * \details  This function will displace the entire mesh by a scalar value.
     *   This function is a blocking call for the mesh communicator, and requires
     *   the same value on all processors.  The displacement vector should be the
     *   size of the physical dimension.
     * \param x  Displacement vector
     */
    virtual void displaceMesh( const std::vector<double> &x ) override;

#ifdef USE_AMP_VECTORS
    /**
     * \brief    Displace the entire mesh
     * \details  This function will displace the entire mesh by displacing
     *   each node by the values provided in the vector.  This function is
     *   a blocking call for the mesh communicator
     * \param x  Displacement vector.  Must have N DOFs per node where N
     *           is the physical dimension of the mesh.
     */
    virtual void displaceMesh( AMP::shared_ptr<const AMP::LinearAlgebra::Vector> x ) override;
#endif

    //! Virtual function to copy the mesh (allows use to proply copy the derived class)
    virtual AMP::shared_ptr<Mesh> copy() const override;

    /**
     * \brief    Return a mesh element's coordinates given it's id.
     * \details  This function queries the mesh to get an element's coordinates given the mesh id.
     *    Ideally, this should be done in O(1) time, but the implimentation is up to
     *    the underlying mesh.  
     * \param[in] index     Mesh element index we are requesting.
     * \param[out] pos      Mesh element coordinates
     */
    virtual void coord( const MeshElementIndex &index, double *pos ) const override;

private:
    MovableBoxMesh(); // Private empty constructor

    // The coordinates of the nodes
    std::vector<MeshElementIndex> d_index;
    std::vector<double> d_coord[3];

    // Boundary information
    std::vector<int> d_ids;
    std::vector<MeshIterator> d_surface[4];
    std::vector<std::vector<MeshIterator>> d_boundary[4];
};


} // Mesh namespace
} // AMP namespace


#endif
