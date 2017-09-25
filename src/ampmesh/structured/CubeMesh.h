#ifndef included_AMP_CubeBoxMesh
#define included_AMP_CubeBoxMesh

#include "ampmesh/structured/BoxMesh.h"

#include <array>
#include <vector>


namespace AMP {
namespace Mesh {


/**
 * \class CubeMesh
 * \brief A derived version of BoxMesh for a cube
 * \details A concrete implementation of BoxMesh for a physically orthogonal space
 */
class CubeMesh : public AMP::Mesh::BoxMesh
{
public:
    //! Default constructor
    explicit CubeMesh( MeshParameters::shared_ptr params );

    /**
     * \brief   Estimate the number of elements in the mesh
     * \details  This function will estimate the number of elements in the mesh.
     *   This is used so that we can properly balance the meshes across multiple processors.
     *   Ideally this should be both an accurate estimate and very fast.  It should not require
     *   any communication and should not have to actually load a mesh.
     * \param params Parameters for constructing a mesh from an input database
     */
    static std::vector<size_t> estimateLogicalMeshSize( const MeshParameters::shared_ptr &params );

    /**
     * \brief    Is the mesh movable
     * \details  This function will check if the mesh can be displaced.
     *    It will return 0 if the mesh cannont be moved, 1 if it can be displaced,
     *    and 2 if the individual nodes can be moved.
     * @return  The if
     */
    virtual int isMeshMovable() const override;

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


public: // BoxMesh specific functionality
    /**
     * \brief    Return the logical coordinates
     * \details  This function queries the mesh to get the logical coordinates in [0,1]
     *     from the physical coordinates.  Not all meshes support this functionallity.
     * \param[in] x         Physical coordinates
     * @return              Returns the logical coordinates
     */
    virtual std::array<double, 3> physicalToLogical( const double *x ) const override;


private:
    CubeMesh(); // Private empty constructor

    // The coordinates of the nodes
    std::vector<double> d_coord[3];
};


} // namespace Mesh
} // namespace AMP


#endif
