#ifndef included_AMP_CircleFrustumBoxMesh
#define included_AMP_CircleFrustumBoxMesh

#include "AMP/ampmesh/structured/StructuredGeometryMesh.h"

#include <array>


namespace AMP {
namespace Mesh {


/**
 * \class CircleFrustumMesh
 * \brief A derived version of BoxMesh for a square frustrum
 * \details A concrete implementation of BoxMesh for a square frustrum
 */
class CircleFrustumMesh : public AMP::Mesh::StructuredGeometryMesh
{
public:
    //! Default constructor
    explicit CircleFrustumMesh( MeshParameters::shared_ptr params );

    /**
     * \brief   Estimate the number of elements in the mesh
     * \details  This function will estimate the number of elements in the mesh.
     *   This is used so that we can properly balance the meshes across multiple processors.
     *   Ideally this should be both an accurate estimate and very fast.  It should not require
     *   any communication and should not have to actually load a mesh.
     * \param params Parameters for constructing a mesh from an input database
     */
    static std::vector<size_t> estimateLogicalMeshSize( const MeshParameters::shared_ptr &params );

    virtual AMP::shared_ptr<Mesh> clone() const override;

private:
    CircleFrustumMesh(); // Private empty constructor
};


} // namespace Mesh
} // namespace AMP


#endif
