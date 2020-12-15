#ifndef included_AMP_StructuredGeometryMesh
#define included_AMP_StructuredGeometryMesh

#include "AMP/ampmesh/LogicalGeometry.h"
#include "AMP/ampmesh/structured/BoxMesh.h"

#include <array>
#include <vector>


namespace AMP {
namespace Mesh {


/**
 * \class StructuredGeometryMesh
 * \brief A derived version of BoxMesh for a given geometry
 * \details A concrete implementation of BoxMesh for a object which has a
 *     well defined geometry than can be used to define the mapping between
 *     physical and logical coordinates
 */
class StructuredGeometryMesh final : public AMP::Mesh::BoxMesh
{
public:
    //! Default constructor
    explicit StructuredGeometryMesh( std::shared_ptr<MeshParameters> );

    //! Copy constructor
    explicit StructuredGeometryMesh( const StructuredGeometryMesh & );

    //! Assignment operator
    StructuredGeometryMesh &operator=( const StructuredGeometryMesh & ) = delete;

public: // Functions derived from BoxMesh
    Mesh::Movable isMeshMovable() const override;
    uint64_t positionHash() const override;
    void displaceMesh( const std::vector<double> &x ) override;
#ifdef USE_AMP_VECTORS
    void displaceMesh( std::shared_ptr<const AMP::LinearAlgebra::Vector> ) override;
#endif
    AMP::Geometry::Point physicalToLogical( const AMP::Geometry::Point &x ) const override;
    void coord( const MeshElementIndex &index, double *pos ) const override;
    std::unique_ptr<Mesh> clone() const override;

private:
    uint32_t d_pos_hash;
    std::shared_ptr<AMP::Geometry::LogicalGeometry> d_geometry2;
};


} // namespace Mesh
} // namespace AMP


#endif
