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
    explicit StructuredGeometryMesh( MeshParameters::shared_ptr );

    //! Copy constructor
    explicit StructuredGeometryMesh( const StructuredGeometryMesh & );

    //! Assignment operator
    StructuredGeometryMesh &operator=( const StructuredGeometryMesh & ) = delete;

public: // Functions derived from BoxMesh
    virtual Mesh::Movable isMeshMovable() const override;
    virtual void displaceMesh( const std::vector<double> &x ) override;
#ifdef USE_AMP_VECTORS
    virtual void displaceMesh( std::shared_ptr<const AMP::LinearAlgebra::Vector> ) override;
#endif
    virtual AMP::Geometry::Point physicalToLogical( const AMP::Geometry::Point &x ) const override;
    virtual void coord( const MeshElementIndex &index, double *pos ) const override;
    virtual std::unique_ptr<Mesh> clone() const override;

private:
    std::shared_ptr<AMP::Geometry::LogicalGeometry> d_geometry2;
};


} // namespace Mesh
} // namespace AMP


#endif
