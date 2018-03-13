#ifndef included_AMP_StructuredGeometryMesh
#define included_AMP_StructuredGeometryMesh

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
class StructuredGeometryMesh : public AMP::Mesh::BoxMesh
{
public: // Functions derived from BoxMesh
    virtual Mesh::Movable isMeshMovable() const override final;
    virtual void displaceMesh( const std::vector<double> &x ) override;
#ifdef USE_AMP_VECTORS
    virtual void displaceMesh( AMP::shared_ptr<const AMP::LinearAlgebra::Vector> ) override final;
#endif
    virtual AMP::Geometry::Point<double>
    physicalToLogical( const AMP::Geometry::Point<double> &x ) const override;
    virtual void coord( const MeshElementIndex &index, double *pos ) const override;
    virtual AMP::shared_ptr<Mesh> clone() const override = 0;

protected:
    explicit StructuredGeometryMesh( MeshParameters::shared_ptr );

private:
    explicit StructuredGeometryMesh( const BoxMesh & );
};


} // namespace Mesh
} // namespace AMP


#endif
