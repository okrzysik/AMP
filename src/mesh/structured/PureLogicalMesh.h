#ifndef included_AMP_PureLogicalMesh
#define included_AMP_PureLogicalMesh

#include "AMP/geometry/LogicalGeometry.h"
#include "AMP/mesh/structured/BoxMesh.h"

#include <array>
#include <vector>


namespace AMP::Mesh {


/**
 * \class PureLogicalMesh
 * \brief A derived version of BoxMesh for a given geometry
 * \details A concrete implementation of BoxMesh for a object which has a
 *     well defined geometry than can be used to define the mapping between
 *     physical and logical coordinates
 */
class PureLogicalMesh final : public AMP::Mesh::BoxMesh
{
public:
    //! Default constructor
    explicit PureLogicalMesh( std::shared_ptr<const MeshParameters> );

    //! Copy constructor
    explicit PureLogicalMesh( const PureLogicalMesh & );

    //! Assignment operator
    PureLogicalMesh &operator=( const PureLogicalMesh & ) = delete;

    //! Return a string with the mesh class name
    std::string meshClass() const override;

    //! Check if two meshes are equal
    bool operator==( const Mesh &mesh ) const override;

    /**
     * \brief    Write restart data to file
     * \details  This function will write the mesh to an HDF5 file
     * \param fid    File identifier to write
     */
    void writeRestart( int64_t fid ) const override;

public: // Functions derived from BoxMesh
    Mesh::Movable isMeshMovable() const override;
    uint64_t positionHash() const override;
    void displaceMesh( const std::vector<double> &x ) override;
    void displaceMesh( std::shared_ptr<const AMP::LinearAlgebra::Vector> ) override;
    AMP::Geometry::Point physicalToLogical( const AMP::Geometry::Point &x ) const override;
    void coord( const MeshElementIndex &index, double *pos ) const override;
    std::unique_ptr<Mesh> clone() const override;


protected: // Functions derived from BoxMesh
    virtual void createBoundingBox() override;
};


} // namespace AMP::Mesh


#endif
