
#ifndef included_AMP_DTK_AMPMeshEntityImpl
#define included_AMP_DTK_AMPMeshEntityImpl

#include "utils/AMP_MPI.h"

#include "ampmesh/MeshElement.h"

#include <DTK_EntityImpl.hpp>
#include "DTKAMPMeshEntityExtraData.h"

namespace AMP {
namespace Operator {


/**
  * AMP Mesh element implementation for DTK EntityImpl interface.
*/
class AMPMeshEntityImpl : public DataTransferKit::EntityImpl
{
public :

    /**
     * Constructor.
     */
    explicit AMPMeshEntityImpl( const AMP::Mesh::MeshElement& element );

    //! Destructor
    ~AMPMeshEntityImpl() { }

    /*!
     * \brief Get the entity type.
     * \return The entity type.
     */
    DataTransferKit::EntityType entityType() const;

    /*!
     * \brief Get the unique global identifier for the entity.
     * \return A unique global identifier for the entity.
     */
    DataTransferKit::EntityId id() const;
    
    /*!
     * \brief Get the parallel rank that owns the entity.
     * \return The parallel rank that owns the entity.
     */
    int ownerRank() const;

    /*!
     * \brief Return the physical dimension of the entity.
     * \return The physical dimension of the entity. Any physical coordinates
     * describing the entity will be of this dimension.
     */
    int physicalDimension() const;

    /*!
     * \brief Return the Cartesian bounding box around an entity.
     * \param bounds The bounds of the box
     * (x_min,y_min,z_min,x_max,y_max,z_max).
     */
    void boundingBox( Teuchos::Tuple<double,6>& bounds ) const;

    /*!
     * \brief Determine if an entity is in the block with the given id.
     */
    bool inBlock( const int block_id ) const;

    /*!
     * \brief Determine if an entity is on the boundary with the given id.
     */
    bool onBoundary( const int boundary_id ) const;

    /*!
     * \brief Get the extra data on the entity.
     */
    Teuchos::RCP<DataTransferKit::EntityExtraData> extraData() const;

  private:
    
    // Mesh element extra data.
    Teuchos::RCP<AMPMeshEntityExtraData> d_extra_data;

    // Mesh element id.
    DataTransferKit::EntityId d_id;
};


} // end namespace Operator
} // end namespace AMP

#endif // end included_AMP_DTK_AMPMeshEntityImpl


