
#include "DTKAMPMeshEntity.h"
#include "DTKAMPMeshEntityImpl.h"

namespace AMP::Operator {


/**
 * AMP Mesh element implementation for DTK Entity interface.
 */

// Constructor.
AMPMeshEntity::AMPMeshEntity(
    const AMP::Mesh::MeshElement &element,
    const std::unordered_map<int, int> &rank_map,
    const std::map<AMP::Mesh::MeshElementID, DataTransferKit::EntityId> &id_map )
{
    this->b_entity_impl = Teuchos::rcp( new AMPMeshEntityImpl( element, rank_map, id_map ) );
}
} // namespace AMP::Operator
