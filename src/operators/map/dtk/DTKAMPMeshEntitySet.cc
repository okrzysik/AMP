
#include "ampmesh/MeshID.h"

#include "DTKAMPMeshEntity.h"
#include "DTKAMPMeshEntityExtraData.h"
#include "DTKAMPMeshEntityIterator.h"
#include "DTKAMPMeshEntitySet.h"

#ifdef USE_EXT_MPI
#include "Teuchos_DefaultMpiComm.hpp"
#else
#include "Teuchos_DefaultSerialComm.hpp"
#endif

namespace AMP {
namespace Operator {


//---------------------------------------------------------------------------//
// Constructor.
AMPMeshEntitySet::AMPMeshEntitySet( const AMP::shared_ptr<AMP::Mesh::Mesh> &mesh )
    : d_amp_mesh( mesh )
{ /* ... */
}

//---------------------------------------------------------------------------//
// Get the parallel communicator for the entity set.
Teuchos::RCP<const Teuchos::Comm<int>> AMPMeshEntitySet::communicator() const
{
#ifdef USE_EXT_MPI
    return Teuchos::rcp( new Teuchos::MpiComm<int>( d_amp_mesh->getComm().getCommunicator() ) );
#else
    return Teuchos::rcp( new Teuchos::SerialComm<int>() );
#endif
}

//---------------------------------------------------------------------------//
// Return the largest physical dimension of the entities in the set.
int AMPMeshEntitySet::physicalDimension() const { return d_amp_mesh->getDim(); }

//---------------------------------------------------------------------------//
// Given an EntityId, get the entity.
void AMPMeshEntitySet::getEntity( const DataTransferKit::EntityType entity_type,
                                  const DataTransferKit::EntityId entity_id,
                                  DataTransferKit::Entity &entity ) const
{
    // We will only access local owned entities through this interface.
    bool is_local               = true;
    AMP::Mesh::GeomType type_id = getGeomTypeFromEntityType( entity_type );
    unsigned int local_id       = entity_id & 0x00000000FFFFFFFF;
    unsigned int owner_rank     = entity_id >> 32;
    owner_rank                  = ( owner_rank >> 8 ) & 0x007FFFFF;
    AMP::Mesh::MeshID mesh_id   = d_amp_mesh->meshID();
    AMP::Mesh::MeshElementID element_id( is_local, type_id, local_id, owner_rank, mesh_id );
    entity = AMPMeshEntity( d_amp_mesh->getElement( element_id ) );
}

//---------------------------------------------------------------------------//
// Get a iterator of the given entity type that satisfy the given predicate.
DataTransferKit::EntityIterator AMPMeshEntitySet::entityIterator(
    const DataTransferKit::EntityType entity_type,
    const std::function<bool( DataTransferKit::Entity )> &predicate ) const
{
    AMP::Mesh::GeomType type_id = getGeomTypeFromEntityType( entity_type );
    int gcw                     = 1;
    return AMPMeshEntityIterator( d_amp_mesh->getIterator( type_id, gcw ), predicate );
}

//---------------------------------------------------------------------------//
// Given an entity, get the entities of the given type that are adjacent to
// it.
void AMPMeshEntitySet::getAdjacentEntities(
    const DataTransferKit::Entity &entity,
    const DataTransferKit::EntityType entity_type,
    Teuchos::Array<DataTransferKit::Entity> &adjacent_entities ) const
{
    AMP::Mesh::GeomType type_id = getGeomTypeFromEntityType( entity_type );
    std::vector<AMP::Mesh::MeshElement> adjacent_elements =
        Teuchos::rcp_dynamic_cast<AMPMeshEntityExtraData>( entity.extraData() )
            ->d_element.getElements( type_id );
    int num_entities = adjacent_elements.size();
    adjacent_entities.resize( num_entities );
    for ( int i = 0; i < num_entities; ++i ) {
        adjacent_entities[i] = AMPMeshEntity( adjacent_elements[i] );
    }
}

//---------------------------------------------------------------------------//
// Given a DTK entity type, get an AMP GeomType.
AMP::Mesh::GeomType
AMPMeshEntitySet::getGeomTypeFromEntityType( const DataTransferKit::EntityType entity_type ) const
{
    AMP::Mesh::GeomType type_id = AMP::Mesh::Vertex;
    switch ( entity_type ) {
    case DataTransferKit::ENTITY_TYPE_NODE: type_id   = AMP::Mesh::Vertex; break;
    case DataTransferKit::ENTITY_TYPE_EDGE: type_id   = AMP::Mesh::Edge; break;
    case DataTransferKit::ENTITY_TYPE_FACE: type_id   = AMP::Mesh::Face; break;
    case DataTransferKit::ENTITY_TYPE_VOLUME: type_id = AMP::Mesh::Volume; break;
    default: type_id                                  = AMP::Mesh::null; break;
    }
    return type_id;
}

//---------------------------------------------------------------------------//
}
}
