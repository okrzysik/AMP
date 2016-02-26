
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
void AMPMeshEntitySet::getEntity( const DataTransferKit::EntityId entity_id,
				  const int topological_dimension,
                                  DataTransferKit::Entity &entity ) const
{
    // We will only access local owned entities through this interface.
    bool is_local               = true;
    AMP::Mesh::GeomType type_id = getGeomTypeFromEntityType( topological_dimension );
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
    const int topological_dimension,
    const DataTransferKit::PredicateFunction& predicate ) const
{
    AMP::Mesh::GeomType type_id = getGeomTypeFromEntityType( topological_dimension );
    int gcw                     = 1;
    return AMPMeshEntityIterator( d_amp_mesh->getIterator( type_id, gcw ), predicate );
}

//---------------------------------------------------------------------------//
// Given an entity, get the entities of the given type that are adjacent to
// it.
void AMPMeshEntitySet::getAdjacentEntities(
    const DataTransferKit::Entity &entity,
    const int adjacent_dimension,
    Teuchos::Array<DataTransferKit::Entity> &adjacent_entities ) const
{
    AMP::Mesh::GeomType type_id = getGeomTypeFromEntityType( adjacent_dimension );
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
AMPMeshEntitySet::getGeomTypeFromEntityType( const int topological_dimension ) const
{
    AMP::Mesh::GeomType type_id = AMP::Mesh::Vertex;
    switch ( topological_dimension ) {
	case (0):
	    type_id = AMP::Mesh::Vertex;
	    break;
	case (1):
	    type_id = AMP::Mesh::Edge;
	    break;
	case (2):
	    type_id = AMP::Mesh::Face;
	    break;
	case (3):
	    type_id = AMP::Mesh::Volume;
	    break;
	default:
	    type_id = AMP::Mesh::null;
	    break;
    }
    return type_id;
}

//---------------------------------------------------------------------------//
}
}
