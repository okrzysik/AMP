
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
    : d_amp_mesh( mesh ), d_id_maps( 4 )
{
    // Build the rank map.
    d_rank_map        = std::make_shared<std::unordered_map<int, int>>();
    auto global_ranks = d_amp_mesh->getComm().globalRanks();
    int size          = d_amp_mesh->getComm().getSize();
    for ( int n = 0; n < size; ++n ) {
        d_rank_map->emplace( global_ranks[n], n );
    }

    // Map the global ids to DTK ids.
    for ( auto &m : d_id_maps )
        m = std::make_shared<std::map<AMP::Mesh::MeshElementID, DataTransferKit::EntityId>>();
    mapGlobalIds( d_amp_mesh->getIterator( AMP::Mesh::GeomType::Vertex, 0 ), d_id_maps[0] );
    if ( (int) d_amp_mesh->getGeomType() > 0 )
        mapGlobalIds( d_amp_mesh->getIterator( AMP::Mesh::GeomType::Edge, 0 ), d_id_maps[1] );
    if ( (int) d_amp_mesh->getGeomType() > 1 )
        mapGlobalIds( d_amp_mesh->getIterator( AMP::Mesh::GeomType::Face, 0 ), d_id_maps[2] );
    if ( (int) d_amp_mesh->getGeomType() > 2 )
        mapGlobalIds( d_amp_mesh->getIterator( AMP::Mesh::GeomType::Volume, 0 ), d_id_maps[3] );
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
    entity = AMPMeshEntity(
        d_amp_mesh->getElement( element_id ), *d_rank_map, *d_id_maps[topological_dimension] );
}

//---------------------------------------------------------------------------//
// Get a iterator of the given entity type that satisfy the given predicate.
DataTransferKit::EntityIterator
AMPMeshEntitySet::entityIterator( const int topological_dimension,
                                  const DataTransferKit::PredicateFunction &predicate ) const
{
    AMP::Mesh::GeomType type_id = getGeomTypeFromEntityType( topological_dimension );
    int gcw                     = 0;
    return AMPMeshEntityIterator( d_rank_map,
                                  d_id_maps[topological_dimension],
                                  d_amp_mesh->getIterator( type_id, gcw ),
                                  predicate );
}

//---------------------------------------------------------------------------//
// Given an entity, get the entities of the given type that are adjacent to
// it.
void AMPMeshEntitySet::getAdjacentEntities(
    const DataTransferKit::Entity &entity,
    const int adjacent_dimension,
    Teuchos::Array<DataTransferKit::Entity> &adjacent_entities ) const
{
    throw std::runtime_error( "Implementation will not give correct behavior for adjacent entities "
                              "that are not locally owned" );
    /*
    AMP::Mesh::GeomType type_id = getGeomTypeFromEntityType( adjacent_dimension );
    std::vector<AMP::Mesh::MeshElement> adjacent_elements =
        Teuchos::rcp_dynamic_cast<AMPMeshEntityExtraData>( entity.extraData() )
            ->d_element.getElements( type_id );
    int num_entities = adjacent_elements.size();
    adjacent_entities.resize( num_entities );
    for ( int i = 0; i < num_entities; ++i ) {
        adjacent_entities[i] = AMPMeshEntity( adjacent_elements[i],
                          *d_rank_map,
                          *d_id_maps[adjacent_dimension] );
    }
    */
}
//---------------------------------------------------------------------------//
// Map the global ids of an iterator to DTK ids.
void AMPMeshEntitySet::mapGlobalIds(
    AMP::Mesh::MeshIterator it,
    AMP::shared_ptr<std::map<AMP::Mesh::MeshElementID, DataTransferKit::EntityId>> &id_map )
{
    int counter = 0;
    for ( it = it.begin(); it != it.end(); ++it ) {
        if ( it->globalID().is_local() ) {
            id_map->emplace( it->globalID(), counter );
            ++counter;
        }
    }
    int comm_rank = d_amp_mesh->getComm().getRank();
    int comm_size = d_amp_mesh->getComm().getSize();
    std::vector<std::size_t> offsets( comm_size, 0 );
    d_amp_mesh->getComm().allGather( id_map->size(), offsets.data() );
    for ( int n = 1; n < comm_size; ++n ) {
        offsets[n] += offsets[n - 1];
    }
    if ( comm_rank > 0 ) {
        for ( auto &i : *id_map )
            i.second += offsets[comm_rank - 1];
    }
}

//---------------------------------------------------------------------------//
// Given a DTK entity type, get an AMP GeomType.
AMP::Mesh::GeomType
AMPMeshEntitySet::getGeomTypeFromEntityType( const int topological_dimension ) const
{
    AMP::Mesh::GeomType type_id = AMP::Mesh::GeomType::Vertex;
    switch ( topological_dimension ) {
    case ( 0 ):
        type_id = AMP::Mesh::GeomType::Vertex;
        break;
    case ( 1 ):
        type_id = AMP::Mesh::GeomType::Edge;
        break;
    case ( 2 ):
        type_id = AMP::Mesh::GeomType::Face;
        break;
    case ( 3 ):
        type_id = AMP::Mesh::GeomType::Volume;
        break;
    default:
        type_id = AMP::Mesh::GeomType::null;
        break;
    }
    return type_id;
}

//---------------------------------------------------------------------------//
} // namespace Operator
} // namespace AMP
