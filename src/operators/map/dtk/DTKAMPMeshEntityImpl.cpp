#include <cstdint>
#include <limits>

#include "DTKAMPMeshEntityImpl.h"

namespace AMP::Operator {


//---------------------------------------------------------------------------//
// Constructor.
AMPMeshEntityImpl::AMPMeshEntityImpl(
    const AMP::Mesh::MeshElement &element,
    const std::unordered_map<int, int> &rank_map,
    const std::map<AMP::Mesh::MeshElementID, DataTransferKit::EntityId> &id_map )
    : d_extra_data( new AMPMeshEntityExtraData( element ) )
{
    AMP_ASSERT( rank_map.count( element.globalOwnerRank() ) );
    d_owner_rank = rank_map.find( element.globalOwnerRank() )->second;

    AMP_ASSERT( id_map.count( element.globalID() ) );
    d_id = id_map.find( element.globalID() )->second;
}

//---------------------------------------------------------------------------//
/*!
 * \brief Get the entity type.
 * \return The entity type.
 */
int AMPMeshEntityImpl::topologicalDimension() const
{
    int entity_dim;
    Mesh::GeomType geom_type = d_extra_data->d_element.elementType();
    switch ( geom_type ) {
    case Mesh::GeomType::Vertex:
        entity_dim = 0;
        break;
    case Mesh::GeomType::Edge:
        entity_dim = 1;
        break;
    case Mesh::GeomType::Face:
        entity_dim = 2;
        break;
    case Mesh::GeomType::Volume:
        entity_dim = 3;
        break;
    default:
        AMP_INSIST( geom_type == Mesh::GeomType::Vertex || geom_type == Mesh::GeomType::Edge ||
                        geom_type == Mesh::GeomType::Face || geom_type == Mesh::GeomType::Volume,
                    "Invalid geometry type!" );
        entity_dim = 0;
        break;
    }
    return entity_dim;
}

//---------------------------------------------------------------------------//
/*!
 * \brief Get the unique global identifier for the entity.
 * \return A unique global identifier for the entity.
 */
DataTransferKit::EntityId AMPMeshEntityImpl::id() const { return d_id; }

//---------------------------------------------------------------------------//
/*!
 * \brief Get the parallel rank that owns the entity.
 * \return The parallel rank that owns the entity.
 */
int AMPMeshEntityImpl::ownerRank() const { return d_owner_rank; }

//---------------------------------------------------------------------------//
/*!
 * \brief Return the physical dimension of the entity.
 * \return The physical dimension of the entity. Any physical coordinates
 * describing the entity will be of this dimension.
 */
int AMPMeshEntityImpl::physicalDimension() const
{
    // Get the vertices of the element.
    std::vector<Mesh::MeshElement> vertices =
        d_extra_data->d_element.getElements( Mesh::GeomType::Vertex );

    // Get the dimension via the coordinates.
    int space_dim = 0;
    if ( 0 < vertices.size() ) {
        space_dim = vertices[0].coord().size();
    }
    return space_dim;
}

//---------------------------------------------------------------------------//
/*!
 * \brief Return the Cartesian bounding box around an entity.
 * \param bounds The bounds of the box
 * (x_min,y_min,z_min,x_max,y_max,z_max).
 */
void AMPMeshEntityImpl::boundingBox( Teuchos::Tuple<double, 6> &bounds ) const
{
    // Initialize a bounding box.
    double max = std::numeric_limits<double>::max();
    bounds     = Teuchos::tuple( max, max, max, -max, -max, -max );

    // Get the vertices of the element.
    auto vertices = d_extra_data->d_element.getElements( Mesh::GeomType::Vertex );

    // Create a bounding box from the vertex coordinates.
    int num_vertices = vertices.size();
    for ( int i = 0; i < num_vertices; ++i ) {
        auto coords = vertices[i].coord();
        for ( unsigned d = 0; d < coords.size(); ++d ) {
            bounds[d]     = std::min( bounds[d], coords[d] );
            bounds[d + 3] = std::max( bounds[d + 3], coords[d] );
        }
    }
}

//---------------------------------------------------------------------------//
/*!
 * \brief Determine if entity is owned by the calling process.
 */
bool AMPMeshEntityImpl::isLocallyOwned() const
{
    return d_extra_data->d_element.globalID().is_local();
}

//---------------------------------------------------------------------------//
/*!
 * \brief Determine if an entity is in the block with the given id.
 */
bool AMPMeshEntityImpl::inBlock( const int block_id ) const
{
    return d_extra_data->d_element.isInBlock( block_id );
}

//---------------------------------------------------------------------------//
/*!
 * \brief Determine if an entity is on the boundary with the given id.
 */
bool AMPMeshEntityImpl::onBoundary( const int boundary_id ) const
{
    return d_extra_data->d_element.isOnBoundary( boundary_id );
}

//---------------------------------------------------------------------------//
/*!
 * \brief Get the extra data on the entity.
 */
Teuchos::RCP<DataTransferKit::EntityExtraData> AMPMeshEntityImpl::extraData() const
{
    return d_extra_data;
}

//---------------------------------------------------------------------------//
} // namespace AMP::Operator
