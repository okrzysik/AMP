#include <limits>
#include <cstdint>

#include "DTKAMPMeshEntityImpl.h"

namespace AMP {
namespace Operator {


//---------------------------------------------------------------------------//
// Constructor.
AMPMeshEntityImpl::AMPMeshEntityImpl( const AMP::Mesh::MeshElement &element,
				      const std::unordered_map<int,int>& rank_map )
{
    d_extra_data = Teuchos::rcp( new AMPMeshEntityExtraData( element ) );

    // Get the owner rank from the rank map.
    int global_rank = element.globalOwnerRank();
    AMP_INSIST( rank_map.count(global_rank), "Owner rank not mapped." );
    d_owner_rank = rank_map.find( global_rank )->second;
    
    // Make a unique 64-bit id.
    // It is a trade-off:
    // We have 32 bits to pack both the mesh id and the owner rank.
    // We use 18 bits for the rank and keep 14 for the mesh.
    uint32_t mesh_id = element.globalID().meshID().getData();
    AMP_INSIST(
        mesh_id < 16384,
        "Large mesh IDs not supported." );
    unsigned int owner_rank = element.globalID().owner_rank();        
    AMP_INSIST(
        owner_rank < 262144,
        "Too many processes to compose unique ids." );
    
    unsigned int local_id = element.globalID().local_id();

    d_id = ( ( (AMP::Mesh::uint64) mesh_id ) << 50 )
        + ( ( (AMP::Mesh::uint64) owner_rank ) << 32 )
        + ( (AMP::Mesh::uint64) local_id );
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
    case Mesh::Vertex:
        entity_dim = 0;
        break;
    case Mesh::Edge:
        entity_dim = 1;
        break;
    case Mesh::Face:
        entity_dim = 2;
        break;
    case Mesh::Volume:
        entity_dim = 3;
        break;
    default:
        AMP_INSIST( geom_type == Mesh::Vertex || geom_type == Mesh::Edge ||
		    geom_type == Mesh::Face || geom_type == Mesh::Volume,
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
int AMPMeshEntityImpl::ownerRank() const
{
    return d_owner_rank;
}

//---------------------------------------------------------------------------//
/*!
 * \brief Return the physical dimension of the entity.
 * \return The physical dimension of the entity. Any physical coordinates
 * describing the entity will be of this dimension.
 */
int AMPMeshEntityImpl::physicalDimension() const
{
    // Get the vertices of the element.
    std::vector<Mesh::MeshElement> vertices = d_extra_data->d_element.getElements( Mesh::Vertex );

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
    std::vector<Mesh::MeshElement> vertices = d_extra_data->d_element.getElements( Mesh::Vertex );

    // Create a bounding box from the vertex coordinates.
    int num_vertices = vertices.size();
    std::vector<double> coords;
    for ( int i = 0; i < num_vertices; ++i ) {
        coords = vertices[i].coord();
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
bool AMPMeshEntityImpl::isLocallyOwned( ) const
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
}
}
