#include <limits>

#include "DTKAMPMeshEntityImpl.h"

namespace AMP {
namespace Operator {


//---------------------------------------------------------------------------//
// Constructor.
AMPMeshEntityImpl::AMPMeshEntityImpl( const AMP::Mesh::MeshElement &element )
{
    d_extra_data = Teuchos::rcp( new AMPMeshEntityExtraData( element ) );

    // Make a unique 64-bit id.
    // The first 1 bits are 0.
    // The next 23 bits are the owning processor id
    // The next  8 bits are the element type
    // The next 32 bits are the local id
    unsigned int tmp = 0x00000000;

    // Add the owner rank
    int owner_rank = element.globalID().owner_rank();
    tmp += ( 0x007FFFFF & owner_rank ) << 8;

    // Add the type id
    AMP::Mesh::GeomType type = element.globalID().type();
    tmp += ( (unsigned char) type );

    // Add the local id
    unsigned int local_id = element.globalID().local_id();
    d_id = ( ( (AMP::Mesh::uint64) tmp ) << 32 ) + ( (AMP::Mesh::uint64) local_id );
}

//---------------------------------------------------------------------------//
/*!
 * \brief Get the entity type.
 * \return The entity type.
 */
DataTransferKit::EntityType AMPMeshEntityImpl::entityType() const
{
    DataTransferKit::EntityType entity_type;
    Mesh::GeomType geom_type = d_extra_data->d_element.elementType();
    switch ( geom_type ) {
    case Mesh::Vertex:
        entity_type = DataTransferKit::ENTITY_TYPE_NODE;
        break;
    case Mesh::Edge:
        entity_type = DataTransferKit::ENTITY_TYPE_EDGE;
        break;
    case Mesh::Face:
        entity_type = DataTransferKit::ENTITY_TYPE_FACE;
        break;
    case Mesh::Volume:
        entity_type = DataTransferKit::ENTITY_TYPE_VOLUME;
        break;
    default:
        AMP_INSIST( geom_type == Mesh::Vertex || geom_type == Mesh::Edge ||
                        geom_type == Mesh::Face || geom_type == Mesh::Volume,
                    "Invalid geometry type!" );
        entity_type = DataTransferKit::ENTITY_TYPE_NODE;
        break;
    }
    return entity_type;
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
    return static_cast<int>( d_extra_data->d_element.globalID().owner_rank() );
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
