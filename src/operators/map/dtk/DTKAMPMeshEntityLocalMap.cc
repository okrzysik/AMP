
#include "DTKAMPMeshEntityLocalMap.h"
#include "DTKAMPMeshEntityExtraData.h"

#include <Shards_BasicTopologies.hpp>
#include <Shards_CellTopology.hpp>

#include <DTK_IntrepidCellLocalMap.hpp>


namespace AMP {
namespace Operator {


//---------------------------------------------------------------------------//
// Constructor.
AMPMeshEntityLocalMap::AMPMeshEntityLocalMap() : d_inclusion_tol( 1.0e-6 ) { /* ... */}

//---------------------------------------------------------------------------//
// Set parameters for mapping.
void AMPMeshEntityLocalMap::setParameters( const Teuchos::ParameterList &parameters )
{
    if ( parameters.isParameter( "Point Inclusion Tolerance" ) ) {
        d_inclusion_tol = parameters.get<double>( "Point Inclusion Tolerance" );
    }
}

//---------------------------------------------------------------------------//
// Return the entity measure with respect to the parameteric dimension (volume
// for a 3D entity, area for 2D, and length for 1D).
double AMPMeshEntityLocalMap::measure( const DataTransferKit::Entity &entity ) const
{
    Intrepid::FieldContainer<double> entity_coords;
    getElementNodeCoordinates( entity, entity_coords );
    shards::CellTopology entity_topo( shards::getCellTopologyData<shards::Hexahedron<8>>() );
    return DataTransferKit::IntrepidCellLocalMap::measure( entity_topo, entity_coords );
}

//---------------------------------------------------------------------------//
// Return the centroid of the entity.
void AMPMeshEntityLocalMap::centroid( const DataTransferKit::Entity &entity,
                                      const Teuchos::ArrayView<double> &centroid ) const
{
    // If we have a node just return the coordinates.
    if ( 0 == entity.topologicalDimension() ) {
        std::vector<double> coords =
            Teuchos::rcp_dynamic_cast<AMPMeshEntityExtraData>( entity.extraData() )
                ->d_element.coord();
        std::copy( coords.begin(), coords.end(), centroid.begin() );
    }

    // Otherwise get the centroid of the element.
    else {
        Intrepid::FieldContainer<double> entity_coords;
        getElementNodeCoordinates( entity, entity_coords );
        shards::CellTopology entity_topo( shards::getCellTopologyData<shards::Hexahedron<8>>() );
        DataTransferKit::IntrepidCellLocalMap::centroid( entity_topo, entity_coords, centroid );
    }
}

//---------------------------------------------------------------------------//
// (Reverse Map) Map a point to the reference space of an entity. Return the
// parameterized point.
bool AMPMeshEntityLocalMap::mapToReferenceFrame(
    const DataTransferKit::Entity &entity,
    const Teuchos::ArrayView<const double> &point,
    const Teuchos::ArrayView<double> &reference_point ) const
{
    Intrepid::FieldContainer<double> entity_coords;
    getElementNodeCoordinates( entity, entity_coords );
    shards::CellTopology entity_topo( shards::getCellTopologyData<shards::Hexahedron<8>>() );
    return DataTransferKit::IntrepidCellLocalMap::mapToReferenceFrame(
        entity_topo, entity_coords, point, reference_point );
}

//---------------------------------------------------------------------------//
// Determine if a reference point is in the parameterized space of an entity.
bool AMPMeshEntityLocalMap::checkPointInclusion(
    const DataTransferKit::Entity &entity,
    const Teuchos::ArrayView<const double> &reference_point ) const
{
    Intrepid::FieldContainer<double> entity_coords;
    getElementNodeCoordinates( entity, entity_coords );
    shards::CellTopology entity_topo( shards::getCellTopologyData<shards::Hexahedron<8>>() );
    return DataTransferKit::IntrepidCellLocalMap::checkPointInclusion(
        entity_topo, reference_point, d_inclusion_tol );
}

//---------------------------------------------------------------------------//
// (Forward Map) Map a reference point to the physical space of an entity.
void AMPMeshEntityLocalMap::mapToPhysicalFrame(
    const DataTransferKit::Entity &entity,
    const Teuchos::ArrayView<const double> &reference_point,
    const Teuchos::ArrayView<double> &point ) const
{
    Intrepid::FieldContainer<double> entity_coords;
    getElementNodeCoordinates( entity, entity_coords );
    shards::CellTopology entity_topo( shards::getCellTopologyData<shards::Hexahedron<8>>() );
    DataTransferKit::IntrepidCellLocalMap::mapToPhysicalFrame(
        entity_topo, entity_coords, reference_point, point );
}

//---------------------------------------------------------------------------//
// Compute the normal on a face (3D) or edge (2D) at a given reference
// point. A default implementation is provided using a finite difference
// scheme.
void AMPMeshEntityLocalMap::normalAtReferencePoint(
    const DataTransferKit::Entity &entity,
    const DataTransferKit::Entity &parent_entity,
    const Teuchos::ArrayView<const double> &reference_point,
    const Teuchos::ArrayView<double> &normal ) const
{
    // Currently not implemented.
    bool not_implemented = true;
    AMP_ASSERT( !not_implemented );
}

//---------------------------------------------------------------------------//
// Given an entity, extract the node coordinates in canonical order.
void AMPMeshEntityLocalMap::getElementNodeCoordinates(
    const DataTransferKit::Entity &entity, Intrepid::FieldContainer<double> &entity_coords ) const
{
    // Get the vertices.
    std::vector<AMP::Mesh::MeshElement> vertices =
        Teuchos::rcp_dynamic_cast<AMPMeshEntityExtraData>( entity.extraData() )
            ->d_element.getElements( AMP::Mesh::Vertex );

    // Allocate and fill the coordinate container.
    int num_cells = 1;
    int num_nodes = vertices.size();
    int space_dim = entity.physicalDimension();
    entity_coords = Intrepid::FieldContainer<double>( num_cells, num_nodes, space_dim );
    for ( int n = 0; n < num_nodes; ++n ) {
        for ( int d = 0; d < space_dim; ++d ) {
            entity_coords( 0, n, d ) = vertices[n].coord( d );
        }
    }
}

//---------------------------------------------------------------------------//
}
}
