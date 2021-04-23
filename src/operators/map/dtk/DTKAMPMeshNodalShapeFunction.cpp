
#include "DTKAMPMeshNodalShapeFunction.h"
#include "DTKAMPMeshEntity.h"
#include "DTKAMPMeshEntityExtraData.h"

#include <Intrepid_FieldContainer.hpp>
#include <Intrepid_HGRAD_HEX_C1_FEM.hpp>

namespace AMP {
namespace Operator {


//---------------------------------------------------------------------------//
// Constructor.
AMPMeshNodalShapeFunction::AMPMeshNodalShapeFunction(
    std::shared_ptr<AMP::Discretization::DOFManager> dof_manager )
    : d_dof_manager( dof_manager )
{ /* ... */
}


//---------------------------------------------------------------------------//
// Given an entity, get the ids of the degrees of freedom in the vector space
// supporting its shape function.
void AMPMeshNodalShapeFunction::entitySupportIds(
    const DataTransferKit::Entity &entity,
    Teuchos::Array<DataTransferKit::SupportId> &dof_ids ) const
{
    AMP::Mesh::MeshElement element =
        Teuchos::rcp_dynamic_cast<AMPMeshEntityExtraData>( entity.extraData() )->d_element;

    std::vector<DataTransferKit::SupportId> entity_dofs;

    std::vector<AMP::Mesh::MeshElement> vertices =
        element.getElements( AMP::Mesh::GeomType::Vertex );
    int num_nodes = vertices.size();
    dof_ids.resize( num_nodes );
    for ( int n = 0; n < num_nodes; ++n ) {
        d_dof_manager->getDOFs( vertices[n].globalID(), entity_dofs );
        AMP_INSIST( 1 == entity_dofs.size(), "Only 1 DOF id is permitted per node" );
        dof_ids[n] = entity_dofs[0];
    }
}

//---------------------------------------------------------------------------//
// Given an entity and a reference point, evaluate the shape function of the
// entity at that point.
void AMPMeshNodalShapeFunction::evaluateValue(
    const DataTransferKit::Entity &entity,
    const Teuchos::ArrayView<const double> &reference_point,
    Teuchos::Array<double> &values ) const
{
    // Get the basis for the entity. Only Hex-8 is currently supported.
    Intrepid::Basis_HGRAD_HEX_C1_FEM<double, Intrepid::FieldContainer<double>> basis;

    // Wrap the reference point.
    Teuchos::Array<int> point_dims( 2 );
    point_dims[0] = 1;
    point_dims[1] = reference_point.size();
    Intrepid::FieldContainer<double> point_container(
        point_dims, const_cast<double *>( reference_point.getRawPtr() ) );

    // Wrap the evaluations.
    values.resize( basis.getCardinality() );
    Teuchos::Array<int> value_dims( 2 );
    value_dims[0] = basis.getCardinality();
    value_dims[1] = 1;
    Intrepid::FieldContainer<double> value_container( value_dims, values.getRawPtr() );

    // Evaluate the basis function.
    basis.getValues( value_container, point_container, Intrepid::OPERATOR_VALUE );
}

//---------------------------------------------------------------------------//
// Given an entity and a reference point, evaluate the gradient of the shape
// function of the entity at that point.
void AMPMeshNodalShapeFunction::evaluateGradient(
    const DataTransferKit::Entity &entity,
    const Teuchos::ArrayView<const double> &reference_point,
    Teuchos::Array<Teuchos::Array<double>> &gradients ) const
{
    // Get the basis for the entity. Only Hex-8 is currently supported.
    Intrepid::Basis_HGRAD_HEX_C1_FEM<double, Intrepid::FieldContainer<double>> basis;

    // Wrap the reference point.
    int space_dim = reference_point.size();
    Teuchos::Array<int> point_dims( 2 );
    point_dims[0] = 1;
    point_dims[1] = space_dim;
    Intrepid::FieldContainer<double> point_container(
        point_dims, const_cast<double *>( reference_point.getRawPtr() ) );

    // Evaluate the basis function.
    int cardinality = basis.getCardinality();
    Intrepid::FieldContainer<double> grad_container( cardinality, 1, space_dim );
    basis.getValues( grad_container, point_container, Intrepid::OPERATOR_GRAD );

    // Extract the evaluations.
    gradients.resize( cardinality );
    for ( int n = 0; n < cardinality; ++n ) {
        gradients[n].resize( space_dim );
        for ( int d = 0; d < space_dim; ++d ) {
            gradients[n][d] = grad_container( n, 0, d );
        }
    }
}

//---------------------------------------------------------------------------//
} // namespace Operator
} // namespace AMP
