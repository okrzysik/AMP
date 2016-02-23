#include "DTKMapOperator.h"
#include "DTKAMPMeshManager.h"
#include "DTKAMPVectorHelpers.h"

#include <DTK_ConsistentInterpolationOperator.hpp>

#include <functional>

namespace AMP {
namespace Operator {

//---------------------------------------------------------------------------//
// Constructor
DTKMapOperator::DTKMapOperator( const AMP::shared_ptr<OperatorParameters> &params )
{
    // Get the operator parameters.
    AMP::shared_ptr<DTKMapOperatorParameters> dtk_op_params =
        AMP::dynamic_pointer_cast<DTKMapOperatorParameters>( params );
    AMP_ASSERT( dtk_op_params );

    // For now, select all elements in a mesh. This needs to be refactored to
    // select elements based on criteria in the parameters.
    std::function<bool( DataTransferKit::Entity )> select_all = [=]( DataTransferKit::Entity e ) {
        return true;
    };

    // For now, select elements in the domain mesh. This needs to be refactored to
    // select elements based on criteria in the parameters.
    DataTransferKit::EntityType domain_type = DataTransferKit::ENTITY_TYPE_VOLUME;

    // Create DTK domain mesh objects.
    DTKAMPMeshManager dtk_domain_mesh(
        dtk_op_params->d_domain_mesh, dtk_op_params->d_domain_dofs, domain_type, select_all );

    // For now, select nodes in the range mesh. This needs to be refactored to
    // select elements based on criteria in the parameters.
    DataTransferKit::EntityType range_type = DataTransferKit::ENTITY_TYPE_NODE;

    // Create DTK range mesh objects.
    DTKAMPMeshManager dtk_range_mesh(
        dtk_op_params->d_range_mesh, dtk_op_params->d_range_dofs, range_type, select_all );

    // Create a map.
    Teuchos::RCP<const Tpetra::Map<int, std::size_t>> domain_map =
        DTKAMPVectorHelpers::createTpetraMapFromAMPDOFManager( dtk_op_params->d_domain_dofs );

    // Create a map.
    Teuchos::RCP<const Tpetra::Map<int, std::size_t>> range_map =
        DTKAMPVectorHelpers::createTpetraMapFromAMPDOFManager( dtk_op_params->d_range_dofs );

    // Create a DTK map operator.
    d_dtk_operator = AMP::shared_ptr<DataTransferKit::MapOperator<double>>(
        new DataTransferKit::ConsistentInterpolationOperator<double>( domain_map, range_map ) );

    // Setup the map operator. For now we use the default parameters for the
    // map operator. This needs to be refactored to gather the needed
    // parameters from the AMP parameters.
    Teuchos::RCP<Teuchos::ParameterList> dtk_parameters = Teuchos::parameterList();
    d_dtk_operator->setup(
        dtk_domain_mesh.functionSpace(), dtk_range_mesh.functionSpace(), dtk_parameters );
}

//---------------------------------------------------------------------------//
//! Apply function.
void DTKMapOperator::apply( AMP::LinearAlgebra::Vector::const_shared_ptr u,
                            AMP::LinearAlgebra::Vector::shared_ptr r)
{
    // Copy the data from the vectors into tpetra vectors.
    Teuchos::RCP<Tpetra::Vector<double, int, std::size_t>> u_tpetra =
        DTKAMPVectorHelpers::pullTpetraVectorFromAMPVector( u );
    Teuchos::RCP<Tpetra::Vector<double, int, std::size_t>> r_tpetra =
        DTKAMPVectorHelpers::pullTpetraVectorFromAMPVector( r );

    // Do the apply.
    d_dtk_operator->apply( *u_tpetra, *r_tpetra );

    // Copy the data into the result.
    DTKAMPVectorHelpers::pushTpetraVectorToAMPVector( *r_tpetra, r );
}

//---------------------------------------------------------------------------//
}
}
