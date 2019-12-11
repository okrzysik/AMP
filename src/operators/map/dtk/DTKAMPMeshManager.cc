
#include "DTKAMPMeshManager.h"
#include "DTKAMPMeshEntityLocalMap.h"
#include "DTKAMPMeshEntitySet.h"
#include "DTKAMPMeshNodalShapeFunction.h"

namespace AMP {
namespace Operator {


//---------------------------------------------------------------------------//
// Constructor.
DTKAMPMeshManager::DTKAMPMeshManager(
    const std::shared_ptr<AMP::Mesh::Mesh> &mesh,
    const std::shared_ptr<AMP::Discretization::DOFManager> &dof_manager,
    const std::function<bool( DataTransferKit::Entity )> &predicate )
{
    Teuchos::RCP<DataTransferKit::EntitySet> entity_set;
    if ( mesh )
        entity_set = Teuchos::rcp( new AMPMeshEntitySet( mesh ) );

    Teuchos::RCP<DataTransferKit::EntityLocalMap> local_map =
        Teuchos::rcp( new AMPMeshEntityLocalMap() );

    Teuchos::RCP<DataTransferKit::EntityShapeFunction> shape_function =
        Teuchos::rcp( new AMPMeshNodalShapeFunction( dof_manager ) );

    Teuchos::RCP<DataTransferKit::EntityIntegrationRule> integration_rule;

    d_function_space = Teuchos::rcp( new DataTransferKit::FunctionSpace(
        entity_set, local_map, shape_function, integration_rule, predicate ) );

    AMP_ASSERT( Teuchos::nonnull( d_function_space ) );
}

//---------------------------------------------------------------------------//
} // namespace Operator
} // namespace AMP
