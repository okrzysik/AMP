
#include "DTKAMPMeshManager.h"
#include "DTKAMPMeshEntityLocalMap.h"
#include "DTKAMPMeshEntitySet.h"
#include "DTKAMPMeshNodalShapeFunction.h"

#include <DTK_EntitySelector.hpp>

namespace AMP {
namespace Operator {


//---------------------------------------------------------------------------//
// Constructor.
DTKAMPMeshManager::DTKAMPMeshManager(
    const AMP::shared_ptr<AMP::Mesh::Mesh> &mesh,
    const AMP::shared_ptr<AMP::Discretization::DOFManager> &dof_manager,
    const DataTransferKit::EntityType entity_type,
    const std::function<bool( DataTransferKit::Entity )> &predicate )
{
    Teuchos::RCP<DataTransferKit::EntitySelector> entity_selector =
        Teuchos::rcp( new DataTransferKit::EntitySelector( entity_type, predicate ) );

    Teuchos::RCP<DataTransferKit::EntitySet> entity_set =
        Teuchos::rcp( new AMPMeshEntitySet( mesh ) );

    Teuchos::RCP<DataTransferKit::EntityLocalMap> local_map =
        Teuchos::rcp( new AMPMeshEntityLocalMap() );

    Teuchos::RCP<DataTransferKit::EntityShapeFunction> shape_function =
        Teuchos::rcp( new AMPMeshNodalShapeFunction( dof_manager ) );

    d_function_space = Teuchos::rcp( new DataTransferKit::FunctionSpace(
        entity_set, entity_selector, local_map, shape_function ) );

    AMP_ASSERT( Teuchos::nonnull( d_function_space ) );
}

//---------------------------------------------------------------------------//
// Get the function space over which the mesh and its fields are defined.
Teuchos::RCP<DataTransferKit::FunctionSpace> DTKAMPMeshManager::functionSpace() const
{
    return d_function_space;
}

//---------------------------------------------------------------------------//
}
}
