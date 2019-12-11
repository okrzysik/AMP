
#ifndef included_AMP_DTK_AMPMeshManager
#define included_AMP_DTK_AMPMeshManager

#include "AMP/ampmesh/Mesh.h"
#include "AMP/discretization/DOF_Manager.h"

#include "AMP/utils/AMP_MPI.h"

#include <DTK_ClientManager.hpp>
#include <DTK_FunctionSpace.hpp>
#include <DTK_Types.hpp>

namespace AMP {
namespace Operator {


/**
 * AMP Mesh manager for DTK.
 */
class DTKAMPMeshManager : public DataTransferKit::ClientManager
{
public:
    /**
     * Constructor.
     */
    DTKAMPMeshManager( const std::shared_ptr<AMP::Mesh::Mesh> &mesh,
                       const std::shared_ptr<AMP::Discretization::DOFManager> &dof_manager,
                       const std::function<bool( DataTransferKit::Entity )> &predicate );

    /*!
     * \brief Get the function space over which the mesh and its fields are
     * defined.
     */
    Teuchos::RCP<DataTransferKit::FunctionSpace> functionSpace() const { return d_function_space; }

    //@{
    //! ClientManager interface implementation.
    /*!
     * \brief Get the entity set over which the fields are defined.
     */
    Teuchos::RCP<DataTransferKit::EntitySet> entitySet() const override
    {
        return d_function_space->entitySet();
    }

    /*!
     * \brief Get the local map for entities supporting the function.
     */
    Teuchos::RCP<DataTransferKit::EntityLocalMap> localMap() const override
    {
        return d_function_space->localMap();
    }

    /*!
     * \brief Get the shape function for entities supporting the function.
     */
    Teuchos::RCP<DataTransferKit::EntityShapeFunction> shapeFunction() const override
    {
        return d_function_space->shapeFunction();
    }

    /*!
     * \brief Get the integration rule for entities supporting the function.
     */
    Teuchos::RCP<DataTransferKit::EntityIntegrationRule> integrationRule() const override
    {
        return d_function_space->integrationRule();
    }

    /*!
     * \brief Get the selector function.
     */
    DataTransferKit::PredicateFunction selectFunction() const override
    {
        return d_function_space->selectFunction();
    }

    /*!
     * \brief Get the field for the given string key.
     */
    Teuchos::RCP<DataTransferKit::Field> field( const std::string &field_name ) const override
    {
        return Teuchos::null;
    }
    //@}

private:
    // The function space over which the mesh and its fields are defined.
    Teuchos::RCP<DataTransferKit::FunctionSpace> d_function_space;
};
} // namespace Operator
} // namespace AMP

#endif
