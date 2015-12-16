
#ifndef included_AMP_DTK_AMPMeshManager
#define included_AMP_DTK_AMPMeshManager

#include "ampmesh/Mesh.h"
#include "discretization/DOF_Manager.h"

#include "utils/AMP_MPI.h"

#include <DTK_FunctionSpace.hpp>
#include <DTK_Types.hpp>

namespace AMP {
namespace Operator {


/**
  * AMP Mesh manager for DTK.
*/
class DTKAMPMeshManager {
public:
    /**
     * Constructor.
     */
    explicit DTKAMPMeshManager( const AMP::shared_ptr<AMP::Mesh::Mesh> &mesh,
                                const AMP::shared_ptr<AMP::Discretization::DOFManager> &dof_manager,
                                const DataTransferKit::EntityType entity_type,
                                const std::function<bool( DataTransferKit::Entity )> &predicate );

    //! Destructor
    ~DTKAMPMeshManager() {}

    /*!
     * \brief Get the function space over which the mesh and its fields are
     * defined.
     */
    Teuchos::RCP<DataTransferKit::FunctionSpace> functionSpace() const;

private:
    // The function space over which the mesh and its fields are defined.
    Teuchos::RCP<DataTransferKit::FunctionSpace> d_function_space;
};
}
}

#endif
