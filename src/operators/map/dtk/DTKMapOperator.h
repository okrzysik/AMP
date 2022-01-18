
#ifndef included_AMP_DTK_DTKMapOperator
#define included_AMP_DTK_DTKMapOperator

#include "AMP/operators/Operator.h"

#include "AMP/utils/AMP_MPI.h"

#include "AMP/mesh/Mesh.h"

#include "AMP/vectors/Vector.h"

#include "DTKAMPMeshManager.h"

#include <DTK_MapOperator.hpp>

namespace AMP::Operator {

class DTKMapOperatorParameters : public OperatorParameters
{
public:
    // Constructor.
    explicit DTKMapOperatorParameters( std::shared_ptr<AMP::Database> db )
        : OperatorParameters( db )
    { /* ... */
    }

    AMP_MPI d_globalComm;
    // Domain mesh. A manager for the mesh that is the data source.
    std::shared_ptr<AMP::Mesh::Mesh> d_domain_mesh;

    // Range mesh. A manager for the mesh that is the data target.
    std::shared_ptr<AMP::Mesh::Mesh> d_range_mesh;

    // Domain DOF manager. A DOF manager for the source data.
    std::shared_ptr<AMP::Discretization::DOFManager> d_domain_dofs;

    // Range DOF manager. A DOF Manager for the target data.
    std::shared_ptr<AMP::Discretization::DOFManager> d_range_dofs;
};


/**
 * AMP operator element implementation for DTK Map operator.
 */
class DTKMapOperator : public Operator
{
public:
    /**
     * Constructor.
     */
    explicit DTKMapOperator( std::shared_ptr<const OperatorParameters> params );

    //! Destructor
    ~DTKMapOperator() {}

    //! Apply function.
    void apply( AMP::LinearAlgebra::Vector::const_shared_ptr u,
                AMP::LinearAlgebra::Vector::shared_ptr r ) override;

private:
    AMP_MPI d_comm;
    Teuchos::RCP<const Teuchos::Comm<int>> d_TeuchosComm;
    bool d_mapOnThisProc;
    // DTK map operator.
    std::shared_ptr<DataTransferKit::MapOperator> d_dtk_operator;

    // DTK domain mesh.
    std::shared_ptr<DTKAMPMeshManager> d_domain_mesh;

    // DTK range mesh.
    std::shared_ptr<DTKAMPMeshManager> d_range_mesh;
};
} // namespace AMP::Operator

#endif
