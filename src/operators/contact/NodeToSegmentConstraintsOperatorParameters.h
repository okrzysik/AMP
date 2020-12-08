
#ifndef included_AMP_NodeToSegmentConstraintsOperatorParameters
#define included_AMP_NodeToSegmentConstraintsOperatorParameters

#include "AMP/discretization/DOF_Manager.h"
#include "AMP/operators/OperatorParameters.h"

namespace AMP {
namespace Operator {

class NodeToSegmentConstraintsOperatorParameters : public OperatorParameters
{
public:
    NodeToSegmentConstraintsOperatorParameters( std::shared_ptr<AMP::Database> db )
        : OperatorParameters( db )
    {
    }

    virtual ~NodeToSegmentConstraintsOperatorParameters() {}

    AMP::AMP_MPI d_GlobalComm;

    size_t d_DOFsPerNode;
    std::shared_ptr<AMP::Discretization::DOFManager> d_DOFManager;

    AMP::Mesh::MeshID d_MasterMeshID;
    AMP::Mesh::MeshID d_SlaveMeshID;

    int d_MasterBoundaryID;
    int d_SlaveBoundaryID;

protected:
private:
};
} // namespace Operator
} // namespace AMP

#endif
