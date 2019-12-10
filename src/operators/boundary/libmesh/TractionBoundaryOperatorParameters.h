
#ifndef included_AMP_TractionBoundaryOperatorParameters
#define included_AMP_TractionBoundaryOperatorParameters

#include "AMP/operators/OperatorParameters.h"

namespace AMP {
namespace Operator {

class TractionBoundaryOperatorParameters : public OperatorParameters
{
public:
    explicit TractionBoundaryOperatorParameters( AMP::shared_ptr<AMP::Database> db )
        : OperatorParameters( db )
    {
    }

    virtual ~TractionBoundaryOperatorParameters() {}

    std::vector<double> d_traction;
    std::vector<double> d_volumeElements;
    std::vector<unsigned int> d_sideNumbers;
    std::vector<AMP::Mesh::MeshElementID> d_nodeID;
};
} // namespace Operator
} // namespace AMP

#endif
