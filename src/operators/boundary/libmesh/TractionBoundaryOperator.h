
#ifndef included_AMP_TractionBoundaryOperator
#define included_AMP_TractionBoundaryOperator

#include "AMP/operators/boundary/BoundaryOperator.h"
#include "AMP/operators/boundary/libmesh/TractionBoundaryOperatorParameters.h"

namespace AMP {
namespace Operator {

class TractionBoundaryOperator : public BoundaryOperator
{
public:
    explicit TractionBoundaryOperator(
        const std::shared_ptr<TractionBoundaryOperatorParameters> &params );

    virtual ~TractionBoundaryOperator() {}

    void apply( AMP::LinearAlgebra::Vector::const_shared_ptr,
                AMP::LinearAlgebra::Vector::shared_ptr f ) override;

    void addRHScorrection( AMP::LinearAlgebra::Vector::shared_ptr rhs ) override;

protected:
    AMP::LinearAlgebra::Vector::shared_ptr
    mySubsetVector( AMP::LinearAlgebra::Vector::shared_ptr vec,
                    AMP::LinearAlgebra::Variable::shared_ptr var );

    void computeCorrection();

    AMP::LinearAlgebra::Variable::shared_ptr d_var;
    AMP::LinearAlgebra::Vector::shared_ptr d_correction;
    std::vector<double> d_traction;
    std::vector<double> d_volumeElements;
    std::vector<unsigned int> d_sideNumbers;
    std::vector<AMP::Mesh::MeshElementID> d_nodeID;
    bool d_residualMode;
};
} // namespace Operator
} // namespace AMP

#endif
