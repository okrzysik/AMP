
#ifndef included_AMP_CustomConstraintsEliminationOperator
#define included_AMP_CustomConstraintsEliminationOperator

#include "AMP/operators/ConstraintsEliminationOperator.h"

namespace AMP {
namespace Operator {

class CustomConstraintsEliminationOperator : public ConstraintsEliminationOperator
{

public:
    explicit CustomConstraintsEliminationOperator(
        std::shared_ptr<const OperatorParameters> params );

    void addSlaveToMaster( AMP::LinearAlgebra::Vector::shared_ptr u ) override;

    void copyMasterToSlave( AMP::LinearAlgebra::Vector::shared_ptr u ) override;

    void initialize(
        std::vector<size_t> const &slaveIndices,
        std::vector<double> const &slaveShift,
        std::vector<std::vector<size_t>> const &masterIndices = std::vector<std::vector<size_t>>(),
        std::vector<std::vector<double>> const &masterCoefficients =
            std::vector<std::vector<double>>() );

protected:
    std::vector<std::vector<size_t>> d_MasterIndices;
    std::vector<std::vector<double>> d_MasterCoefficients;
};
} // namespace Operator
} // namespace AMP

#endif
