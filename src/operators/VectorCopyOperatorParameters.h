#ifndef included_AMP_VectorCopyOperatorParameters
#define included_AMP_VectorCopyOperatorParameters

#include <vector>

#include "AMP/operators/OperatorParameters.h"
#include "AMP/vectors/Variable.h"
#include "AMP/vectors/Vector.h"

namespace AMP {
namespace Operator {

class VectorCopyOperatorParameters : public OperatorParameters
{
public:
    explicit VectorCopyOperatorParameters( std::shared_ptr<AMP::Database> db )
        : OperatorParameters( db )
    {
    }

    virtual ~VectorCopyOperatorParameters() {}

    std::shared_ptr<AMP::LinearAlgebra::Variable> d_copyVariable;

    AMP::LinearAlgebra::Vector::shared_ptr d_copyVector;

protected:
private:
};
} // namespace Operator
} // namespace AMP

#endif
