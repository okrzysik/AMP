#ifndef included_AMP_VectorCopyOperatorParameters
#define included_AMP_VectorCopyOperatorParameters

#include <vector>

#include "operators/OperatorParameters.h"
#include "vectors/Variable.h"
#include "vectors/Vector.h"

namespace AMP {
namespace Operator {

class VectorCopyOperatorParameters : public OperatorParameters
{
public:
    VectorCopyOperatorParameters( const AMP::shared_ptr<AMP::Database> &db )
        : OperatorParameters( db )
    {
    }

    virtual ~VectorCopyOperatorParameters() {}

    AMP::LinearAlgebra::Variable::shared_ptr d_copyVariable;

    AMP::LinearAlgebra::Vector::shared_ptr d_copyVector;

protected:
private:
};
}
}

#endif
