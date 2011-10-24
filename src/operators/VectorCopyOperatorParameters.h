#ifndef included_AMP_VectorCopyOperatorParameters
#define included_AMP_VectorCopyOperatorParameters

#include <vector>

#include "vectors/Vector.h"
#include "operators/OperatorParameters.h"
#include "vectors/Variable.h"

namespace AMP {
namespace Operator {

class VectorCopyOperatorParameters: public OperatorParameters {
public:

    VectorCopyOperatorParameters(const boost::shared_ptr<
            AMP::Database> &db) :
        OperatorParameters(db)
    {
    }

    ~VectorCopyOperatorParameters() {
    }

    AMP::LinearAlgebra::Variable::shared_ptr d_copyVariable;

    AMP::LinearAlgebra::Vector::shared_ptr d_copyVector;

protected:

private:

};

}
}

#endif

