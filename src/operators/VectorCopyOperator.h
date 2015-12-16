#ifndef included_AMP_VectorCopyOperator
#define included_AMP_VectorCopyOperator

#include <algorithm>
#include <iostream>
#include <string>
#include <vector>

#include "operators/Operator.h"
#include "operators/VectorCopyOperatorParameters.h"
#include "utils/shared_ptr.h"
#include "vectors/Vector.h"

namespace AMP {
namespace Operator {

class VectorCopyOperator : public Operator
{
public:
    explicit VectorCopyOperator( const AMP::shared_ptr<VectorCopyOperatorParameters> &params );

    virtual ~VectorCopyOperator() {}

    virtual void apply( AMP::LinearAlgebra::Vector::const_shared_ptr u,
                        AMP::LinearAlgebra::Vector::shared_ptr f ) override;

    AMP::LinearAlgebra::Variable::shared_ptr getOutputVariable() override;

    AMP::LinearAlgebra::Variable::shared_ptr getInputVariable() override;

private:
    // vector to copy into
    AMP::shared_ptr<AMP::LinearAlgebra::Vector> d_copyVector;
    AMP::shared_ptr<AMP::LinearAlgebra::Variable> d_copyVariable;
};

} // namespace Operator
} // namespace AMP

#endif
