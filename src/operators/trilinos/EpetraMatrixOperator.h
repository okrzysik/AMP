#ifndef included_AMP_EpetraMatrixOperator_h
#define included_AMP_EpetraMatrixOperator_h

#include "AMP/matrices/trilinos/ManagedEpetraMatrix.h"
#include "AMP/operators/LinearOperator.h"
#include "AMP/operators/trilinos/EpetraMatrixOperatorParameters.h"

namespace AMP {
namespace Operator {

class EpetraMatrixOperator : public LinearOperator
{
private:
    AMP::LinearAlgebra::Variable::shared_ptr d_Input, d_Output;

public:
    explicit EpetraMatrixOperator( const std::shared_ptr<EpetraMatrixOperatorParameters> &params )
        : LinearOperator( params )
    {
        AMP::LinearAlgebra::Matrix::shared_ptr t(
            new AMP::LinearAlgebra::ManagedEpetraMatrix( params->d_Matrix ) );
        setMatrix( t );
    }

    virtual ~EpetraMatrixOperator() {}
};
} // namespace Operator
} // namespace AMP

#endif
