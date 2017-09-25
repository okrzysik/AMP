#ifndef included_EpetraMatrixOperator_h
#define included_EpetraMatrixOperator_h

#include "matrices/trilinos/ManagedEpetraMatrix.h"
#include "operators/LinearOperator.h"
#include "operators/trilinos/EpetraMatrixOperatorParameters.h"

namespace AMP {
namespace Operator {

class EpetraMatrixOperator : public LinearOperator
{
private:
    AMP::LinearAlgebra::Variable::shared_ptr d_Input, d_Output;

public:
    explicit EpetraMatrixOperator( const AMP::shared_ptr<EpetraMatrixOperatorParameters> &params )
        : LinearOperator( params )
    {
        AMP::LinearAlgebra::Matrix::shared_ptr t(
            new AMP::LinearAlgebra::ManagedEpetraMatrix( params->d_Matrix ) );
        setMatrix( t );
    }

    virtual ~EpetraMatrixOperator() {}
};
}
}

#endif
