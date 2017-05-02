#ifndef included_EpetraMatrixOperator_h
#define included_EpetraMatrixOperator_h

#include "operators/trilinos/EpetraMatrixOperatorParameters.h"
#include "operators/LinearOperator.h"
#include "matrices/trilinos/ManagedEpetraMatrix.h"

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

    void setVariables( AMP::LinearAlgebra::Variable::shared_ptr in,
                       AMP::LinearAlgebra::Variable::shared_ptr out )
    {
        d_Input  = in;
        d_Output = out;
    }

    virtual AMP::LinearAlgebra::Variable::shared_ptr getInputVariable() { return d_Input; }
    virtual AMP::LinearAlgebra::Variable::shared_ptr getOutputVariable() { return d_Output; }

    virtual ~EpetraMatrixOperator() {}
};
}
}

#endif
