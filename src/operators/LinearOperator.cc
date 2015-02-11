#include "LinearOperator.h"
#include "utils/Utilities.h"


namespace AMP {
namespace Operator {


LinearOperator :: LinearOperator (const AMP::shared_ptr<OperatorParameters> & params)
    : Operator (params) 
{
    d_matrix.reset();
}


LinearOperator :: LinearOperator ()
    : Operator () 
{
    d_matrix.reset();
}


AMP::shared_ptr<AMP::LinearAlgebra::Matrix> LinearOperator :: getMatrix() 
{
    return d_matrix;
}


void LinearOperator :: setMatrix( AMP::shared_ptr<AMP::LinearAlgebra::Matrix> in_mat) 
{
    d_matrix = in_mat;
}


void LinearOperator :: apply(AMP::LinearAlgebra::Vector::const_shared_ptr f, 
    AMP::LinearAlgebra::Vector::const_shared_ptr u,
    AMP::LinearAlgebra::Vector::shared_ptr r, const double a, const double b)
{
    AMP_INSIST( ((u.get()) != NULL), "NULL Solution Vector" );
    AMP_INSIST( ((r.get()) != NULL), "NULL Residual Vector" );
    AMP_INSIST( ((d_matrix.get()) != NULL), "NULL Matrix" );

    AMP::LinearAlgebra::Vector::const_shared_ptr uInternal = subsetInputVector(u);
    AMP::LinearAlgebra::Vector::shared_ptr rInternal = subsetOutputVector(r);

    AMP_INSIST( (uInternal.get() != NULL), "uInternal is NULL" );
    AMP_INSIST( (rInternal.get() != NULL), "rInternal is NULL" );

    d_matrix->mult(uInternal, rInternal);

    if(f.get() == NULL) {
        rInternal->scale(a);
    } else {
        AMP::LinearAlgebra::Vector::const_shared_ptr fInternal = subsetOutputVector(f);
        if(fInternal.get() == NULL) {
            rInternal->scale(a);
        } else {
            rInternal->axpby(b, a, fInternal);
        }
    }
    rInternal->makeConsistent(AMP::LinearAlgebra::Vector::CONSISTENT_SET);
}


}
}

