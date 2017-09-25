#include "LinearOperator.h"
#include "utils/Utilities.h"


namespace AMP {
namespace Operator {


LinearOperator::LinearOperator( const AMP::shared_ptr<OperatorParameters> &params )
    : Operator( params )
{
    d_matrix.reset();
}


LinearOperator::LinearOperator() : Operator() { d_matrix.reset(); }


AMP::shared_ptr<AMP::LinearAlgebra::Matrix> LinearOperator::getMatrix() { return d_matrix; }


void LinearOperator::setMatrix( AMP::shared_ptr<AMP::LinearAlgebra::Matrix> in_mat )
{
    d_matrix = in_mat;
}


void LinearOperator::apply( AMP::LinearAlgebra::Vector::const_shared_ptr u,
                            AMP::LinearAlgebra::Vector::shared_ptr f )
{
    AMP_INSIST( ( ( u.get() ) != nullptr ), "NULL Solution Vector" );
    AMP_INSIST( ( ( f.get() ) != nullptr ), "NULL Residual Vector" );
    AMP_INSIST( ( ( d_matrix.get() ) != nullptr ), "NULL Matrix" );

    AMP::LinearAlgebra::Vector::const_shared_ptr uInternal = subsetInputVector( u );
    AMP::LinearAlgebra::Vector::shared_ptr fInternal       = subsetOutputVector( f );

    AMP_INSIST( ( uInternal.get() != nullptr ), "uInternal is NULL" );
    AMP_INSIST( ( fInternal.get() != nullptr ), "fInternal is NULL" );

    d_matrix->mult( uInternal, fInternal );

    fInternal->makeConsistent( AMP::LinearAlgebra::Vector::ScatterType::CONSISTENT_SET );
}
} // namespace Operator
} // namespace AMP
