
#include "AMP/operators/MoveMeshOperator.h"

namespace AMP {
namespace Operator {

MoveMeshOperator::MoveMeshOperator( std::shared_ptr<const OperatorParameters> params )
    : Operator( params )
{
    d_prevDisp.reset();
    d_var.reset();
}

void MoveMeshOperator::setVariable( AMP::LinearAlgebra::Variable::shared_ptr var ) { d_var = var; }

AMP::LinearAlgebra::Variable::shared_ptr MoveMeshOperator::getInputVariable() { return d_var; }

void MoveMeshOperator::apply( AMP::LinearAlgebra::Vector::const_shared_ptr u,
                              AMP::LinearAlgebra::Vector::shared_ptr )
{
    AMP::LinearAlgebra::Vector::const_shared_ptr dispVec = u->constSubsetVectorForVariable( d_var );

    if ( d_prevDisp == nullptr ) {
        d_prevDisp = dispVec->cloneVector();
        d_prevDisp->zero();
    }

    AMP::LinearAlgebra::Vector::shared_ptr deltaDisp = dispVec->cloneVector();
    deltaDisp->subtract( *dispVec, *d_prevDisp );

    d_Mesh->displaceMesh( deltaDisp );

    d_prevDisp->copyVector( dispVec );
}
} // namespace Operator
} // namespace AMP
