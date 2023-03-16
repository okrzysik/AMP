#include "AMP/operators/MoveMeshOperator.h"
#include "AMP/mesh/Mesh.h"

namespace AMP::Operator {

MoveMeshOperator::MoveMeshOperator( std::shared_ptr<const OperatorParameters> params )
    : Operator( params )
{
    d_prevDisp.reset();
    d_var.reset();
}

void MoveMeshOperator::setVariable( std::shared_ptr<AMP::LinearAlgebra::Variable> var )
{
    d_var = var;
}

std::shared_ptr<AMP::LinearAlgebra::Variable> MoveMeshOperator::getInputVariable() { return d_var; }

void MoveMeshOperator::apply( AMP::LinearAlgebra::Vector::const_shared_ptr u,
                              AMP::LinearAlgebra::Vector::shared_ptr )
{
    AMP::LinearAlgebra::Vector::const_shared_ptr dispVec = u->subsetVectorForVariable( d_var );

    if ( d_prevDisp == nullptr ) {
        d_prevDisp = dispVec->clone();
        d_prevDisp->zero();
    }

    AMP::LinearAlgebra::Vector::shared_ptr deltaDisp = dispVec->clone();
    deltaDisp->subtract( *dispVec, *d_prevDisp );

    d_Mesh->displaceMesh( deltaDisp );

    d_prevDisp->copyVector( dispVec );
}
} // namespace AMP::Operator
