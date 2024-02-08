#include "CoupledOperator.h"
#include "AMP/utils/Utilities.h"
#include "AMP/vectors/Vector.h"

#include "ProfilerApp.h"

#include <vector>


namespace AMP::Operator {


CoupledOperator::CoupledOperator( std::shared_ptr<const OperatorParameters> params )
    : ColumnOperator( params )
{
    auto myparams = std::dynamic_pointer_cast<const CoupledOperatorParameters>( params );
    d_operators.push_back( myparams->d_NodeToGaussPointOperator );
    d_operators.push_back( myparams->d_CopyOperator );
    d_operators.push_back( myparams->d_MapOperator );
    d_operators.push_back( myparams->d_BVPOperator );
}


void CoupledOperator::apply( AMP::LinearAlgebra::Vector::const_shared_ptr u,
                             AMP::LinearAlgebra::Vector::shared_ptr r )
{
    PROFILE_START( "apply" );
    // Fill the gauss-point vector if necessary
    if ( d_operators[0] ) {
        d_operators[0]->apply( u, d_frozenGaussPointVector );
    }
    // Call copy vector
    if ( d_operators[1] ) {
        if ( d_operators[0] ) {
            d_operators[1]->apply( d_frozenGaussPointVector, r );
        } else {
            d_operators[1]->apply( u, r );
        }
    }
    // Call the map
    if ( d_operators[2] ) {
        d_operators[2]->apply( u, r );
    }
    // Call the operator
    d_operators[3]->apply( u, r );
    PROFILE_STOP( "apply" );
}

void CoupledOperator::residual( AMP::LinearAlgebra::Vector::const_shared_ptr f,
                                AMP::LinearAlgebra::Vector::const_shared_ptr u,
                                AMP::LinearAlgebra::Vector::shared_ptr r )
{
    this->apply( u, r );

    //    AMP::LinearAlgebra::Vector::shared_ptr rInternal = subsetOutputVector( r );
    //    AMP_INSIST( ( rInternal  ), "rInternal is NULL" );

    // the rhs can be NULL
    if ( f ) {
        //        AMP::LinearAlgebra::Vector::const_shared_ptr fInternal = subsetOutputVector( f );
        r->subtract( *f, *r );
    } else {
        r->scale( -1.0 );
    }

    r->makeConsistent( AMP::LinearAlgebra::ScatterType::CONSISTENT_SET );
}
} // namespace AMP::Operator
