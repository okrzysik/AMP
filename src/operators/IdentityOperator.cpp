#include "AMP/operators/IdentityOperator.h"
#include "AMP/utils/Utilities.h"


namespace AMP::Operator {


IdentityOperator::IdentityOperator() : LinearOperator() {}

IdentityOperator::IdentityOperator( std::shared_ptr<const OperatorParameters> params )
    : LinearOperator( params )
{
    reset( params );
}

void IdentityOperator::reset( std::shared_ptr<const OperatorParameters> params )
{
    d_memory_location = params->d_memory_location;
    if ( params->d_db ) {
        if ( params->d_db->keyExists( "InputVariable" ) ) {
            std::string inpVar = params->d_db->getString( "InputVariable" );
            d_inputVariable.reset( new AMP::LinearAlgebra::Variable( inpVar ) );
        }
        if ( params->d_db->keyExists( "OutputVariable" ) ) {
            std::string outVar = params->d_db->getString( "OutputVariable" );
            d_outputVariable.reset( new AMP::LinearAlgebra::Variable( outVar ) );
        }
    }
}

void IdentityOperator::setMatrix( std::shared_ptr<AMP::LinearAlgebra::Matrix> )
{
    AMP_ERROR( "setMatrix is invalid for the Identity operator" );
}


void IdentityOperator::apply( AMP::LinearAlgebra::Vector::const_shared_ptr u,
                              AMP::LinearAlgebra::Vector::shared_ptr r )
{
    AMP_INSIST( ( ( u.get() ) != nullptr ), "NULL Solution Vector" );
    AMP_INSIST( ( ( r.get() ) != nullptr ), "NULL Residual Vector" );

    AMP::LinearAlgebra::Vector::const_shared_ptr uInternal = subsetInputVector( u );
    AMP::LinearAlgebra::Vector::shared_ptr rInternal       = subsetOutputVector( r );

    AMP_INSIST( ( uInternal ), "uInternal is NULL" );
    AMP_INSIST( ( rInternal ), "rInternal is NULL" );

    rInternal->copyVector( uInternal );

    rInternal->makeConsistent( AMP::LinearAlgebra::ScatterType::CONSISTENT_SET );
}
} // namespace AMP::Operator
