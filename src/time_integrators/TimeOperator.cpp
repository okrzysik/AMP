#include "TimeOperator.h"
#include "AMP/utils/Database.h"
#include "TimeOperatorParameters.h"


namespace AMP::TimeIntegrator {

TimeOperator::TimeOperator( std::shared_ptr<const AMP::Operator::OperatorParameters> in_params )
    : Operator( in_params )

{
    auto params = std::dynamic_pointer_cast<const TimeOperatorParameters>( in_params );

    d_bLinearMassOperator = false;
    d_bLinearRhsOperator  = false;
    d_bAlgebraicComponent = false;

    d_dCurrentDt = 0.0;

    d_pRhsOperator  = params->d_pRhsOperator;
    d_pMassOperator = params->d_pMassOperator;

    d_pSourceTerm        = params->d_pSourceTerm;
    d_pAlgebraicVariable = params->d_pAlgebraicVariable;

    reset( in_params );
}

TimeOperator::~TimeOperator() = default;

void TimeOperator::getFromInput( std::shared_ptr<AMP::Database> db )
{
    d_bLinearMassOperator = db->getWithDefault<bool>( "bLinearMassOperator", false );
    d_bLinearRhsOperator  = db->getWithDefault<bool>( "bLinearRhsOperator", false );

    AMP_INSIST( db->keyExists( "CurrentDt" ), "key CurrentDt missing in input" );

    d_dCurrentDt = db->getScalar<double>( "CurrentDt" );

    d_bAlgebraicComponent = db->getWithDefault<bool>( "bAlgebraicComponent", false );
}

void TimeOperator::reset( std::shared_ptr<const AMP::Operator::OperatorParameters> in_params )
{
    auto params = std::dynamic_pointer_cast<const TimeOperatorParameters>( in_params );

    AMP_INSIST( params, "Error: NULL TimeOperatorParameters object" );

    getFromInput( params->d_db );

    if ( params->d_pRhsOperatorParameters ) {
        d_pRhsOperator->reset( params->d_pRhsOperatorParameters );
    }

    if ( params->d_pMassOperatorParameters ) {
        d_pMassOperator->reset( params->d_pMassOperatorParameters );
    }
}

} // namespace AMP::TimeIntegrator
