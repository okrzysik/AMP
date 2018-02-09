#include "TimeOperator.h"
#include "AMP/utils/InputDatabase.h"
#include "TimeOperatorParameters.h"


namespace AMP {
namespace TimeIntegrator {

TimeOperator::TimeOperator( AMP::shared_ptr<AMP::Operator::OperatorParameters> in_params )
    : Operator( in_params )

{
    AMP::shared_ptr<TimeOperatorParameters> params =
        AMP::dynamic_pointer_cast<TimeOperatorParameters>( in_params );

    d_bLinearMassOperator = false;
    d_bLinearRhsOperator  = false;
    d_bAlgebraicComponent = false;

    d_dCurrentDt = 0.0;

    d_pRhsOperator  = params->d_pRhsOperator;
    d_pMassOperator = params->d_pMassOperator;

    d_pSourceTerm           = params->d_pSourceTerm;
    d_pPreviousTimeSolution = params->d_pPreviousTimeSolution;
    d_pAlgebraicVariable    = params->d_pAlgebraicVariable;

    reset( in_params );
}

TimeOperator::~TimeOperator() = default;

void TimeOperator::getFromInput( const AMP::shared_ptr<AMP::Database> &db )
{
    d_bLinearMassOperator = db->getBoolWithDefault( "bLinearMassOperator", false );
    d_bLinearRhsOperator  = db->getBoolWithDefault( "bLinearRhsOperator", false );

    AMP_INSIST( db->keyExists( "CurrentDt" ), "key CurrentDt missing in input" );

    d_dCurrentDt = db->getDouble( "CurrentDt" );

    d_bAlgebraicComponent = db->getBoolWithDefault( "bAlgebraicComponent", false );
}

void TimeOperator::reset( const AMP::shared_ptr<AMP::Operator::OperatorParameters> &in_params )
{
    AMP::shared_ptr<TimeOperatorParameters> params =
        AMP::dynamic_pointer_cast<TimeOperatorParameters>( in_params );

    AMP_INSIST( params.get() != nullptr, "Error: NULL TimeOperatorParameters object" );

    getFromInput( params->d_db );

    if ( params->d_pRhsOperatorParameters.get() != nullptr ) {
        d_pRhsOperator->reset( params->d_pRhsOperatorParameters );
    }

    if ( params->d_pMassOperatorParameters.get() != nullptr ) {
        d_pMassOperator->reset( params->d_pMassOperatorParameters );
    }
}

void TimeOperator::apply( AMP::LinearAlgebra::Vector::const_shared_ptr u,
                          AMP::LinearAlgebra::Vector::shared_ptr r )
{

    // this routine evaluates a*[ ( M(u))/dt+fRhs(u) ] +b*f
    // where the time operator is given by u_t = fRhs(u)

    AMP::shared_ptr<AMP::LinearAlgebra::Vector> fTmp;

    AMP_INSIST( d_pMassOperator.get() != nullptr,
                "ERROR: "
                "AMP::TimeIntegrator::TimeIntegrator::TimeOperator::"
                "apply, the mass operator is NULL!" );
    AMP_INSIST( d_pRhsOperator.get() != nullptr,
                "ERROR: "
                "AMP::TimeIntegrator::TimeIntegrator::TimeOperator::"
                "apply, the rhs operator is NULL!" );

    if ( u.get() != nullptr )
        AMP_ASSERT( u->getUpdateStatus() == AMP::LinearAlgebra::Vector::UpdateState::UNCHANGED );

    AMP_INSIST( ( r != nullptr ), "NULL Residual/Output Vector" );
    AMP::LinearAlgebra::Vector::shared_ptr rInternal = this->subsetOutputVector( r );
    AMP_INSIST( ( rInternal != nullptr ), "NULL Residual/Output Vector" );

    d_pScratchVector = rInternal->cloneVector();
    d_pScratchVector->zero();

    d_pMassOperator->apply( u, rInternal );
    rInternal->scale( 1.0 / d_dCurrentDt );

    d_pRhsOperator->apply( u, d_pScratchVector );

    rInternal->add( *rInternal, *d_pScratchVector );

    rInternal->makeConsistent( AMP::LinearAlgebra::Vector::ScatterType::CONSISTENT_SET );
}

AMP::shared_ptr<AMP::Operator::OperatorParameters>
TimeOperator::getParameters( const std::string &type,
                             AMP::LinearAlgebra::Vector::const_shared_ptr u,
                             AMP::shared_ptr<AMP::Operator::OperatorParameters> params )
{
    AMP::shared_ptr<AMP::InputDatabase> timeOperator_db(
        new AMP::InputDatabase( "TimeOperatorDatabase" ) );
    timeOperator_db->putDouble( "CurrentDt", d_dCurrentDt );
    timeOperator_db->putString( "name", "TimeOperator" );
    timeOperator_db->putBool( "bLinearMassOperator", d_bLinearMassOperator );
    timeOperator_db->putBool( "bLinearRhsOperator", d_bLinearRhsOperator );
    timeOperator_db->putBool( "bAlgebraicComponent", d_bAlgebraicComponent );
    timeOperator_db->putDouble( "ScalingFactor", 1.0 / d_dCurrentDt );

    AMP::shared_ptr<TimeOperatorParameters> timeOperatorParameters(
        new AMP::TimeIntegrator::TimeOperatorParameters( timeOperator_db ) );

    timeOperatorParameters->d_Mesh = d_Mesh;
    // if we have a linear rhs operator then just pass the pointer to the rhs operator instead of
    // the parameter object
    if ( d_bLinearRhsOperator ) {
        timeOperatorParameters->d_pRhsOperator = d_pRhsOperator;
    } else {
        timeOperatorParameters->d_pRhsOperatorParameters =
            d_pRhsOperator->getParameters( type, u, params );
    }

    if ( !d_bAlgebraicComponent ) {
        // if we have a linear mass operator then just pass the pointer to the mass operator instead
        // of the parameter
        // object
        if ( d_bLinearMassOperator ) {
            timeOperatorParameters->d_pMassOperator = d_pMassOperator;
        } else {
            timeOperatorParameters->d_pMassOperatorParameters =
                d_pMassOperator->getParameters( type, u, params );
        }
    }

    return timeOperatorParameters;
}
} // namespace TimeIntegrator
} // namespace AMP
