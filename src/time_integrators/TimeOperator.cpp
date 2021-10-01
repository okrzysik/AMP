#include "TimeOperator.h"
#include "AMP/utils/Database.h"
#include "TimeOperatorParameters.h"


namespace AMP {
namespace TimeIntegrator {

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

    d_pSourceTerm           = params->d_pSourceTerm;
    d_pPreviousTimeSolution = params->d_pPreviousTimeSolution;
    d_pAlgebraicVariable    = params->d_pAlgebraicVariable;

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

void TimeOperator::apply( AMP::LinearAlgebra::Vector::const_shared_ptr u,
                          AMP::LinearAlgebra::Vector::shared_ptr r )
{

    // this routine evaluates a*[ ( M(u))/dt+fRhs(u) ] +b*f
    // where the time operator is given by u_t = fRhs(u)

    std::shared_ptr<AMP::LinearAlgebra::Vector> fTmp;

    AMP_INSIST( d_pMassOperator,
                "ERROR: "
                "AMP::TimeIntegrator::TimeIntegrator::TimeOperator::"
                "apply, the mass operator is NULL!" );
    AMP_INSIST( d_pRhsOperator,
                "ERROR: "
                "AMP::TimeIntegrator::TimeIntegrator::TimeOperator::"
                "apply, the rhs operator is NULL!" );

    if ( u )
        AMP_ASSERT( u->getUpdateStatus() ==
                    AMP::LinearAlgebra::VectorData::UpdateState::UNCHANGED );

    AMP_INSIST( ( r != nullptr ), "NULL Residual/Output Vector" );
    AMP::LinearAlgebra::Vector::shared_ptr rInternal = this->subsetOutputVector( r );
    AMP_INSIST( ( rInternal != nullptr ), "NULL Residual/Output Vector" );

    d_pScratchVector = rInternal->cloneVector();
    d_pScratchVector->zero();

    d_pMassOperator->apply( u, rInternal );
    rInternal->scale( 1.0 / d_dCurrentDt );

    d_pRhsOperator->apply( u, d_pScratchVector );

    rInternal->add( *rInternal, *d_pScratchVector );

    rInternal->makeConsistent( AMP::LinearAlgebra::VectorData::ScatterType::CONSISTENT_SET );
}

std::shared_ptr<AMP::Operator::OperatorParameters>
TimeOperator::getParameters( const std::string &type,
                             AMP::LinearAlgebra::Vector::const_shared_ptr u,
                             std::shared_ptr<AMP::Operator::OperatorParameters> params )
{
    auto timeOperator_db = std::make_shared<AMP::Database>( "TimeOperatorDatabase" );
    timeOperator_db->putScalar( "CurrentDt", d_dCurrentDt );
    timeOperator_db->putScalar( "name", "TimeOperator" );
    timeOperator_db->putScalar( "bLinearMassOperator", d_bLinearMassOperator );
    timeOperator_db->putScalar( "bLinearRhsOperator", d_bLinearRhsOperator );
    timeOperator_db->putScalar( "bAlgebraicComponent", d_bAlgebraicComponent );
    timeOperator_db->putScalar( "ScalingFactor", 1.0 / d_dCurrentDt );

    auto timeOperatorParameters =
        std::make_shared<AMP::TimeIntegrator::TimeOperatorParameters>( timeOperator_db );

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
