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

    //    AMP_INSIST( db->keyExists( "CurrentDt" ), "key CurrentDt missing in input" );

    d_dCurrentDt = db->getWithDefault<double>( "CurrentDt", 0.0 );

    d_bAlgebraicComponent = db->getWithDefault<bool>( "bAlgebraicComponent", false );
}

void TimeOperator::reset( std::shared_ptr<const AMP::Operator::OperatorParameters> in_params )
{
    auto params = std::dynamic_pointer_cast<const TimeOperatorParameters>( in_params );

    if ( params ) {

        getFromInput( params->d_db );

        if ( params->d_pRhsOperatorParameters ) {
            d_pRhsOperator->reset( params->d_pRhsOperatorParameters );
        }

        if ( d_pMassOperator && params->d_pMassOperatorParameters ) {
            d_pMassOperator->reset( params->d_pMassOperatorParameters );
        }
    } else {
        d_pRhsOperator->reset( nullptr );
    }
}

void TimeOperator::applyRhs( std::shared_ptr<const AMP::LinearAlgebra::Vector> x,
                             std::shared_ptr<AMP::LinearAlgebra::Vector> f )
{
    AMP_INSIST( d_pRhsOperator, "RHS Operator is NULL" );
    d_pRhsOperator->apply( x, f );
}

void TimeOperator::apply( AMP::LinearAlgebra::Vector::const_shared_ptr u,
                          AMP::LinearAlgebra::Vector::shared_ptr r )
{
    // fRhs(x^{n+1})
    applyRhs( u, r );

    if ( !d_pMassOperator ) {

        // f =  x^{n+1} - \gamma*fRhs(x^{n+1}) (identity mass matrix)
        r->axpy( -d_dGamma, *r, *u );

    } else {

        if ( !d_pScratchVector ) {
            d_pScratchVector = r->cloneVector();
            d_pScratchVector->zero();
        }

        // f =  M x^{n+1}
        d_pMassOperator->apply( u, d_pScratchVector );

        if ( d_iDebugPrintInfoLevel > 4 ) {
            AMP::pout << "Output of M * yp in TimeOperator" << std::endl;
            AMP::pout << d_pScratchVector << std::endl;
        }

        if ( d_pAlgebraicVariable ) {
            auto algebraicComponent =
                d_pScratchVector->subsetVectorForVariable( d_pAlgebraicVariable );
            algebraicComponent->zero();
        }

        // f =  M x^{n+1} - \gamma*fRhs(x^{n+1})
        r->axpy( -d_dGamma, *r, *d_pScratchVector );
    }

    if ( d_iDebugPrintInfoLevel > 5 ) {
        AMP::pout << "Output of M * yp-fRhs(y,t) in TimeOperator" << std::endl;
        AMP::pout << r << std::endl;
    }

    if ( d_pSourceTerm ) {
        r->axpy( d_dGamma, *d_pSourceTerm, *r );
    }

    if ( d_iDebugPrintInfoLevel > 6 ) {
        AMP::pout << "Output of M * yp-gamma * ( fRhs(y,t)-g ) in TimeOperator" << std::endl;
        AMP::pout << r << std::endl;
    }
}

void TimeOperator::residual( std::shared_ptr<const AMP::LinearAlgebra::Vector> f,
                             std::shared_ptr<const AMP::LinearAlgebra::Vector> u,
                             std::shared_ptr<AMP::LinearAlgebra::Vector> r )
{
    apply( u, r );

    if ( f ) {
        r->axpy( -1.0, *r, *f );
    } else {
        r->scale( -1.0, *r );
    }
}

std::shared_ptr<AMP::Operator::OperatorParameters>
TimeOperator::getParameters( const std::string &type,
                             AMP::LinearAlgebra::Vector::const_shared_ptr u,
                             std::shared_ptr<AMP::Operator::OperatorParameters> params )
{
    std::shared_ptr<AMP::Operator::OperatorParameters> nestedParams;
    if ( params ) {
        nestedParams = params;
    } else {
        auto db      = std::make_shared<AMP::Database>( "nestedTimeOpParams" );
        nestedParams = std::make_shared<AMP::Operator::OperatorParameters>( db );
    }

    auto nested_db = nestedParams->d_db;
    nested_db->putScalar<bool>( "time_dependent_Jacobian", true );
    nested_db->putScalar<double>( "gamma", d_dGamma );

    AMP_INSIST( !d_pMassOperator, "Not implemented for mass operator" );
    return d_pRhsOperator->getParameters( type, u, nestedParams );
}

} // namespace AMP::TimeIntegrator
