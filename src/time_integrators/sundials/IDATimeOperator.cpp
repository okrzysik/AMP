#include "AMP/time_integrators/sundials/IDATimeOperator.h"


namespace AMP {
namespace TimeIntegrator {

IDATimeOperator::IDATimeOperator( std::shared_ptr<AMP::Operator::OperatorParameters> in_params )
    : TimeOperator( in_params ), d_current_time( 0 )
{
    d_cloningHappened = false;
    // BP, commenting out because this class does not define
    // a reset and the reset called is the reset for the base
    // class which already has been called once
    //        reset(in_params);

    // JL
    d_beta = 1.0;
}

IDATimeOperator::~IDATimeOperator() = default;


/*
void
IDATimeOperator::reset( std::shared_ptr<OperatorParameters> in_params)
{
    auto params = std::dynamic_pointer_cast<IDATimeOperatorParameters>(in_params);
    getFromInput(params->d_db);
}
*/

void IDATimeOperator::apply( AMP::LinearAlgebra::Vector::const_shared_ptr u,
                             AMP::LinearAlgebra::Vector::shared_ptr r )
{
    if ( d_cloningHappened == 0 ) {
        d_pScratchVector = r->cloneVector();
        d_pScratchVector->zero();
        d_cloningHappened = true;
    }

    d_pMassOperator->apply( d_pIDATimeDerivative, d_pScratchVector );

    if ( d_iDebugPrintInfoLevel > 4 ) {
        AMP::pout << "Output of M * yp in IDATimeOperator" << std::endl;
        AMP::pout << d_pScratchVector << std::endl;
    }

    if ( d_pAlgebraicVariable ) {
        auto algebraicComponent = d_pScratchVector->subsetVectorForVariable( d_pAlgebraicVariable );
        algebraicComponent->zero();
    }

    d_pRhsOperator->apply( u, r );
    r->axpby( 1.0, 1.0, *d_pScratchVector );

    bool dpSourceTermNull = ( d_pSourceTerm.get() == nullptr );

    if ( d_iDebugPrintInfoLevel > 5 ) {
        AMP::pout << "Output of M * yp-frhs(y,t) in IDATimeOperator" << std::endl;
        AMP::pout << r << std::endl;
    }

    if ( !dpSourceTermNull ) {
        r->axpby( -1.0, 1.0, *d_pSourceTerm );
    }

    if ( d_iDebugPrintInfoLevel > 6 ) {
        AMP::pout << "Output of M * yp-frhs(y,t)-g in IDATimeOperator" << std::endl;
        AMP::pout << r << std::endl;
    }
}
} // namespace TimeIntegrator
} // namespace AMP
