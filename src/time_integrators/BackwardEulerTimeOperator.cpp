#include "BackwardEulerTimeOperator.h"
#include "AMP/utils/Database.h"
#include "AMP/utils/Utilities.h"

namespace AMP {
namespace TimeIntegrator {

BackwardEulerTimeOperator::BackwardEulerTimeOperator(
    std::shared_ptr<AMP::Operator::OperatorParameters> params )
    : TimeOperator( params )
{
}

void BackwardEulerTimeOperator::apply( AMP::LinearAlgebra::Vector::const_shared_ptr u,
                                       AMP::LinearAlgebra::Vector::shared_ptr r )
{

    // this routine evaluates a*[ ( M(u)-M(uOld) )/dt-fRhs(u) -source_term] +b*f
    // where the time operator is given by u_t = fRhs(u)

    std::shared_ptr<AMP::LinearAlgebra::Vector> fTmp;

    AMP_INSIST( d_pRhsOperator.get() != nullptr,
                "ERROR: "
                "AMP::TimeIntegrator::TimeIntegrator::TimeOperator::"
                "apply, the rhs operator is NULL!" );

    if ( d_pMassOperator.get() != nullptr ) {
        if ( d_bLinearMassOperator ) {
            d_pScratchVector->subtract( *u, *d_pPreviousTimeSolution );
            d_pMassOperator->apply( d_pScratchVector, r );
            r->scale( 1.0 / d_dCurrentDt );
        } else {
            d_pMassOperator->apply( d_pPreviousTimeSolution, d_pScratchVector );
            d_pMassOperator->apply( u, r );
            r->scale( -1.0 );
            r->add( *r, *d_pScratchVector );
            r->scale( 1.0 / d_dCurrentDt );
        }
    } else {
        r->subtract( *u, *d_pPreviousTimeSolution );
        r->scale( 1.0 / d_dCurrentDt );
    }

    d_pRhsOperator->residual( fTmp, u, d_pScratchVector );

    r->add( *r, *d_pScratchVector );

    // subtract any source sink terms and boundary corrections

    r->subtract( *r, *d_pSourceTerm );
}
} // namespace TimeIntegrator
} // namespace AMP
