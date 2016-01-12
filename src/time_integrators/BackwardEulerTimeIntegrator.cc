#include "utils/Utilities.h"

#ifndef included_TimeIntegratorParameters
#include "TimeIntegratorParameters.h"
#endif
#ifndef included_BackwardEulerTimeIntegrator
#include "BackwardEulerTimeIntegrator.h"
#endif

#include "BackwardEulerTimeOperator.h"
#include "TimeOperatorParameters.h"

/*Design-By-Contract Macros*/
#include "utils/Utilities.h"

namespace AMP {
namespace TimeIntegrator {


/***********************************************************************
*  Constructor.                                                        *
***********************************************************************/
BackwardEulerTimeIntegrator::BackwardEulerTimeIntegrator(
    AMP::shared_ptr<TimeIntegratorParameters> parameters )
    : ImplicitTimeIntegrator( parameters )
{
    initialize( parameters );
}


/***********************************************************************
*  Destructor.                                                         *
***********************************************************************/
BackwardEulerTimeIntegrator::~BackwardEulerTimeIntegrator() {}


/***********************************************************************
* Initialize.                                                          *
***********************************************************************/
void BackwardEulerTimeIntegrator::initialize( AMP::shared_ptr<TimeIntegratorParameters> parameters )
{
    AMP_ASSERT( parameters.get() != nullptr );

    /*
     * Initialize data members from input.
     */
    getFromInput( parameters->d_db );

    ImplicitTimeIntegrator::initialize( parameters );

    // This call must take place in the constructor
    initializeTimeOperator( parameters );

    d_solver->registerOperator( d_operator );
}
void BackwardEulerTimeIntegrator::reset( AMP::shared_ptr<TimeIntegratorParameters> parameters )
{
    AMP_ASSERT( parameters.get() != nullptr );

    AMP_ERROR( "Not Finished" );
}


/***********************************************************************
*  Calculate an approximate time advanced solution.  We use FE as the  *
*  predictor.                                                          *
***********************************************************************/
void BackwardEulerTimeIntegrator::setInitialGuess( const bool,
                                                   const double,
                                                   const double,
                                                   const double )
{
    // lousy initial guess - just to get things moving...
    d_solution->setToScalar( (double) 0.0 );
}


/***********************************************************************
*  Update internal state to reflect time advanced solution.            *
***********************************************************************/
void BackwardEulerTimeIntegrator::updateSolution( void )
{
    // we can figure out a swap later
    d_pPreviousTimeSolution->copyVector( d_solution );
}


/***********************************************************************
*                                                                      *
* Read input from database.                                            *
*                                                                      *
***********************************************************************/
void BackwardEulerTimeIntegrator::getFromInput( AMP::shared_ptr<AMP::Database> input_db )
{
    if ( input_db->keyExists( "initial_timestep" ) ) {
        d_initial_dt = input_db->getDouble( "initial_timestep" );
    } else {
        AMP_ERROR( d_object_name << " -- Key data `initial_timestep'"
                                 << " missing in input." );
    }
}
double BackwardEulerTimeIntegrator::getNextDt( const bool ) { return d_current_dt; }


/***********************************************************************
*  Check whether time advanced solution is acceptable.                 *
***********************************************************************/
bool BackwardEulerTimeIntegrator::checkNewSolution( void ) const
{
    /*
     * Ordinarily we would check the actual error in the solution
     * (proportional to the size of d_corrector) against a specified
     * tolerance.  For now, accept everything.
     */
    return ( true );
}
void BackwardEulerTimeIntegrator::initializeTimeOperator(
    AMP::shared_ptr<TimeIntegratorParameters> parameters )
{
    d_pTimeOperatorParameters.reset( new TimeOperatorParameters( parameters->d_db ) );

    d_pTimeOperatorParameters->d_pRhsOperator = parameters->d_operator;

    d_pTimeOperatorParameters->d_pMassOperator = parameters->d_pMassOperator;

    d_pTimeOperatorParameters->d_pPreviousTimeSolution = d_pPreviousTimeSolution;

    d_pTimeOperatorParameters->d_pSourceTerm = d_pSourceTerm;

    // note that we are resetting the d_operator pointer which
    // may initially have pointer to the rhs operator, now replacing
    // it with a pointer to a TimeOperator
    d_operator.reset( new BackwardEulerTimeOperator( d_pTimeOperatorParameters ) );
}
}
}
