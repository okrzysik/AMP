#include "utils/Utilities.h"

#ifndef included_TimeIntegratorParameters
#include "TimeIntegratorParameters.h"
#endif
#ifndef included_ExplicitEuler
#include "ExplicitEuler.h"
#endif

/*Design-By-Contract Macros*/
#include "utils/Utilities.h"

namespace AMP {
namespace TimeIntegrator {

/*
************************************************************************
*                                                                      *
*  Constructor.                                                        *
*                                                                      *
************************************************************************
*/
ExplicitEuler::ExplicitEuler( AMP::shared_ptr<TimeIntegratorParameters> parameters )
    : TimeIntegrator( parameters )
{
    initialize( parameters );
}

/*
************************************************************************
*                                                                      *
*  Destructor.                                                         *
*                                                                      *
************************************************************************
*/
ExplicitEuler::~ExplicitEuler() {}

/*
************************************************************************
*                                                                      *
* Initialize.                                                          *
*                                                                      *
************************************************************************
*/
void ExplicitEuler::initialize( AMP::shared_ptr<TimeIntegratorParameters> parameters )
{
    AMP_ASSERT( parameters.get() != nullptr );

    TimeIntegrator::initialize( parameters );

    setupVectors();

    /*
     * Initialize data members from input.
     */
    getFromInput( parameters->d_db );
}

void ExplicitEuler::reset( AMP::shared_ptr<TimeIntegratorParameters> parameters )
{
    AMP_ASSERT( parameters.get() != nullptr );

    AMP_ERROR( "Not Finished" );
}

void ExplicitEuler::setupVectors( void )
{

    // clone vectors so they have the same data layout as d_solution
    d_new_solution = d_solution->cloneVector( "new solution" );
    d_f_vec        = d_solution->cloneVector( "f term" );

    /* allocateVectorData is no longer necessary
    d_new_solution->allocateVectorData();
    d_f_vec->allocateVectorData();
    */

    /*
     * Set initial value of vectors to 0.
     */
    d_new_solution->setToScalar( (double) 0.0 );
    d_f_vec->setToScalar( (double) 0.0 );
}

int ExplicitEuler::advanceSolution( const double dt, const bool first_step )
{
    AMP::shared_ptr<AMP::LinearAlgebra::Vector> f;

    if ( first_step ) {
        d_current_dt = d_initial_dt;
    } else {
        d_current_dt = dt;
    }

    if ( stepsRemaining() && ( d_current_time < d_final_time ) ) {
        // f_vec = f(tn,un)
        d_operator->apply( d_solution, d_f_vec );
        // u* = un+dt*f
        d_new_solution->axpy( d_current_dt, *d_f_vec, *d_solution );
    }

    return ( 0 );
}

/*
************************************************************************
*                                                                      *
*  Check whether time advanced solution is acceptable.                 *
*                                                                      *
************************************************************************
*/
bool ExplicitEuler::checkNewSolution( void ) const
{
    /*
     * Ordinarily we would check the actual error in the solution
     * (proportional to the size of d_corrector) against a specified
     * tolerance.  For now, accept everything.
     */
    return ( true );
}

/*
************************************************************************
*                                                                      *
*  Update internal state to reflect time advanced solution.            *
*                                                                      *
************************************************************************
*/
void ExplicitEuler::updateSolution( void ) { 
    d_current_time += d_current_dt;
    d_solution->swapVectors( *d_new_solution ); 
}

/*
************************************************************************
*                                                                      *
* Read input from database.                                            *
*                                                                      *
************************************************************************
*/
void ExplicitEuler::getFromInput( AMP::shared_ptr<AMP::Database> input_db )
{
    if ( input_db->keyExists( "initial_timestep" ) ) {
        d_initial_dt = input_db->getDouble( "initial_timestep" );
    } else {
        AMP_ERROR( d_object_name << " -- Key data `initial_timestep'"
                                 << " missing in input." );
    }
}

double ExplicitEuler::getNextDt( const bool ) { return d_current_dt; }
}
}
