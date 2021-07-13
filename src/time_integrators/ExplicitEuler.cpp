#include "AMP/time_integrators/ExplicitEuler.h"
#include "AMP/time_integrators/TimeIntegratorParameters.h"
#include "AMP/utils/Utilities.h"

namespace AMP {
namespace TimeIntegrator {

/*
************************************************************************
*                                                                      *
*  Constructor.                                                        *
*                                                                      *
************************************************************************
*/
ExplicitEuler::ExplicitEuler( std::shared_ptr<TimeIntegratorParameters> parameters )
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
ExplicitEuler::~ExplicitEuler() = default;

/*
************************************************************************
*                                                                      *
* Initialize.                                                          *
*                                                                      *
************************************************************************
*/
void ExplicitEuler::initialize( std::shared_ptr<TimeIntegratorParameters> parameters )
{
    AMP_ASSERT( parameters );

    TimeIntegrator::initialize( parameters );

    setupVectors();

    /*
     * Initialize data members from input.
     */
    getFromInput( parameters->d_db );
}

void ExplicitEuler::reset( std::shared_ptr<const TimeIntegratorParameters> parameters )
{
    AMP_ASSERT( parameters );

    AMP_ERROR( "Not Finished" );
}

void ExplicitEuler::setupVectors()
{

    // clone vectors so they have the same data layout as d_solution_vector
    d_new_solution = d_solution_vector->cloneVector( "new solution" );
    d_f_vec        = d_solution_vector->cloneVector( "f term" );

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
    std::shared_ptr<AMP::LinearAlgebra::Vector> f;

    if ( first_step ) {
        d_current_dt = d_initial_dt;
    } else {
        d_current_dt = dt;
    }

    if ( stepsRemaining() && ( d_current_time < d_final_time ) ) {
        // f_vec = f(tn,un)
        d_operator->apply( d_solution_vector, d_f_vec );
        // u* = un+dt*f
        d_new_solution->axpy( d_current_dt, *d_f_vec, *d_solution_vector );
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
bool ExplicitEuler::checkNewSolution()
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
void ExplicitEuler::updateSolution()
{
    d_current_time += d_current_dt;
    d_solution_vector->swapVectors( *d_new_solution );
}

/*
************************************************************************
*                                                                      *
* Read input from database.                                            *
*                                                                      *
************************************************************************
*/
void ExplicitEuler::getFromInput( std::shared_ptr<AMP::Database> input_db )
{
    if ( input_db->keyExists( "initial_timestep" ) ) {
        d_initial_dt = input_db->getScalar<double>( "initial_timestep" );
    } else {
        AMP_ERROR( d_object_name << " -- Key data `initial_timestep'"
                                 << " missing in input." );
    }
}

double ExplicitEuler::getNextDt( const bool ) { return d_current_dt; }
} // namespace TimeIntegrator
} // namespace AMP
