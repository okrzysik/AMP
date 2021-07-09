#include "AMP/time_integrators/RK23TimeIntegrator.h"
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
RK23TimeIntegrator::RK23TimeIntegrator( std::shared_ptr<TimeIntegratorParameters> parameters )
    : TimeIntegrator( parameters )
{
    d_safety_factor = 0.0;
    d_atol          = 0.0;

    initialize( parameters );
}

/*
************************************************************************
*                                                                      *
*  Destructor.                                                         *
*                                                                      *
************************************************************************
*/
RK23TimeIntegrator::~RK23TimeIntegrator() = default;

/*
************************************************************************
*                                                                      *
* Initialize.                                                          *
*                                                                      *
************************************************************************
*/
void RK23TimeIntegrator::initialize( std::shared_ptr<TimeIntegratorParameters> parameters )
{
    AMP_ASSERT( parameters.get() != (TimeIntegratorParameters *) nullptr );

    TimeIntegrator::initialize( parameters );

    setupVectors();

    /*
     * Initialize data members from input.
     */
    getFromInput( parameters->d_db );
}

void RK23TimeIntegrator::reset( std::shared_ptr<const TimeIntegratorParameters> parameters )
{
    AMP_ASSERT( parameters.get() != (TimeIntegratorParameters *) nullptr );

    AMP_ERROR( "Not Finished" );
}

void RK23TimeIntegrator::setupVectors()
{

    // clone vectors so they have the same data layout as d_solution_vector
    d_new_solution = d_solution_vector->cloneVector( "new solution" );
    d_k1_vec       = d_solution_vector->cloneVector( "k1 term" );
    d_k2_vec       = d_solution_vector->cloneVector( "k2 term" );
    d_k3_vec       = d_solution_vector->cloneVector( "k3 term" );
    d_k4_vec       = d_solution_vector->cloneVector( "k4 term" );
    d_z_vec        = d_solution_vector->cloneVector( "z term" );

    /* Allocate vector data no longer necessary
    d_new_solution->allocateVectorData();
    d_k1_vec->allocateVectorData();
    d_k2_vec->allocateVectorData();
    d_k3_vec->allocateVectorData();
    d_k4_vec->allocateVectorData();
    d_z_vec->allocateVectorData();
    */

    /*
     * Set initial value of vectors to 0.
     */
    d_new_solution->setToScalar( (double) 0.0 );
    d_k1_vec->setToScalar( (double) 0.0 );
    d_k2_vec->setToScalar( (double) 0.0 );
    d_k3_vec->setToScalar( (double) 0.0 );
    d_k4_vec->setToScalar( (double) 0.0 );
    d_z_vec->setToScalar( (double) 0.0 );
}

int RK23TimeIntegrator::advanceSolution( const double dt, const bool first_step )
{
    if ( first_step ) {
        // k1 = f(tn,un)
        d_operator->apply( d_solution_vector, d_k1_vec );
    } else {
        d_current_dt = dt;
        d_k1_vec->swapVectors( *d_k4_vec );
    }

    // u* = un+k1*dt/2
    d_new_solution->axpy( d_current_dt / 2.0, *d_k1_vec, *d_solution_vector );
    // k2 = f(t+dt/2, u*)
    d_operator->apply( d_new_solution, d_k2_vec );
    // u* = un+0.75*k2*dt
    d_new_solution->axpy( 0.75 * d_current_dt, *d_k2_vec, *d_solution_vector );
    // k3 = f(t+0.75dt, u*)
    d_operator->apply( d_new_solution, d_k3_vec );

    // first we calculate the 3rd order solution in d_new_solution
    // u* = un+k1*2dt/9
    d_new_solution->axpy( 2.0 * d_current_dt / 9.0, *d_k1_vec, *d_solution_vector );
    // u* = u*+k2*dt/3
    d_new_solution->axpy( d_current_dt / 3.0, *d_k2_vec, *d_new_solution );
    // u* = u*+k3*4dt/9
    d_new_solution->axpy( 4.0 * d_current_dt / 9.0, *d_k3_vec, *d_new_solution );

    // k4 = f(t+dt, u*)
    d_operator->apply( d_new_solution, d_k4_vec );

    // now we calculate the 2nd order solution in d_z_vec for adapting the timestep
    // z = un+ dt*(7k1/24+k2/4+k3/3+k4/8)
    d_z_vec->axpy( 7.0 / 24.0 * d_current_dt, *d_k1_vec, *d_solution_vector );
    d_z_vec->axpy( 0.25 * d_current_dt, *d_k2_vec, *d_z_vec );
    d_z_vec->axpy( 1.0 / 3.0 * d_current_dt, *d_k3_vec, *d_z_vec );
    d_z_vec->axpy( 0.125 * d_current_dt, *d_k4_vec, *d_z_vec );

    // store the difference in d_z_vec
    d_z_vec->subtract( *d_new_solution, *d_z_vec );

    return ( 0 );
}

/*
************************************************************************
*                                                                      *
*  Check whether time advanced solution is acceptable.                 *
*                                                                      *
************************************************************************
*/
bool RK23TimeIntegrator::checkNewSolution()
{
    bool retcode = false;

    double l2NormOfEstimatedError = d_z_vec->L2Norm().get<double>();

    // we flag the solution as being acceptable if the l2 norm of the error
    // is less than the required tolerance or we are at the minimum time step
    if ( ( l2NormOfEstimatedError < d_atol ) || ( fabs( d_current_dt - d_min_dt ) < 1.0e-10 ) ) {
        retcode = true;
    }

    return ( retcode );
}

/*
************************************************************************
*                                                                      *
*  Update internal state to reflect time advanced solution.            *
*                                                                      *
************************************************************************
*/
void RK23TimeIntegrator::updateSolution()
{
    d_solution_vector->swapVectors( *d_new_solution );
    d_current_time += d_current_dt;
}

/*
************************************************************************
*                                                                      *
* Read input from database.                                            *
*                                                                      *
************************************************************************
*/
void RK23TimeIntegrator::getFromInput( std::shared_ptr<AMP::Database> input_db )
{
    if ( input_db->keyExists( "initial_timestep" ) ) {
        d_initial_dt = input_db->getScalar<double>( "initial_timestep" );
        d_current_dt = d_initial_dt;
    } else {
        AMP_ERROR( d_object_name << " -- Key data `initial_timestep'"
                                 << " missing in input." );
    }

    d_safety_factor = input_db->getWithDefault<double>( "safety_factor", 0.9 );

    d_atol = input_db->getWithDefault<double>( "absolute_tolerance", 1.0e-09 );
}

double RK23TimeIntegrator::getNextDt( const bool )
{
    double d_next_dt;

    double l2NormOfEstimatedError( d_z_vec->L2Norm() );

    d_next_dt =
        d_safety_factor * d_current_dt * pow( ( d_atol / l2NormOfEstimatedError ), 1.0 / 3.0 );

    // check to make sure the timestep is not too small
    if ( d_next_dt < d_min_dt )
        d_next_dt = d_min_dt;

    // check to make sure the timestep is not too large
    if ( d_next_dt > d_max_dt )
        d_next_dt = d_max_dt;

    // check to make sure the time step does not result in stepping beyond final time
    if ( d_current_time + d_next_dt > d_final_time )
        d_next_dt = d_final_time - d_current_time;

    return d_next_dt;
}
} // namespace TimeIntegrator
} // namespace AMP
