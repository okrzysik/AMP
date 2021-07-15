#include "AMP/time_integrators/RK2TimeIntegrator.h"
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
RK2TimeIntegrator::RK2TimeIntegrator( std::shared_ptr<TimeIntegratorParameters> parameters )
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
RK2TimeIntegrator::~RK2TimeIntegrator() = default;

/*
************************************************************************
*                                                                      *
* Initialize.                                                          *
*                                                                      *
************************************************************************
*/
void RK2TimeIntegrator::initialize( std::shared_ptr<TimeIntegratorParameters> parameters )
{
    AMP_ASSERT( parameters.get() != (TimeIntegratorParameters *) nullptr );

    TimeIntegrator::initialize( parameters );

    setupVectors();

    /*
     * Initialize data members from input.
     */
    getFromInput( parameters->d_db );
}

void RK2TimeIntegrator::reset( std::shared_ptr<const TimeIntegratorParameters> parameters )
{
    AMP_ASSERT( parameters.get() != (TimeIntegratorParameters *) nullptr );

    AMP_ERROR( "Not Finished" );
}

void RK2TimeIntegrator::setupVectors()
{

    // clone vectors so they have the same data layout as d_solution_vector
    d_new_solution = d_solution_vector->cloneVector( "new solution" );
    d_k1_vec       = d_solution_vector->cloneVector( "k1 term" );
    d_k2_vec       = d_solution_vector->cloneVector( "k2 term" );

    /* allocateVectorData no longer necessary
    d_new_solution->allocateVectorData();
    d_k1_vec->allocateVectorData();
    d_k2_vec->allocateVectorData();
    */

    /*
     * Set initial value of vectors to 0.
     */
    d_new_solution->setToScalar( (double) 0.0 );
    d_k1_vec->setToScalar( (double) 0.0 );
    d_k2_vec->setToScalar( (double) 0.0 );
}

int RK2TimeIntegrator::advanceSolution( const double dt,
                                        const bool,
                                        std::shared_ptr<AMP::LinearAlgebra::Vector> in,
                                        std::shared_ptr<AMP::LinearAlgebra::Vector> out )
{
    d_solution_vector = in;
    // k1 = f(tn,un)
    d_operator->apply( d_solution_vector, d_k1_vec );
    // u* = un+dt*k1
    d_new_solution->axpy( dt, *d_k1_vec, *d_solution_vector );
    // k2 = f(t+dt, u*)
    d_operator->apply( d_new_solution, d_k2_vec );
    // u_new = un+ dt*(k1+k2)/2
    d_k2_vec->add( *d_k1_vec, *d_k2_vec );
    d_new_solution->axpy( dt / 2.0, *d_k2_vec, *d_solution_vector );
    out->copyVector( d_new_solution );
    return ( 0 );
}

/*
************************************************************************
*                                                                      *
*  Check whether time advanced solution is acceptable.                 *
*                                                                      *
************************************************************************
*/
bool RK2TimeIntegrator::checkNewSolution()
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
void RK2TimeIntegrator::updateSolution()
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
void RK2TimeIntegrator::getFromInput( std::shared_ptr<AMP::Database> input_db )
{
    if ( input_db->keyExists( "initial_timestep" ) ) {
        d_initial_dt = input_db->getScalar<double>( "initial_timestep" );
    } else {
        AMP_ERROR( d_object_name << " -- Key data `initial_timestep'"
                                 << " missing in input." );
    }
}

double RK2TimeIntegrator::getNextDt( const bool ) { return d_current_dt; }
} // namespace TimeIntegrator
} // namespace AMP
