//
// $Id: ExplicitEuler.C,v 1.2 2006/02/07 17:36:01 philipb Exp $
// $Revision: 1.2 $
// $Date: 2006/02/07 17:36:01 $
//
// File:  ExplicitEuler.h
// Copyright:  (c) 2005 The Regents of the University of California
// Description:  Concrete time integrator using explicit Euler
//
#include "AMP/time_integrators/ExplicitEuler.h"

#include "AMP/time_integrators/TimeIntegratorParameters.h"
#include "AMP/vectors/Vector.h"

#include "ProfilerApp.h"

namespace AMP::TimeIntegrator {

/*
************************************************************************
*                                                                      *
*  Constructor.                                                        *
*                                                                      *
************************************************************************
*/
ExplicitEuler::ExplicitEuler(
    std::shared_ptr<AMP::TimeIntegrator::TimeIntegratorParameters> parameters )
    : AMP::TimeIntegrator::TimeIntegrator( parameters )
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
void ExplicitEuler::initialize(
    std::shared_ptr<AMP::TimeIntegrator::TimeIntegratorParameters> parameters )
{

    AMP::TimeIntegrator::TimeIntegrator::initialize( parameters );

    setupVectors();

    /*
     * Initialize data members from input.
     */
    getFromInput( parameters->d_db );
}

void ExplicitEuler::reset(
    std::shared_ptr<const AMP::TimeIntegrator::TimeIntegratorParameters> parameters )
{
    AMP_ASSERT( parameters != nullptr );

    abort();
}

void ExplicitEuler::setupVectors()
{

    // clone vectors so they have the same data layout as d_solution_vector
    d_new_solution = d_solution_vector->clone( "new solution" );
    d_f_vec        = d_solution_vector->clone( "f term" );

    /*
     * Set initial value of vectors to 0.
     */
    d_new_solution->setToScalar( (double) 0.0 );
    d_f_vec->setToScalar( (double) 0.0 );
}

int ExplicitEuler::advanceSolution( const double dt,
                                    const bool first_step,
                                    std::shared_ptr<AMP::LinearAlgebra::Vector> in,
                                    std::shared_ptr<AMP::LinearAlgebra::Vector> out )
{
    PROFILE_START( "ExplicitEuler::advanceSolution1" );
    d_solution_vector = in;
    if ( first_step ) {
        d_current_dt = d_initial_dt;
    } else {
        d_current_dt = dt;
    }
    PROFILE_STOP( "ExplicitEuler::advanceSolution1" );
    if ( stepsRemaining() && ( d_current_time < d_final_time ) ) {

        if ( d_iDebugPrintInfoLevel > 0 ) {
            AMP::pout << "ExplicitEuler::advanceSolution "
                      << ", timestep # " << d_integrator_step << ", current_time " << d_current_time
                      << ", dt " << d_current_dt << std::endl;
            AMP::pout << "L2 Norm of d_solution_vector " << d_solution_vector->L2Norm()
                      << std::endl;
        }
        // f_vec = f(tn,un)
        PROFILE_START( "ExplicitEuler::advanceSolution apply" );
        d_operator->apply( d_solution_vector, d_f_vec );
        PROFILE_STOP( "ExplicitEuler::advanceSolution apply" );
        if ( d_iDebugPrintInfoLevel > 0 ) {
            AMP::pout << "L2 Norm of f(d_solution_vector) " << d_f_vec->L2Norm() << std::endl;
        }


        // u* = un+dt*f
        PROFILE_START( "ExplicitEuler::advanceSolution axby" );
        d_new_solution->axpy( d_current_dt, *d_f_vec, *d_solution_vector );
        PROFILE_STOP( "ExplicitEuler::advanceSolution axby" );

        if ( d_iDebugPrintInfoLevel > 0 ) {
            AMP::pout << "L2 Norm of d_new_solution " << d_new_solution->L2Norm() << std::endl;
        }
    }

    out->copyVector( d_new_solution );

    return ( 1 );
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
    bool retcode = true;
    /*
     * Ordinarily we would check the actual error in the solution
     * (proportional to the size of d_corrector) against a specified
     * tolerance.  For now, accept everything.
     */
    if ( ( d_iDebugPrintInfoLevel > 0 ) && retcode ) {
        AMP::pout << "\n++++++++++++++++++++++++++++++++++++++++++++++++" << std::endl;
        AMP::pout << "End of timestep # " << d_integrator_step << std::endl;
        AMP::pout << "Advance solution past " << d_current_time << std::endl;
        AMP::pout << "++++++++++++++++++++++++++++++++++++++++++++++++\n" << std::endl;
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
void ExplicitEuler::updateSolution()
{
    d_solution_vector->swapVectors( d_new_solution );
    d_current_time += d_current_dt;
    ++d_integrator_step;

    if ( d_iDebugPrintInfoLevel > 0 ) {

        AMP::pout << "\n++++++++++++++++++++++++++++++++++++++++++++++++" << std::endl;
        AMP::pout << "End of timestep # " << d_integrator_step - 1 << std::endl;
        AMP::pout << "Simulation time is " << d_current_time << std::endl;
        AMP::pout << "++++++++++++++++++++++++++++++++++++++++++++++++\n" << std::endl;
    }
}

} // namespace AMP::TimeIntegrator
