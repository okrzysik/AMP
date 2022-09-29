//
// $Id: RK12TimeIntegrator.C,v 1.2 2006/02/07 17:36:01 philipb Exp $
// $Revision: 1.2 $
// $Date: 2006/02/07 17:36:01 $
//
// File:  RK12TimeIntegrator.h
// Copyright:  (c) 2005 The Regents of the University of California
// Description:  Concrete time integrator using Runge-Kutta method of order 2 -
// RK12 method
//
#include "AMP/time_integrators/RK12TimeIntegrator.h"


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
RK12TimeIntegrator::RK12TimeIntegrator(
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
RK12TimeIntegrator::~RK12TimeIntegrator() = default;

/*
************************************************************************
*                                                                      *
* Initialize.                                                          *
*                                                                      *
************************************************************************
*/
void RK12TimeIntegrator::initialize(
    std::shared_ptr<AMP::TimeIntegrator::TimeIntegratorParameters> parameters )
{

    AMP::TimeIntegrator::TimeIntegrator::initialize( parameters );

    setupVectors();

    /*
     * Initialize data members from input.
     */
    getFromInput( parameters->d_db );
}

void RK12TimeIntegrator::reset(
    std::shared_ptr<const AMP::TimeIntegrator::TimeIntegratorParameters> )
{
    // AMP_ASSERT(parameters!=nullptr);
    d_new_solution->getVectorData()->reset();
    d_k1_vec->getVectorData()->reset();
    d_k2_vec->getVectorData()->reset();
    d_z_vec->getVectorData()->reset();
}

void RK12TimeIntegrator::setupVectors()
{

    // clone vectors so they have the same data layout as d_solution_vector
    d_new_solution = d_solution_vector->cloneVector();
    d_k1_vec       = d_solution_vector->cloneVector();
    d_k2_vec       = d_solution_vector->cloneVector();
    d_z_vec        = d_solution_vector->cloneVector();

    /*
     * Set initial value of vectors to 0.
     */
    d_new_solution->setToScalar( (double) 0.0 );
    d_k1_vec->setToScalar( (double) 0.0 );
    d_k2_vec->setToScalar( (double) 0.0 );
    d_z_vec->setToScalar( (double) 0.0 );
}

int RK12TimeIntegrator::advanceSolution( const double dt,
                                         const bool,
                                         std::shared_ptr<AMP::LinearAlgebra::Vector> in,
                                         std::shared_ptr<AMP::LinearAlgebra::Vector> out )
{
    PROFILE_START( "advanceSolution" );

    d_solution_vector = in;
    d_current_dt      = dt;

    if ( d_iDebugPrintInfoLevel > 5 ) {
        AMP::pout << "*****************************************" << std::endl;
        AMP::pout << "Current solution vector: " << std::endl;
        d_solution_vector->getVectorData()->print( AMP::pout );
        AMP::pout << "*****************************************" << std::endl;
    }

    // k1 = f(tn,un)
    d_operator->apply( d_solution_vector, d_k1_vec );

    if ( d_iDebugPrintInfoLevel > 5 ) {
        AMP::pout << "*****************************************" << std::endl;
        AMP::pout << "Current K1 vector: " << std::endl;
        d_k1_vec->getVectorData()->print( AMP::pout );
        AMP::pout << "*****************************************" << std::endl;
    }

    if ( d_iDebugPrintInfoLevel > 4 ) {
        std::cout << "L2 norm of k1 " << d_k1_vec->L2Norm().get<double>() << std::endl;
    }
    // u* = un+k1*dt
    d_z_vec->axpy( d_current_dt, *d_k1_vec, *d_solution_vector );

    if ( d_iDebugPrintInfoLevel > 4 ) {
        std::cout << "L2 norm of  (un+k1*dt/2) " << d_new_solution->L2Norm().get<double>()
                  << std::endl;
    }

    // k2 = f(t+dt, u*)
    d_operator->apply( d_z_vec, d_k2_vec );

    if ( d_iDebugPrintInfoLevel > 4 ) {
        std::cout << "L2 norm of k2 " << d_k2_vec->L2Norm().get<double>() << std::endl;
    }

    // u_{n+1} = un+0.5*dt (k1+k2)
    d_new_solution->axpy( 0.5 * d_current_dt, *d_k2_vec, *d_solution_vector );
    d_new_solution->axpy( 0.5 * d_current_dt, *d_k1_vec, *d_new_solution );

    // store the difference in d_z_vec
    d_z_vec->subtract( *d_new_solution, *d_z_vec );

    if ( d_iDebugPrintInfoLevel > 4 ) {
        std::cout << "L2 norm of (z-u*) " << d_z_vec->L2Norm().get<double>() << std::endl;
    }

    out->copyVector( d_new_solution );
    PROFILE_STOP( "advanceSolution" );
    return ( 1 );
}

/*
************************************************************************
*                                                                      *
*  Check whether time advanced solution is acceptable.                 *
*                                                                      *
************************************************************************
*/
bool RK12TimeIntegrator::checkNewSolution()
{
    bool retcode = false;

    auto l2Norm                 = d_z_vec->L2Norm();
    auto l2NormOfEstimatedError = l2Norm.get<double>();

    // we flag the solution as being acceptable if the l2 norm of the error
    // is less than the required tolerance or we are at the minimum time step
    if ( ( l2NormOfEstimatedError < d_atol ) || ( fabs( d_current_dt - d_min_dt ) < 1.0e-10 ) ) {
        retcode = true;
    }

    if ( ( d_iDebugPrintInfoLevel > 0 ) && ( !retcode ) ) {
        AMP::pout << "\n++++++++++++++++++++++++++++++++++++++++++++++++" << std::endl;
        AMP::pout << "End of timestep # " << d_integrator_step << std::endl;
        AMP::pout << "Failed to advance solution past " << d_current_time << std::endl;
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
void RK12TimeIntegrator::updateSolution()
{
    d_solution_vector->copyVector( d_new_solution );
    d_current_time += d_current_dt;
    ++d_integrator_step;

    if ( d_iDebugPrintInfoLevel > 0 ) {

        AMP::pout << "\n++++++++++++++++++++++++++++++++++++++++++++++++" << std::endl;
        AMP::pout << "End of timestep # " << d_integrator_step - 1 << std::endl;
        AMP::pout << "Simulation time is " << d_current_time << std::endl;
        AMP::pout << "++++++++++++++++++++++++++++++++++++++++++++++++\n" << std::endl;
    }
}

/*
************************************************************************
*                                                                      *
* Read input from database.                                            *
*                                                                      *
************************************************************************
*/
void RK12TimeIntegrator::getFromInput( std::shared_ptr<AMP::Database> input_db )
{

    d_safety_factor = input_db->getWithDefault<double>( "safety_factor", 0.9 );

    d_atol = input_db->getWithDefault<double>( "absolute_tolerance", 1.0e-09 );

    d_use_fixed_dt = input_db->getWithDefault<bool>( "use_fixed_dt", false );
}

double RK12TimeIntegrator::getNextDt( const bool good_solution )
{
    double next_dt;

    if ( d_use_fixed_dt ) {
        next_dt = std::min( d_current_dt, d_final_time - d_current_time );
    } else {

        if ( good_solution ) {
            auto l2NormOfEstimatedError = d_z_vec->L2Norm().get<double>();

            next_dt = d_safety_factor * d_current_dt *
                      pow( ( d_atol / l2NormOfEstimatedError ), 1.0 / 2.0 );

            // check to make sure the timestep is not too small or large
            next_dt = std::min( std::max( next_dt, d_min_dt ), d_max_dt );
            // check to make sure we don't step past final time
            next_dt = std::min( next_dt, d_final_time - d_current_time );

            if ( d_iDebugPrintInfoLevel > 0 ) {
                AMP::pout << "Timestep # " << d_integrator_step << ", dt: " << next_dt << std::endl;
                AMP::pout << "++++++++++++++++++++++++++++++++++++++++++++++++\n" << std::endl;
            }
        } else {
            next_dt = d_safety_factor * d_current_dt;
            ++d_total_steprejects;
            if ( d_iDebugPrintInfoLevel > 0 ) {
                AMP::pout << "Failed to advance timestep # " << d_integrator_step
                          << ", new dt: " << next_dt << std::endl;
                AMP::pout << "++++++++++++++++++++++++++++++++++++++++++++++++\n" << std::endl;
            }
        }
    }

    return next_dt;
}

} // namespace AMP::TimeIntegrator
