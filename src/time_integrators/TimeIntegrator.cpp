#include "AMP/time_integrators/TimeIntegrator.h"
#include "AMP/operators/Operator.h"
#include "AMP/time_integrators/TimeIntegratorParameters.h"
#include "AMP/utils/AMP_MPI.h"
#include "AMP/utils/Utilities.h"

#include <cstdio>
#include <cstdlib>
#include <fstream>


namespace AMP::TimeIntegrator {

/************************************************************************
 *                                                                       *
 * Constructor and destructor for TimeIntegrator.  The                   *
 * constructor sets default values for data members, then overrides      *
 * them with values read from input or restart.  The destructor does     *
 * nothing interesting.                                      s            *
 *                                                                       *
 ************************************************************************/

TimeIntegrator::TimeIntegrator(
    std::shared_ptr<AMP::TimeIntegrator::TimeIntegratorParameters> parameters )
{
    AMP_INSIST( parameters, "Null parameter" );

    initialize( parameters );
}

TimeIntegrator::~TimeIntegrator() = default;

/*
*************************************************************************
*                                                                       *
* Initialize integrator and nonlinear solver:                           *
*                                                                       *
* (1) Create vector containing solution state advanced in time.         *
*                                                                       *
* (2) Equation class registers data components with solution vector.    *
*                                                                       *
* (3) Initialize nonlinear solver.                                      *
*                                                                       *
*************************************************************************
*/

void TimeIntegrator::initialize( std::shared_ptr<TimeIntegratorParameters> parameters )
{
    d_object_name = parameters->d_object_name;

    // for now the solution is set to the initial conditions by Jungho
    d_ic_vector = parameters->d_ic_vector;
    AMP_ASSERT( d_ic_vector );

    // for now the solution is set to the initial conditions
    //    d_solution_vector = d_ic_vector->clone( "current solution" );
    d_solution_vector = d_ic_vector->clone();
    d_solution_vector->copyVector( d_ic_vector );

    d_pSourceTerm = parameters->d_pSourceTerm;

    if ( d_solution_vector.get() == nullptr ) {
        AMP_ERROR( "TimeIntegrator::TimeIntegrator()::TimeIntegrators must be initialized with non "
                   "null initial "
                   "condition_vectors" );
    }

    // initialize the mass operator
    d_pMassOperator = parameters->d_pMassOperator;
    // initialize the rhs operator
    d_operator = parameters->d_operator;

    getFromInput( parameters->d_db );

    d_current_time = d_initial_time;
}

/*
*************************************************************************
*                                                                       *
* Get next dt from user-supplied equation class.  Timestep selection    *
* is generally based on whether the nonlinear solution iteration        *
* converged and, if so, whether the solution meets some user-defined    *
* criteria.  It is assumed that, before this routine is called, the     *
* routine checkNewSolution() has been called.  The boolean argument     *
* is the return value from that call.  The integer argument is          *
* that which is returned by the particular nonlinear solver package     *
* that generated the solution.                                          *
*                                                                       *
*************************************************************************
*/

double TimeIntegrator::getNextDt( const bool )
{
    return std::min( std::min( d_current_dt, d_max_dt ), d_final_time - d_current_time );
}

/*
*************************************************************************
*                                                                       *
* If simulation is not from restart, read data from input database.     *
* Otherwise, override restart values for a subset of the data members   *
* with those found in input.                                            *
*                                                                       *
*************************************************************************
*/

void TimeIntegrator::getFromInput( std::shared_ptr<const AMP::Database> db )
{
    AMP_ASSERT( db );

    if ( db->keyExists( "initial_time" ) ) {
        d_initial_time = db->getScalar<double>( "initial_time" );
    } else {
        AMP_ERROR( d_object_name + " -- Key data `initial_time' missing in input" );
    }

    if ( db->keyExists( "final_time" ) ) {
        d_final_time = db->getScalar<double>( "final_time" );
        if ( d_final_time < d_initial_time ) {
            AMP_ERROR( d_object_name + " -- Error in input data  final_time < initial_time" );
        }
    } else {
        AMP_ERROR( d_object_name + " -- Key data `final_time' missing in input." );
    }

    if ( db->keyExists( "max_integrator_steps" ) ) {
        d_max_integrator_steps = db->getScalar<int>( "max_integrator_steps" );
        if ( d_max_integrator_steps < 0 ) {
            AMP_ERROR( d_object_name + " -- Error in input data max_integrator_steps < 0" );
        }
    } else {
        AMP_ERROR( d_object_name + " -- Key data `max_integrator_steps' missing in input" );
    }

    d_max_dt = db->getWithDefault<double>( "max_dt", std::numeric_limits<double>::max() );

    if ( d_max_dt < 0.0 ) {
        AMP_ERROR( d_object_name + " -- Error in input data max_dt < 0." );
    }

    d_min_dt = db->getWithDefault<double>( "min_dt", std::numeric_limits<double>::min() );

    if ( d_min_dt < 0.0 ) {
        AMP_ERROR( d_object_name + " -- Error in input data min_dt < 0." );
    }

    if ( db->keyExists( "initial_dt" ) ) {
        d_initial_dt = db->getWithDefault<double>( "initial_dt", 0.0 );
    }

    d_iDebugPrintInfoLevel = db->getWithDefault<int>( "print_info_level", 0 );

    if ( db->keyExists( "name" ) ) {
        d_object_name = db->getString( "name" );
    } else {
        AMP_ERROR( " -- Key data `name' missing in input." );
    }

    d_current_dt = d_initial_dt;
    d_old_dt     = d_initial_dt;
}

/*
*************************************************************************
*                                                                       *
* Write out class version number and data members to database.          *
*                                                                       *
*************************************************************************
*/

void TimeIntegrator::putToDatabase( std::shared_ptr<AMP::Database> db )
{
    AMP_ASSERT( !db.use_count() );

    db->putScalar( "d_initial_time", d_initial_time );
    db->putScalar( "d_final_time", d_final_time );
    db->putScalar( "d_current_time", d_current_time );
    db->putScalar( "d_current_dt", d_current_dt );
    db->putScalar( "d_old_dt", d_old_dt );
    db->putScalar( "d_integrator_step", d_integrator_step );
    db->putScalar( "d_max_integrator_steps", d_max_integrator_steps );
}
/*
*************************************************************************
*                                                                       *
* Print class data members to given output stream.                      *
*                                                                       *
*************************************************************************
*/

void TimeIntegrator::printClassData( std::ostream &os ) const
{
    os << "\nTimeIntegrator::printClassData..." << std::endl;
    os << "TimeIntegrator: this = " << this << std::endl;
    os << "d_object_name = " << d_object_name << std::endl;
    os << "d_initial_time = " << d_initial_time << std::endl;
    os << "d_final_time = " << d_final_time << std::endl;
    os << "d_current_time = " << d_current_time << std::endl;
    os << "d_current_dt = " << d_current_dt << std::endl;
    os << "d_old_dt = " << d_old_dt << std::endl;
    os << "d_integrator_step = " << d_integrator_step << std::endl;
    os << "d_max_integrator_steps = " << d_max_integrator_steps << std::endl;
}
} // namespace AMP::TimeIntegrator
