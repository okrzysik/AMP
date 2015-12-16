#ifndef included_ImplicitTimeIntegrator_C
#define included_ImplicitTimeIntegrator_C

#ifndef included_ImplicitTimeIntegrator
#include "ImplicitTimeIntegrator.h"
#endif

#ifndef included_Impl0icitTimeIntegratorParameters
#include "ImplicitTimeIntegratorParameters.h"
#endif

#include "TimeOperatorParameters.h"

/*Design-By-Contract Macros*/
#include "utils/Utilities.h"

#include "utils/AMP_MPI.h"
#include "utils/Utilities.h"
#include <fstream>
#include <stdio.h>
#include <stdlib.h>

#ifndef NULL
#define NULL ( 0 )
#endif

namespace AMP {
namespace TimeIntegrator {

/*
*************************************************************************
*                                                                       *
* Constructor and destructor for ImplicitTimeIntegrator.  The         *
* constructor sets default values for data members, then overrides      *
* them with values read from input or restart.  The destructor does     *
* nothing interesting.                                                  *
*                                                                       *
*************************************************************************
*/

ImplicitTimeIntegrator::ImplicitTimeIntegrator(
    AMP::shared_ptr<TimeIntegratorParameters> parameters )
    : TimeIntegrator( parameters )
{
    initialize( parameters );
}

ImplicitTimeIntegrator::~ImplicitTimeIntegrator() {}

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

void ImplicitTimeIntegrator::initialize( AMP::shared_ptr<TimeIntegratorParameters> parameters )
{
    AMP::shared_ptr<ImplicitTimeIntegratorParameters> params =
        AMP::dynamic_pointer_cast<ImplicitTimeIntegratorParameters>( parameters );

    if ( params.get() != NULL ) {
        d_solver = params->d_solver;

        /*
        * Initialize object with data read from input and restart databases.
        */
        getFromInput( params->d_db );
    } else {
        AMP_ERROR( "ImplicitTimeIntegrator::ImplicitTimeIntegrator: TimeIntegratorParameters "
                   "argument must be of "
                   "derived type ImplicitTimeIntegratorParameters" );
    }

    initializeTimeOperator( parameters );

    d_solver->registerOperator( d_operator );
}

/*
*************************************************************************
*                                                                       *
* Integrate solution through given time increment:                      *
*                                                                       *
* (1) Construct initial guess at new solution by extrapolation.         *
*                                                                       *
* (2) Call the equation advance set up routine.                         *
*                                                                       *
* (3) Compute the new solution using the nonlinear solver.              *
*                                                                       *
* (4) Return integer return code define by nonlinear solver.            *
*                                                                       *
*************************************************************************
*/

int ImplicitTimeIntegrator::advanceSolution( const double dt, const bool first_step )
{
    int retcode = -1;

    if ( stepsRemaining() && ( d_current_time < d_final_time ) ) {

        d_current_dt = dt;

        d_pTimeOperatorParameters->d_db->putDouble( "CurrentDt", dt );

        d_operator->reset( d_pTimeOperatorParameters );

        setInitialGuess( first_step, d_current_time, d_current_dt, d_old_dt );

        AMP::shared_ptr<AMP::LinearAlgebra::Vector> rhs;
        rhs.reset();

        d_solver->setInitialGuess( d_solution );
        d_solver->solve( rhs, d_solution );
    }

    return ( retcode );
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

void ImplicitTimeIntegrator::getFromInput( AMP::shared_ptr<AMP::Database> db )
{
    AMP_ASSERT( db.get() != NULL );
}

/*
*************************************************************************
*                                                                       *
* Write out class version number and data members to database.          *
*                                                                       *
*************************************************************************
*/

void ImplicitTimeIntegrator::putToDatabase( AMP::shared_ptr<AMP::Database> db )
{
    TimeIntegrator::putToDatabase( db );
}

/*
*************************************************************************
*                                                                       *
* Check to make sure that the version number of the class is that same  *
* as the version number in the restart file.  If these values are equal *
* then read values for data members from the restart file.              *
*                                                                       *
*************************************************************************
*/

void ImplicitTimeIntegrator::getFromRestart() {}

/*
*************************************************************************
*                                                                       *
* Print class data members to given output stream.                      *
*                                                                       *
*************************************************************************
*/

void ImplicitTimeIntegrator::printClassData( std::ostream &os ) const
{
    TimeIntegrator::printClassData( os );
}

void ImplicitTimeIntegrator::initializeTimeOperator( AMP::shared_ptr<TimeIntegratorParameters> ) {}
}
}

#endif
