#ifndef included_ExplicitEuler
#define included_ExplicitEuler

#include <string>

#ifndef included_TimeIntegrator
#include "TimeIntegrator.h"
#endif

namespace AMP{
namespace TimeIntegrator{

/** \class ExplicitEuler
 * 
 * Class ExplicitEuler is a concrete time integrator
 * that implements the explicit Runge-Kutta second order (RK2) method.
 */
class ExplicitEuler : public TimeIntegrator
{
public:
   /**
    * Constructor that accepts parameter list.
    */
   ExplicitEuler(  AMP::shared_ptr<TimeIntegratorParameters> parameters );

   /**
    * Destructor.
    */
   virtual ~ExplicitEuler();

   /**
    * Initialize from parameter list.
    */
   void initialize( AMP::shared_ptr<TimeIntegratorParameters> parameters );

   /**
   * Resets the internal state of the time integrator as needed.
   * A parameter argument is passed to allow for general flexibility
   * in determining what needs to be reset Typically used after a regrid.
   */
   void reset( AMP::shared_ptr<TimeIntegratorParameters> parameters);

   /**
    * Specify initial time step.
    */
   double getInitialDt();

   /**
    * Specify next time step to use.
    */
   double getNextDt( const bool good_solution );
   
   /**
   * Determine whether time advanced solution is satisfactory.
    */
   bool checkNewSolution( void ) const;
   
   /**
   * Update state of the solution. 
   */
   void updateSolution( void );

   int advanceSolution( const double dt, const bool first_step );

private:
   /**
    * Constructor.
    */
   ExplicitEuler();

   /**
    * Read data from input database.
    */
   void getFromInput( AMP::shared_ptr<AMP::Database> input_db );

   /**
   * setup the vectors used by BE
   */
   void setupVectors(void);

   int d_number_regrid_states;

   AMP::shared_ptr<AMP::LinearAlgebra::Vector> d_new_solution;
   AMP::shared_ptr<AMP::LinearAlgebra::Vector> d_f_vec;

};

}
}

#endif
