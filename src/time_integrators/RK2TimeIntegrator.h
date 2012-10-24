#ifndef included_RK2TimeIntegrator
#define included_RK2TimeIntegrator

#include <string>

#ifndef included_TimeIntegrator
#include "TimeIntegrator.h"
#endif

namespace AMP{
namespace TimeIntegrator{

/** \class RK2TimeIntegrator
 * 
 * Class RK2TimeIntegrator is a concrete time integrator
 * that implements the explicit Runge-Kutta second order (RK2) method.
 */
class RK2TimeIntegrator : public TimeIntegrator
{
public:
   /**
    * Constructor that accepts parameter list.
    */
   RK2TimeIntegrator( boost::shared_ptr<TimeIntegratorParameters> parameters );

   /**
    * Destructor.
    */
   virtual ~RK2TimeIntegrator();

   /**
    * Initialize from parameter list.
    */
   void initialize( boost::shared_ptr<TimeIntegratorParameters> parameters );

   /**
   * Resets the internal state of the time integrator as needed.
   * A parameter argument is passed to allow for general flexibility
   * in determining what needs to be reset Typically used after a regrid.
   */
   void reset( boost::shared_ptr<TimeIntegratorParameters> parameters);

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
   RK2TimeIntegrator();

   /**
    * Read data from input database.
    */
   void getFromInput( boost::shared_ptr<AMP::Database> input_db );

   /**
   * setup the vectors used by RK2
   */
   void setupVectors(void);

   int d_number_regrid_states;

   double d_initial_dt;

   boost::shared_ptr<AMP::LinearAlgebra::Vector> d_new_solution;
   boost::shared_ptr<AMP::LinearAlgebra::Vector> d_k1_vec;
   boost::shared_ptr<AMP::LinearAlgebra::Vector> d_k2_vec;

};

}
}

#endif
