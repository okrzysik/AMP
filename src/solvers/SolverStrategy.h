#ifndef included_SolverStrategy
#define included_SolverStrategy

#include "boost/shared_ptr.hpp"
#include "operators/Operator.h"
#include "SolverStrategyParameters.h"
#include "vectors/Vector.h"


namespace AMP {
namespace Solver {

/**
 * Class SolverStrategy is a base class for methods to solve 
 * equations of the form \f$A(u) = f\f$. $A$ may be a nonlinear
 * or linear operator.
 */

class SolverStrategy
{
 public:

  /**
   * Default constructor
   */
   SolverStrategy();

   /**
    * Main constructor for the base class.
      @param [in] parameters The parameters object contains a database object which must contain the
      following fields:
      
     1. type: integer, name: max_iterations (required)
     acceptable values (non-negative integer values)
     
     2. type: double, name: max_error, (required)
     acceptable values (non-negative real values)
     
     3. type: integer, name: print_info_level, default value: 0, 
     acceptable values (non negative integer values, the higher the value the more verbose the debugging information provided)
     
     4. type: bool, name: zero_initial_guess, default value: false, 
        acceptable values (TRUE, FALSE)
    */
   SolverStrategy(boost::shared_ptr<SolverStrategyParameters> parameters);

   /**
    * Default destructor. Currently does not do anything.
    */
   virtual ~SolverStrategy();
   
   /**
    * Solve the system \f$A(u) = f\f$.  This is a pure virtual function that the derived classes
    * need to provide an implementation of.
    @param [in] f : shared pointer to right hand side vector
    @param [out] u : shared pointer to approximate computed solution 
    */
   virtual void solve(boost::shared_ptr<AMP::LinearAlgebra::Vector>  f,
		      boost::shared_ptr<AMP::LinearAlgebra::Vector>  u) = 0;

   /**
    * Initialize the solution vector and potentially create internal vectors needed for solution
    @param [in] parameters The parameters object
    contains a database object. Currently there are no required fields for the database object.
    */
   virtual void initialize(boost::shared_ptr<SolverStrategyParameters> const parameters);

   /**
    * Provide the initial guess for the solver. This is a pure virtual function that the derived classes
    * need to provide an implementation of.
    * @param [in] initialGuess: shared pointer to the initial guess vector.
    */
   virtual void setInitialGuess( boost::shared_ptr<AMP::LinearAlgebra::Vector> initialGuess ) { }
   
   /**
    * Specify stopping criteria.
    @param [in] max_iterations : maximum number of iterations
    @param [in] max_error : error tolerance (l2 error)
    */
   virtual void setConvergenceTolerance( const int max_iterations,
					 const double max_error);
   /**
    * Specify level of diagnostic information printed during iterations.
    @param [in] print_level : integer level value with permissible values 0 and higher. Setting
    to zero should provide minimial debugging information with higher values resulting in increasingly
    verbose information being printed out.
    */
   virtual void setDebugPrintInfoLevel(int print_level){d_iDebugPrintInfoLevel = print_level;}

   /**
    * Get level of diagnostic information printed during iterations.
    */
    int getDebugPrintInfoLevel( void ) { return d_iDebugPrintInfoLevel; }

   /**
    * Return the number of iterations taken by the solver to converge.
    */
   virtual int getIterations(void) const { return(d_iNumberIterations);}

   /**
    * Tells the solver to use an initial guess of zero and not try to
    * copy an initial guess into the solution vector
    @param [in] use_zero_guess : boolean to specify whether zero initial guess should
    be used or not.
    */
   virtual void setZeroInitialGuess(bool use_zero_guess){d_bUseZeroInitialGuess = use_zero_guess;}

   /**
    * Register the operator that the solver will use during solves
    @param [in] op shared pointer to operator \f$A()\f$ for equation \f$A(u) = f\f$ 
   */
   virtual void registerOperator(const boost::shared_ptr<AMP::Operator::Operator> op){d_pOperator = op;}

   /**
   * Resets the operator registered with the solver with new parameters if necessary
   * @param parameters
   *        OperatorParameters object that is NULL by default
   */
   virtual void resetOperator(const boost::shared_ptr<AMP::Operator::OperatorParameters> parameters);

   /**
   * Resets the solver internally with new parameters if necessary
   * @param parameters
   *        SolverStrategyParameters object that is NULL by default
   */
   virtual void reset(boost::shared_ptr<SolverStrategyParameters> parameters){ }

   /**
    * Return a shared pointer to the operator registered with the solver.
    */
   virtual boost::shared_ptr<AMP::Operator::Operator> getOperator(void) { return d_pOperator; }

protected:
   
   void getFromInput(const boost::shared_ptr<AMP::Database>& db);

   int d_iNumberIterations;             // iterations in solver
   
   double d_dResidualNorm;

   int d_iMaxIterations;
   
   double d_dMaxRhs;

   double d_dMaxError;

   int d_iDebugPrintInfoLevel;

   bool d_bUseZeroInitialGuess;

   int d_iObjectId;

   static int d_iInstanceId;       // used to differentiate between different instances of the class

   boost::shared_ptr<AMP::Operator::Operator> d_pOperator;

private:

};

}
}

#endif


