#ifndef included_ImplicitTimeIntegratorParameters
#define included_ImplicitTimeIntegratorParameters

#ifndef included_AMP_config

#endif

#ifndef included_tbox_Pointer
#include "boost/shared_ptr.hpp"
#endif

#ifndef included_tbox_InputDatabase
#include "utils/InputDatabase.h"
#endif

#ifndef included_SolverStrategy
#include "solvers/SolverStrategy.h"
#endif

#ifndef included_TimeIntegratorParameters
#include "TimeIntegratorParameters.h"
#endif

namespace AMP{
namespace TimeIntegrator{

/*!
  @brief Parameter class for implicit time integrators

  Class ImplicitTimeIntegratorParameters contains the parameters to 
  initialize an implicit time integrator class. It contains a Database
  object and a pointer to a SolverStrategy object.
  
  @param d_solver pointer to SolverStrategy

  @see SolverStrategy
*/
class ImplicitTimeIntegratorParameters: public TimeIntegratorParameters
{
public:
  ImplicitTimeIntegratorParameters( boost::shared_ptr<AMP::Database> db);
  
  ~ImplicitTimeIntegratorParameters();
  
  /**
   * Pointers to implicit equation and solver strategy objects
   * The strategies provide nonlinear equation and solver 
   * routines for treating the nonlinear problem on the hierarchy.
   */
  boost::shared_ptr<AMP::Solver::SolverStrategy> d_solver;
 protected:
  
 private:
  // not implemented
  ImplicitTimeIntegratorParameters();
  ImplicitTimeIntegratorParameters(const ImplicitTimeIntegratorParameters&);
  void operator=(const ImplicitTimeIntegratorParameters&);
};

}
}

#endif

