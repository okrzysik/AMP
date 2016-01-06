#ifndef included_ImplicitTimeIntegratorParameters
#define included_ImplicitTimeIntegratorParameters

#include "TimeIntegratorParameters.h"
#include "solvers/SolverStrategy.h"
#include "utils/InputDatabase.h"
#include "utils/shared_ptr.h"


namespace AMP {
namespace TimeIntegrator {

/*!
  @brief Parameter class for implicit time integrators

  Class ImplicitTimeIntegratorParameters contains the parameters to
  initialize an implicit time integrator class. It contains a Database
  object and a pointer to a SolverStrategy object.

  @param d_solver pointer to SolverStrategy

  @see SolverStrategy
*/
class ImplicitTimeIntegratorParameters : public TimeIntegratorParameters
{
public:
    explicit ImplicitTimeIntegratorParameters( AMP::shared_ptr<AMP::Database> db );

    virtual ~ImplicitTimeIntegratorParameters();

    /**
     * Pointers to implicit equation and solver strategy objects
     * The strategies provide nonlinear equation and solver
     * routines for treating the nonlinear problem on the hierarchy.
     */
    AMP::shared_ptr<AMP::Solver::SolverStrategy> d_solver;

protected:
private:
    // not implemented
    ImplicitTimeIntegratorParameters();
    explicit ImplicitTimeIntegratorParameters( const ImplicitTimeIntegratorParameters & );
    void operator=( const ImplicitTimeIntegratorParameters & );
};
}
}

#endif
