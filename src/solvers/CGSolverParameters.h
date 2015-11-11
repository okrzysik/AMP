#ifndef included_AMP_CGSolverParameters
#define included_AMP_CGSolverParameters

#include "utils/shared_ptr.h"
#include "utils/Database.h"
#include "solvers/SolverStrategyParameters.h"
#include "solvers/SolverStrategy.h"

namespace AMP {
namespace Solver {

  /**
   * Class CGSolverParameters provides a uniform mechanism to pass
   * initialization parameters to the CGSolver solver. It contains
   * shared pointers to a database object and a preconditioner. All member variables are public.
   */
  class CGSolverParameters: public SolverStrategyParameters{
  public:
    CGSolverParameters(){}
    CGSolverParameters(const AMP::shared_ptr<AMP::Database> db);
    virtual ~CGSolverParameters(){}

    AMP_MPI d_comm;

    AMP::shared_ptr<AMP::Solver::SolverStrategy> d_pPreconditioner;

  protected:
  private:
    
  };

}
}

#endif
