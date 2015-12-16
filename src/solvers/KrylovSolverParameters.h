#ifndef included_AMP_KrylovSolverParameters
#define included_AMP_KrylovSolverParameters

#include "solvers/SolverStrategy.h"
#include "solvers/SolverStrategyParameters.h"
#include "utils/Database.h"
#include "utils/shared_ptr.h"

namespace AMP {
namespace Solver {

/**
 * Class KrylovSolverParameters provides a uniform mechanism to pass
 * initialization parameters to Krylov solvers. It contains
 * shared pointers to a database object and a preconditioner (which could be NULL).
 * All member variables are public.
 */
class KrylovSolverParameters : public SolverStrategyParameters
{
public:
    KrylovSolverParameters() {}
    explicit KrylovSolverParameters( const AMP::shared_ptr<AMP::Database> db );
    virtual ~KrylovSolverParameters() {}

    AMP_MPI d_comm;

    AMP::shared_ptr<AMP::Solver::SolverStrategy> d_pPreconditioner;

protected:
private:
};
}
}

#endif
