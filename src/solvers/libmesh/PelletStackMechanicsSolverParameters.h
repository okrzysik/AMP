#ifndef included_AMP_PelletStackMechanicsSolverParameters
#define included_AMP_PelletStackMechanicsSolverParameters

#include "solvers/ColumnSolver.h"

namespace AMP {
namespace Solver {

class PelletStackMechanicsSolverParameters : public SolverStrategyParameters
{
public:
    PelletStackMechanicsSolverParameters() {}
    explicit PelletStackMechanicsSolverParameters( const AMP::shared_ptr<AMP::Database> db )
        : SolverStrategyParameters( db )
    {
    }
    virtual ~PelletStackMechanicsSolverParameters() {}

    AMP::shared_ptr<AMP::Solver::ColumnSolver> d_columnSolver;

protected:
private:
};
}
}

#endif
