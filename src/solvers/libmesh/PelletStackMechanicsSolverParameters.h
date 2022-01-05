#ifndef included_AMP_PelletStackMechanicsSolverParameters
#define included_AMP_PelletStackMechanicsSolverParameters

#include "AMP/solvers/ColumnSolver.h"

namespace AMP::Solver {

class PelletStackMechanicsSolverParameters : public SolverStrategyParameters
{
public:
    PelletStackMechanicsSolverParameters() {}
    explicit PelletStackMechanicsSolverParameters( std::shared_ptr<AMP::Database> db )
        : SolverStrategyParameters( db )
    {
    }
    virtual ~PelletStackMechanicsSolverParameters() {}

    std::shared_ptr<AMP::Solver::ColumnSolver> d_columnSolver;

protected:
private:
};
} // namespace AMP::Solver

#endif
