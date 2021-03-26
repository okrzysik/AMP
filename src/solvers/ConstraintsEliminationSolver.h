#ifndef included_AMP_ConstraintsEliminationSolver
#define included_AMP_ConstraintsEliminationSolver

#include "AMP/solvers/SolverStrategy.h"

namespace AMP {
namespace Solver {

typedef SolverStrategyParameters ConstraintsEliminationSolverParameters;

class ConstraintsEliminationSolver : public SolverStrategy
{
public:
    explicit ConstraintsEliminationSolver(
        std::shared_ptr<ConstraintsEliminationSolverParameters> params );
    virtual void apply( std::shared_ptr<const AMP::LinearAlgebra::Vector> f,
                        std::shared_ptr<AMP::LinearAlgebra::Vector> u ) override;
};
} // namespace Solver
} // namespace AMP

#endif
