#ifndef included_AMP_PelletStackMechanicsSolver
#define included_AMP_PelletStackMechanicsSolver

#include "operators/libmesh/PelletStackOperator.h"
#include "solvers/ColumnSolver.h"
#include "solvers/SolverStrategy.h"
#include "solvers/libmesh/PelletStackMechanicsSolverParameters.h"

namespace AMP {
namespace Solver {


class PelletStackMechanicsSolver : public SolverStrategy
{
public:
    explicit PelletStackMechanicsSolver(
        AMP::shared_ptr<PelletStackMechanicsSolverParameters> params );

    virtual ~PelletStackMechanicsSolver() {}

    void resetOperator( const AMP::shared_ptr<AMP::Operator::OperatorParameters> params );

    void solve( AMP::shared_ptr<const AMP::LinearAlgebra::Vector> f,
                AMP::shared_ptr<AMP::LinearAlgebra::Vector>
                    u );

protected:
    void solveSerial( AMP::shared_ptr<const AMP::LinearAlgebra::Vector> f,
                      AMP::shared_ptr<AMP::LinearAlgebra::Vector>
                          u );

    void solveScan( AMP::shared_ptr<const AMP::LinearAlgebra::Vector> f,
                    AMP::shared_ptr<AMP::LinearAlgebra::Vector>
                        u );

    AMP::shared_ptr<AMP::Operator::PelletStackOperator> d_pelletStackOp;
    AMP::shared_ptr<AMP::Solver::ColumnSolver> d_columnSolver;
    AMP::shared_ptr<AMP::LinearAlgebra::Vector> d_fbuffer1;
    AMP::shared_ptr<AMP::LinearAlgebra::Vector> d_fbuffer2;
};
}
}

#endif
