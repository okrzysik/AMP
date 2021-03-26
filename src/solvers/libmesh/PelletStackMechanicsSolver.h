#ifndef included_AMP_PelletStackMechanicsSolver
#define included_AMP_PelletStackMechanicsSolver

#include "AMP/operators/libmesh/PelletStackOperator.h"
#include "AMP/solvers/ColumnSolver.h"
#include "AMP/solvers/SolverStrategy.h"
#include "AMP/solvers/libmesh/PelletStackMechanicsSolverParameters.h"

namespace AMP {
namespace Solver {


class PelletStackMechanicsSolver : public SolverStrategy
{
public:
    explicit PelletStackMechanicsSolver(
        std::shared_ptr<PelletStackMechanicsSolverParameters> params );

    virtual ~PelletStackMechanicsSolver() {}

    virtual void
    resetOperator( const std::shared_ptr<AMP::Operator::OperatorParameters> params ) override;

    virtual void apply( std::shared_ptr<const AMP::LinearAlgebra::Vector> f,
                        std::shared_ptr<AMP::LinearAlgebra::Vector> u ) override;

protected:
    void solveSerial( std::shared_ptr<const AMP::LinearAlgebra::Vector> f,
                      std::shared_ptr<AMP::LinearAlgebra::Vector> u );

    void solveScan( std::shared_ptr<const AMP::LinearAlgebra::Vector> f,
                    std::shared_ptr<AMP::LinearAlgebra::Vector> u );

    std::shared_ptr<AMP::Operator::PelletStackOperator> d_pelletStackOp;
    std::shared_ptr<AMP::Solver::ColumnSolver> d_columnSolver;
    std::shared_ptr<AMP::LinearAlgebra::Vector> d_fbuffer1;
    std::shared_ptr<AMP::LinearAlgebra::Vector> d_fbuffer2;
};
} // namespace Solver
} // namespace AMP

#endif
