#ifndef included_AMP_CoupledFlow1DSolver
#define included_AMP_CoupledFlow1DSolver

#include "AMP/operators/map/MapOperator.h"
#include "AMP/solvers/SolverStrategy.h"
#include "AMP/solvers/SolverStrategyParameters.h"
#include "AMP/solvers/libmesh/CoupledFlow1DSolverParameters.h"
#include "AMP/solvers/libmesh/Flow1DSolver.h"


namespace AMP {
namespace Solver {

class CoupledFlow1DSolver : public SolverStrategy
{
public:
    explicit CoupledFlow1DSolver( std::shared_ptr<SolverStrategyParameters> parameters );

    virtual ~CoupledFlow1DSolver();

    void solve( std::shared_ptr<const AMP::LinearAlgebra::Vector> f,
                std::shared_ptr<AMP::LinearAlgebra::Vector> u ) override;

    void setInitialGuess( std::shared_ptr<AMP::LinearAlgebra::Vector> initialGuess ) override;

    void reset( std::shared_ptr<SolverStrategyParameters> ) override;

    void resetOperator( const std::shared_ptr<AMP::Operator::OperatorParameters> params ) override;

protected:
    int d_numpoints; /**< Number of points in z direction */

    std::vector<double> zPoints;

    std::shared_ptr<AMP::LinearAlgebra::Variable> d_inpVariable;
    std::shared_ptr<AMP::LinearAlgebra::Variable> d_outVariable;

private:
    std::shared_ptr<AMP::LinearAlgebra::Variable> d_SimpleVariable; /**< Simple Input Variable */

    AMP::LinearAlgebra::Vector::const_shared_ptr d_Rhs;
    AMP::LinearAlgebra::Vector::shared_ptr d_Sol;
    AMP::LinearAlgebra::Vector::shared_ptr d_flowInput;
    AMP::LinearAlgebra::Vector::shared_ptr d_flowOutput;

    AMP::Operator::MapOperator::shared_ptr d_flowInternal3to1;
    AMP::Operator::MapOperator::shared_ptr d_flowInternal1to3;

    std::shared_ptr<AMP::Solver::Flow1DSolver> d_flow1DSolver;
};
} // namespace Solver
} // namespace AMP

#endif
