#ifndef included_AMP_CoupledFlow1DSolver
#define included_AMP_CoupledFlow1DSolver

#include "operators/map/MapOperator.h"
#include "solvers/SolverStrategy.h"
#include "solvers/SolverStrategyParameters.h"
#include "solvers/libmesh/CoupledFlow1DSolverParameters.h"
#include "solvers/libmesh/Flow1DSolver.h"


namespace AMP {
namespace Solver {

class CoupledFlow1DSolver : public SolverStrategy
{
public:
    explicit CoupledFlow1DSolver( AMP::shared_ptr<SolverStrategyParameters> parameters );

    virtual ~CoupledFlow1DSolver();

    void solve( AMP::shared_ptr<const AMP::LinearAlgebra::Vector> f,
                AMP::shared_ptr<AMP::LinearAlgebra::Vector>
                    u );

    void setInitialGuess( AMP::shared_ptr<AMP::LinearAlgebra::Vector> initialGuess );

    void reset( AMP::shared_ptr<SolverStrategyParameters> );

    void resetOperator( const AMP::shared_ptr<AMP::Operator::OperatorParameters> params );

protected:
    int d_numpoints; /**< Number of points in z direction */

    std::vector<double> zPoints;

    AMP::shared_ptr<AMP::LinearAlgebra::Variable> d_inpVariable;
    AMP::shared_ptr<AMP::LinearAlgebra::Variable> d_outVariable;

private:
    AMP::shared_ptr<AMP::LinearAlgebra::Variable> d_SimpleVariable; /**< Simple Input Variable */

    AMP::LinearAlgebra::Vector::const_shared_ptr d_Rhs;
    AMP::LinearAlgebra::Vector::shared_ptr d_Sol;
    AMP::LinearAlgebra::Vector::shared_ptr d_flowInput;
    AMP::LinearAlgebra::Vector::shared_ptr d_flowOutput;

    AMP::Operator::MapOperator::shared_ptr d_flowInternal3to1;
    AMP::Operator::MapOperator::shared_ptr d_flowInternal1to3;

    AMP::shared_ptr<AMP::Solver::Flow1DSolver> d_flow1DSolver;
};
}
}

#endif
