#ifndef included_AMP_Flow1DSolver
#define included_AMP_Flow1DSolver

#include "AMP/solvers/SolverStrategy.h"
#include "AMP/solvers/SolverStrategyParameters.h"

namespace AMP {
namespace Solver {

typedef SolverStrategyParameters Flow1DSolverParameters;


class Flow1DSolver : public SolverStrategy
{
public:
    explicit Flow1DSolver( AMP::shared_ptr<Flow1DSolverParameters> parameters );

    virtual ~Flow1DSolver();

    virtual void solve( AMP::shared_ptr<const AMP::LinearAlgebra::Vector> f,
                        AMP::shared_ptr<AMP::LinearAlgebra::Vector> u ) override;

    virtual void
    setInitialGuess( AMP::shared_ptr<AMP::LinearAlgebra::Vector> initialGuess ) override;

    virtual void initialize( AMP::shared_ptr<SolverStrategyParameters> const parameters ) override;

    virtual void reset( AMP::shared_ptr<SolverStrategyParameters> ) override;

    virtual void
    resetOperator( const AMP::shared_ptr<AMP::Operator::OperatorParameters> params ) override;

    AMP::LinearAlgebra::Variable::shared_ptr getInputVariable( int varId = -1 );

protected:
    int d_numpoints; /**< Number of points in z direction */

    std::vector<double> zPoints;

    double d_Cp;

    double d_dCp;

    AMP::LinearAlgebra::Vector::shared_ptr d_cladVec;

    double d_De;

    double d_G;

    double d_K;

    double d_Re;

    double d_Pr;

    AMP::shared_ptr<AMP::LinearAlgebra::Variable> d_inpVariable;
    AMP::shared_ptr<AMP::LinearAlgebra::Variable> d_outVariable;

private:
    AMP::LinearAlgebra::Vector::shared_ptr d_Rhs;
    AMP::LinearAlgebra::Vector::shared_ptr d_Sol;
};
} // namespace Solver
} // namespace AMP

#endif
