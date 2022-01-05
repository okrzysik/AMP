#ifndef included_AMP_Flow1DSolver
#define included_AMP_Flow1DSolver

#include "AMP/solvers/SolverStrategy.h"
#include "AMP/solvers/SolverStrategyParameters.h"

namespace AMP::Solver {

typedef SolverStrategyParameters Flow1DSolverParameters;


class Flow1DSolver : public SolverStrategy
{
public:
    explicit Flow1DSolver( std::shared_ptr<Flow1DSolverParameters> parameters );

    virtual ~Flow1DSolver();

    void apply( std::shared_ptr<const AMP::LinearAlgebra::Vector> f,
                std::shared_ptr<AMP::LinearAlgebra::Vector> u ) override;

    void setInitialGuess( std::shared_ptr<AMP::LinearAlgebra::Vector> initialGuess ) override;

    void initialize( std::shared_ptr<const SolverStrategyParameters> parameters ) override;

    void reset( std::shared_ptr<SolverStrategyParameters> ) override;

    void resetOperator( std::shared_ptr<const AMP::Operator::OperatorParameters> params ) override;

    std::shared_ptr<AMP::LinearAlgebra::Variable> getInputVariable( int varId = -1 );

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

    std::shared_ptr<AMP::LinearAlgebra::Variable> d_inpVariable;
    std::shared_ptr<AMP::LinearAlgebra::Variable> d_outVariable;

private:
    AMP::LinearAlgebra::Vector::shared_ptr d_Rhs;
    AMP::LinearAlgebra::Vector::shared_ptr d_Sol;
};
} // namespace AMP::Solver

#endif
