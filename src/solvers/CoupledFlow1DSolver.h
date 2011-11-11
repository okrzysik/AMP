#ifndef included_AMP_CoupledFlow1DSolver
#define included_AMP_CoupledFlow1DSolver

#include "SolverStrategy.h"
#include "Flow1DSolver.h"
#include "SolverStrategyParameters.h"
#include "CoupledFlow1DSolverParameters.h"

namespace AMP {
namespace Solver {

  class CoupledFlow1DSolver: public SolverStrategy {
    public:

      CoupledFlow1DSolver(boost::shared_ptr<SolverStrategyParameters> parameters);
      
      ~CoupledFlow1DSolver();

    void solve(boost::shared_ptr<AMP::LinearAlgebra::Vector>  f,
	       boost::shared_ptr<AMP::LinearAlgebra::Vector>  u);
    
    void setInitialGuess( boost::shared_ptr<AMP::LinearAlgebra::Vector>  initialGuess );

    void reset(boost::shared_ptr<SolverStrategyParameters> );

    void resetOperator(const boost::shared_ptr<AMP::Operator::OperatorParameters> params);

    protected:

    int d_numpoints; /**< Number of points in z direction */

    std::vector<double> zPoints; 
    
    boost::shared_ptr<AMP::LinearAlgebra::Variable> d_inpVariable; 
    boost::shared_ptr<AMP::LinearAlgebra::Variable> d_outVariable; 

    private:

    boost::shared_ptr<AMP::LinearAlgebra::Variable> d_SimpleVariable; /**< Simple Input Variable */

    AMP::LinearAlgebra::Vector::shared_ptr d_Rhs;
    AMP::LinearAlgebra::Vector::shared_ptr d_Sol;
    AMP::LinearAlgebra::Vector::shared_ptr d_flowInput;
    AMP::LinearAlgebra::Vector::shared_ptr d_flowOutput;
    
    //AMP::Operator::MapOperator::shared_ptr d_flowInternal3to1;
    //AMP::Operator::MapOperator::shared_ptr d_flowInternal1to3;

    boost::shared_ptr<AMP::Solver::Flow1DSolver> d_flow1DSolver;
  };

}
}

#endif


