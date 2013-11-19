#ifndef included_AMP_Flow1DSolver
#define included_AMP_Flow1DSolver

#include "SolverStrategy.h"
#include "SolverStrategyParameters.h"

namespace AMP {
namespace Solver {

  typedef SolverStrategyParameters Flow1DSolverParameters;

  class Flow1DSolver: public SolverStrategy {
    public:

      Flow1DSolver(boost::shared_ptr<Flow1DSolverParameters> parameters);
      
     virtual ~Flow1DSolver();

    void solve(boost::shared_ptr<const AMP::LinearAlgebra::Vector>  f,
	       boost::shared_ptr<AMP::LinearAlgebra::Vector>  u);
    
    void setInitialGuess( boost::shared_ptr<AMP::LinearAlgebra::Vector>  initialGuess );

    void initialize(boost::shared_ptr<SolverStrategyParameters> const parameters);

    void reset(boost::shared_ptr<SolverStrategyParameters> );
    
    void resetOperator(const boost::shared_ptr<AMP::Operator::OperatorParameters> params);

    AMP::LinearAlgebra::Variable::shared_ptr getInputVariable(int varId = -1) {
      return d_inpVariable;
    }

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

    boost::shared_ptr<AMP::LinearAlgebra::Variable> d_inpVariable; 
    boost::shared_ptr<AMP::LinearAlgebra::Variable> d_outVariable; 

    private:

    AMP::LinearAlgebra::Vector::shared_ptr d_Rhs;
    AMP::LinearAlgebra::Vector::shared_ptr d_Sol;
    
  };

}
}

#endif


