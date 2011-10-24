#ifndef included_AMP_CoupledFlowFrapconParameters
#define included_AMP_CoupledFlowFrapconParameters

#ifndef included_Pointer
#include "boost/shared_ptr.hpp"
#endif

#ifndef included_Database
#include "utils/Database.h"
#endif

#include "operators/CoupledFlowFrapconOperator.h"

#ifndef included_AMP_SolverStrategyParameters
#include "SolverStrategyParameters.h"
#endif

#ifndef included_AMP_SolverStrategy
#include "SolverStrategy.h"
#endif

namespace AMP {
namespace Solver {

  class CoupledFlow1DSolverParameters: public SolverStrategyParameters{
  public:
    CoupledFlow1DSolverParameters(){}
    CoupledFlow1DSolverParameters(const boost::shared_ptr<AMP::Database> &db): 
      SolverStrategyParameters(db){ }
    ~CoupledFlow1DSolverParameters(){}

    boost::shared_ptr<AMP::Solver::SolverStrategy> d_flow1DSolver;
    boost::shared_ptr<AMP::Operator::Operator> d_pOperator;

  protected:
  private:
    
  };

}
}

#endif
