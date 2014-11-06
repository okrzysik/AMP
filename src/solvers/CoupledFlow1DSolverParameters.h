#ifndef included_AMP_CoupledFlowFrapconParameters
#define included_AMP_CoupledFlowFrapconParameters

#include "utils/shared_ptr.h"
#include "utils/Database.h"
#include "operators/subchannel/CoupledFlowFrapconOperator.h"
#include "SolverStrategyParameters.h"
#include "SolverStrategy.h"

namespace AMP {
namespace Solver {

  class CoupledFlow1DSolverParameters: public SolverStrategyParameters{
  public:
    CoupledFlow1DSolverParameters(){}
    CoupledFlow1DSolverParameters(const AMP::shared_ptr<AMP::Database> &db): 
      SolverStrategyParameters(db){ }
    virtual ~CoupledFlow1DSolverParameters(){}

    AMP::shared_ptr<AMP::Solver::SolverStrategy> d_flow1DSolver;
    AMP::shared_ptr<AMP::Operator::Operator> d_pOperator;

  protected:
  private:
    
  };

}
}

#endif
