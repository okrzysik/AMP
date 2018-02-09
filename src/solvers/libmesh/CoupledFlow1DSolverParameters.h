#ifndef included_AMP_CoupledFlowFrapconParameters
#define included_AMP_CoupledFlowFrapconParameters

#include "AMP/operators/subchannel/CoupledFlowFrapconOperator.h"
#include "AMP/solvers/SolverStrategy.h"
#include "AMP/solvers/SolverStrategyParameters.h"
#include "AMP/utils/Database.h"
#include "AMP/utils/shared_ptr.h"

namespace AMP {
namespace Solver {

class CoupledFlow1DSolverParameters : public SolverStrategyParameters
{
public:
    CoupledFlow1DSolverParameters() {}
    explicit CoupledFlow1DSolverParameters( const AMP::shared_ptr<AMP::Database> &db )
        : SolverStrategyParameters( db )
    {
    }
    virtual ~CoupledFlow1DSolverParameters() {}

    AMP::shared_ptr<AMP::Solver::SolverStrategy> d_flow1DSolver;
    using SolverStrategyParameters::d_pOperator;

protected:
private:
};
} // namespace Solver
} // namespace AMP

#endif
