#ifndef included_AMP_SolverStrategyParameters
#define included_AMP_SolverStrategyParameters

#include "AMP/operators/Operator.h"
#include "AMP/utils/Database.h"
#include "AMP/utils/ParameterBase.h"
#include <memory>


namespace AMP {
namespace Solver {

/**\class SolverStrategyParameters
 *
 * SolverStrategyParameters encapsulates parameters used to initialize
 * SolverStrategy objects
 */

class SolverStrategyParameters : public ParameterBase
{
public:
    /**
     * Empty constructor.
     */
    SolverStrategyParameters();

    /**
     * Construct and initialize a parameter list according to input
     * data.  Guess what the required and optional keywords are.
     */
    explicit SolverStrategyParameters( std::shared_ptr<AMP::Database> db );

    /**
     * Destructor.
     */
    virtual ~SolverStrategyParameters();

    /**
     *  Pointer to database object which needs to be initialized specific to the solver.
     *  Documentation for parameters required by each solver can be found in the
     *  documentation for the solver.
     */
    std::shared_ptr<AMP::Database> d_db = nullptr;

    std::shared_ptr<AMP::Operator::Operator> d_pOperator = nullptr;

protected:
private:
};
} // namespace Solver
} // namespace AMP

#endif
