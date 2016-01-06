#ifndef included_SolverStrategyParameters
#define included_SolverStrategyParameters

#include "operators/Operator.h"
#include "utils/Database.h"
#include "utils/ParameterBase.h"
#include "utils/shared_ptr.h"


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
    explicit SolverStrategyParameters( const AMP::shared_ptr<AMP::Database> &db );

    /**
     * Destructor.
     */
    virtual ~SolverStrategyParameters();

    /**
    *  Pointer to database object which needs to be initialized specific to the solver.
    *  Documentation for parameters required by each solver can be found in the
    *  documentation for the solver.
    */
    AMP::shared_ptr<AMP::Database> d_db;

    AMP::shared_ptr<AMP::Operator::Operator> d_pOperator;

protected:
private:
};
}
}

#endif //  included_SolverStrategyParameters
