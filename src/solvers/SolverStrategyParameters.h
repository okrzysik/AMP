#ifndef included_SolverStrategyParameters
#define included_SolverStrategyParameters

#ifndef included_AMP_config

#endif

#ifndef included_Pointer
#include "utils/shared_ptr.h"
#endif

#ifndef included_Database
#include "utils/Database.h"
#endif

#ifndef included_ParameterBase
#include "utils/ParameterBase.h"
#endif

#ifndef included_AMP_Operator
#include "operators/Operator.h"
#endif


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
