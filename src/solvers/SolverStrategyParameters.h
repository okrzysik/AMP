#ifndef included_AMP_SolverStrategyParameters
#define included_AMP_SolverStrategyParameters

#include "AMP/operators/Operator.h"
#include "AMP/utils/Database.h"
#include "AMP/utils/ParameterBase.h"
#include <memory>


namespace AMP::Solver {

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

    /**
     * List of vectors to be used during solver initialization
     */
    std::vector<std::shared_ptr<AMP::LinearAlgebra::Vector>> d_vectors;

    /** Pointer to global database
     *  This is temporary fix and eventually either d_global_db or d_db should go away
     *  This is introduced to allow for solver factories to access databases in the global
     *  database for the construction of nested solvers
     */
    std::shared_ptr<AMP::Database> d_global_db;

protected:
private:
};
} // namespace AMP::Solver

#endif
