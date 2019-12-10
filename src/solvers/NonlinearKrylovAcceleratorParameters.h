#ifndef included_AMP_NonlinearKrylovAcceleratorParameters
#define included_AMP_NonlinearKrylovAcceleratorParameters

#include <string>

#include "AMP/utils/Database.h"

#include "AMP/utils/ParameterBase.h"

#ifndef included_AMP_SolverStrategy
#include "SolverStrategy.h"
#endif

#include "AMP/vectors/Vector.h"

namespace AMP {
namespace Solver {

/** \class NonlinearKrylovAcceleratorParameters
 *
 * Class NonlinearKrylovAcceleratorParameters provides a uniform mechanism to pass
 * initialization parameters to the NonlinearKrylovAccelerator solver. It contains
 * shared pointers to a database object, a preconditioner, and a solution vector for
 * storing initial guesses. All member variables are public.
 */
class NonlinearKrylovAcceleratorParameters : public SolverStrategyParameters
{
public:
    /**
     * Empty constructor.
     */
    NonlinearKrylovAcceleratorParameters();

    /**
     * Construct and initialize a parameter list according to input
     * data.  See Application for a list of required and optional keywords.
     */
    explicit NonlinearKrylovAcceleratorParameters( AMP::shared_ptr<AMP::Database> database );

    /**
     * Destructor.
     */
    virtual ~NonlinearKrylovAcceleratorParameters();

    AMP::shared_ptr<SolverStrategy> d_pPreconditioner;

    AMP::shared_ptr<AMP::LinearAlgebra::Vector> d_pInitialGuess;
};
} // namespace Solver
} // namespace AMP

#endif
