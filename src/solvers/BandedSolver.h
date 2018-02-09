#ifndef included_BandedSolver
#define included_BandedSolver

#include "AMP/operators/Operator.h"
#include "AMP/solvers/SolverStrategy.h"
#include "AMP/vectors/Vector.h"


namespace AMP {
namespace Solver {

/**
 * Class BandedSolver is a derived class to solve
 * equations of the form \f$A(u) = f\f$ where $A$ is a
 * banded matrix
 */

class BandedSolver : public SolverStrategy
{
public:
    /**
     *  Main constructor for the base class.
     *  @param[in] parameters   The parameters object contains a database object which must contain
     * the
     *                          following fields:
     *                          1. type: integer, name: KL (required)
     *                             acceptable values (non-negative integer values)
     *                          1. type: integer, name: KU (required)
     *                             acceptable values (non-negative integer values)
     */
    explicit BandedSolver( AMP::shared_ptr<SolverStrategyParameters> parameters );

    /**
     * Default destructor. Currently does not do anything.
     */
    virtual ~BandedSolver();

    /**
     * Solve the system \f$A(u) = f\f$.  This is a pure virtual function that the derived classes
     * need to provide an implementation of.
     * @param[in]  f    shared pointer to right hand side vector
     * @param[out] u    shared pointer to approximate computed solution
     */
    virtual void solve( AMP::shared_ptr<const AMP::LinearAlgebra::Vector> f,
                        AMP::shared_ptr<AMP::LinearAlgebra::Vector> u );

    /**
     * Resets the operator registered with the solver with new parameters if necessary
     * @param parameters
     *        OperatorParameters object that is NULL by default
     */
    virtual void
    resetOperator( const AMP::shared_ptr<AMP::Operator::OperatorParameters> parameters );

    /**
     * Resets the solver internally with new parameters if necessary
     * @param parameters
     *        BandedSolverParameters object that is NULL by default
     */
    virtual void reset( AMP::shared_ptr<SolverStrategyParameters> parameters );


protected:
    BandedSolver();


private:
    int N, M, KL, KU;
    double *AB;
    int *IPIV;
    AMP::Discretization::DOFManager::shared_ptr rightDOF;
    AMP::Discretization::DOFManager::shared_ptr leftDOF;
};
} // namespace Solver
} // namespace AMP

#endif
