#ifndef included_AMP_BandedSolver
#define included_AMP_BandedSolver

#include "AMP/operators/Operator.h"
#include "AMP/solvers/SolverStrategy.h"
#include "AMP/vectors/Vector.h"


namespace AMP::Solver {

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
     *  @param[in] params   The parameters object contains a database object which must contain
     * the
     *                          following fields:
     *                          1. type: integer, name: KL (required)
     *                             acceptable values (non-negative integer values)
     *                          1. type: integer, name: KU (required)
     *                             acceptable values (non-negative integer values)
     */
    explicit BandedSolver( std::shared_ptr<SolverStrategyParameters> params );

    /**
     * static create routine that is used by SolverFactory
     @param [in] params The parameters object
     contains a database objects with the fields listed for the constructor above
     */
    static std::unique_ptr<SolverStrategy>
    createSolver( std::shared_ptr<SolverStrategyParameters> params )
    {
        return std::make_unique<BandedSolver>( params );
    }

    /**
     * Default destructor. Currently does not do anything.
     */
    virtual ~BandedSolver();

    std::string type() const override { return "BandedSolver"; }

    /**
     * Solve the system \f$A(u) = f\f$.  This is a pure virtual function that the derived classes
     * need to provide an implementation of.
     * @param[in]  f    shared pointer to right hand side vector
     * @param[out] u    shared pointer to approximate computed solution
     */
    void apply( std::shared_ptr<const AMP::LinearAlgebra::Vector> f,
                std::shared_ptr<AMP::LinearAlgebra::Vector> u ) override;

    /**
     * Resets the operator registered with the solver with new parameters if necessary
     * @param params  OperatorParameters object that is NULL by default
     */
    void resetOperator( std::shared_ptr<const AMP::Operator::OperatorParameters> params ) override;

    /**
     * Resets the solver internally with new parameters if necessary
     * @param params  BandedSolverParameters object that is NULL by default
     */
    void reset( std::shared_ptr<SolverStrategyParameters> params ) override;


protected:
    BandedSolver();


private:
    int N, M, KL, KU;
    double *AB;
    int *IPIV;
    std::shared_ptr<AMP::Discretization::DOFManager> rightDOF;
    std::shared_ptr<AMP::Discretization::DOFManager> leftDOF;
};
} // namespace AMP::Solver

#endif
