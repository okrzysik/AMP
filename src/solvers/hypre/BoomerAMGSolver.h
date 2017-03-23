#ifndef included_AMP_BoomerAMGSolver
#define included_AMP_BoomerAMGSolver

#include "solvers/SolverStrategy.h"
#include "solvers/SolverStrategyParameters.h"

#include "HYPRE.h"
#include "HYPRE_IJ_mv.h"

namespace AMP {
namespace Solver {


typedef SolverStrategyParameters BoomerAMGSolverParameters;


/**
 * The BoomerAMGSolver is a wrapper to the HYPRE BoomerAMG solver. BoomerAMG provides implementations of
 * various algebraic multigrid methods. The wrapper at present simply provides an adaptor
 * to enable AMP users to use the black box BoomerAMG preconditioner.
 */

class BoomerAMGSolver : public SolverStrategy
{

public:
    /**
     * Default constructor
     */
    BoomerAMGSolver();

    /**
     * Main constructor.
     @param [in] parameters The parameters object contains a database object which must contain the
     following fields in addition to the fields expected by the base class SolverStrategy class:

    */
    explicit BoomerAMGSolver( AMP::shared_ptr<BoomerAMGSolverParameters> parameters );

    /**
     * Default destructor
     */
    virtual ~BoomerAMGSolver();

    /**
     * Solve the system \f$Au = f\f$.
     @param [in] f : shared pointer to right hand side vector
     @param [out] u : shared pointer to approximate computed solution
     */
    void solve( AMP::shared_ptr<const AMP::LinearAlgebra::Vector> f,
                AMP::shared_ptr<AMP::LinearAlgebra::Vector> u );

    /**
     * Initialize the solution vector and potentially create internal vectors needed for solution
     @param [in] parameters The parameters object
     contains a database object. Refer to the documentation for the constructor to see what fields
     are required.
     This routine assumes that a non-NULL operator of type LinearOperator has been registered with
     the solver.
     The LinearOperator currently is assumed to contain a pointer to an EpetraMatrix object.
     */
    void initialize( AMP::shared_ptr<SolverStrategyParameters> const parameters );

    /**
     * Register the operator that the solver will use during solves
     @param [in] op shared pointer to the linear operator $A$ for equation \f$A u = f\f$
     */
    void registerOperator( const AMP::shared_ptr<AMP::Operator::Operator> op );

    /**
     * Resets the associated operator internally with new parameters if necessary
     * @param [in] params
     *        OperatorParameters object that is NULL by default
     */
    void resetOperator( const AMP::shared_ptr<AMP::Operator::OperatorParameters> params );

    /**
     * Resets the solver internally with new parameters if necessary
     * @param [in] params
     *        SolverStrategyParameters object that is NULL by default
     * Currently every call to reset destroys the BoomerAMG preconditioner object
     * and recreates it based on the parameters object. See constructor for
     * fields required for parameter object.
     */
    void reset( AMP::shared_ptr<SolverStrategyParameters> params );

    void getFromInput( const AMP::shared_ptr<AMP::Database> &db );

private:

    /**
     * create the internal HYPRE_IJMatrix based on the AMP matrix
     */
    void createHYPREMatrix( const AMP::shared_ptr<AMP::LinearAlgebra::Matrix> matrix );

    AMP_MPI d_comm;

    HYPRE_IJMatrix d_ijMatrix;  //! pointer to HYPRE matrix struct

    bool d_bCreationPhase; /**< set to true if the PC is not ready and false otherwise. */
};
}
}

#endif
