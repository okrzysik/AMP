
#ifndef included_AMP_TrilinosMueLuSolver
#define included_AMP_TrilinosMueLuSolver

#include "matrices/trilinos/EpetraMatrix.h"
#include "solvers/SolverStrategy.h"
#include "solvers/SolverStrategyParameters.h"
#include "solvers/trilinos/ml/MLoptions.h"

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-parameter"
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wextra"
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-local-typedefs"
#include "Teuchos_ParameterList.hpp"
//#include <MueLu.hpp>
#pragma GCC diagnostic pop
#pragma GCC diagnostic pop
#pragma GCC diagnostic pop

namespace MueLu{
  class EpetraOperator;
}


namespace AMP {
namespace Solver {


typedef SolverStrategyParameters TrilinosMueLuSolverParameters;


/**
 * The TrilinosMueLuSolver is a wrapper to the Trilinos ML solver. ML provides implementations of
 * various algebraic multigrid methods. The wrapper at present simply provides an adaptor
 * to enable AMP users to use the black box ML preconditioner.
 */

class TrilinosMueLuSolver : public SolverStrategy
{

public:
    /**
     * Default constructor
     */
    TrilinosMueLuSolver();

    /**
     * Main constructor.
     @param [in] parameters The parameters object contains a database object which must contain the
     following fields in addition to the fields expected by the base class SolverStrategy class:

    */
    explicit TrilinosMueLuSolver( AMP::shared_ptr<TrilinosMueLuSolverParameters> parameters );

    /**
     * Default destructor
     */
    virtual ~TrilinosMueLuSolver();

    /**
     * Solve the system \f$Au = f\f$.
     @param [in] f : shared pointer to right hand side vector
     @param [out] u : shared pointer to approximate computed solution
     */
    void solve( AMP::shared_ptr<const AMP::LinearAlgebra::Vector> f,
                AMP::shared_ptr<AMP::LinearAlgebra::Vector> u );

    /**
     * Return a shared pointer to the ML_Epetra::MultiLevelPreconditioner object
     */
    inline const AMP::shared_ptr<MueLu::EpetraOperator> getMLSolver( void )
    {
        return d_mueluSolver;
    }

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
     * Currently every call to reset destroys the ML preconditioner object
     * and recreates it based on the parameters object. See constructor for
     * fields required for parameter object.
     */
    void reset( AMP::shared_ptr<SolverStrategyParameters> params );

protected:
    void reSolveWithLU( AMP::shared_ptr<const AMP::LinearAlgebra::Vector> f,
                        AMP::shared_ptr<AMP::LinearAlgebra::Vector>
                            u );

    void getFromInput( const AMP::shared_ptr<AMP::Database> &db );

private:
    bool d_bUseEpetra;

    AMP_MPI d_comm;

    bool d_bCreationPhase; /**< set to true if the PC is not ready and false otherwise. */
    bool d_bRobustMode;

    AMP::shared_ptr<MueLu::EpetraOperator> d_mueluSolver;

    AMP::shared_ptr<AMP::LinearAlgebra::EpetraMatrix> d_matrix;
    Teuchos::ParameterList d_MueLuParameterList;

};
}
}

#endif
