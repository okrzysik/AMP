#ifndef included_AMP_TrilinosNOXSolver
#define included_AMP_TrilinosNOXSolver

// AMP includes
#include "AMP/solvers/SolverStrategy.h"
#include "AMP/solvers/trilinos/nox/TrilinosNOXSolverParameters.h"
#include "AMP/solvers/trilinos/thyra/TrilinosThyraModelEvaluator.h"


// Trilinos includes
DISABLE_WARNINGS
#include "NOX_Solver_Factory.H"
#include "NOX_StatusTest_Combo.H"
#include "NOX_Thyra.H"
#include <NOX_Solver_Generic.H>
ENABLE_WARNINGS


namespace AMP::Solver {


/**
 * The TrilinosNOXSolver is a wrapper to the PETSc SNES solver which provides an implementation
 * of the inexact Newton method.
 */
class TrilinosNOXSolver : public SolverStrategy
{
public:
    /**
     * default constructor, sets default values for member variables
     */
    TrilinosNOXSolver();

    /**
     * main constructor
     @param [in] parameters The parameters object
     contains a database objects containing the following fields:

     1. type: string, name: SNESOptions, default value: "",

     2. type: bool, name: usesJacobian, default value: false,
     acceptable values (true, false)

     3. name: MFFDDifferencingStrategy, type:string , default value: MATMFFD_WP,
     acceptable values ()

     4. name: MFFDFunctionDifferencingError, type: double, default value: PETSC_DEFAULT,
     acceptable values ()

     5. name: maximumFunctionEvals, type: integer, default value: none
     acceptable values ()

     6. name: absoluteTolerance, type: double, default value: none
     acceptable values ()

     7. name: relativeTolerance, type: double, default value: none
     acceptable values ()

     8. name: stepTolerance, type: double, default value: none
     acceptable values ()

     9. name: enableLineSearchPreCheck, type: bool, default value: FALSE
     acceptable values ()

     10. name: numberOfLineSearchPreCheckAttempts, type: integer, default value: 5
     acceptable values (non negative integer values)

     11. name: enableMFFDBoundsCheck, type: bool, default value: FALSE
     acceptable values ()

     12. name: operatorComponentToEnableBoundsCheck, type: integer, default value: none
     acceptable values ()
    */
    explicit TrilinosNOXSolver( std::shared_ptr<TrilinosNOXSolverParameters> parameters );

    /**
     * Default destructor.
     */
    virtual ~TrilinosNOXSolver();

    /**
     * Resets the solver internally with new parameters if necessary
     * @param parameters
     *        SolverStrategyParameters object that is NULL by default
     */
    void reset( std::shared_ptr<SolverStrategyParameters> parameters ) override;

    /**
     * Solve the system \f$Au = 0\f$.
     * @param [in] f : shared pointer to right hand side vector
     * @param [out] u : shared pointer to approximate computed solution
     */
    void apply( std::shared_ptr<const AMP::LinearAlgebra::Vector> f,
                std::shared_ptr<AMP::LinearAlgebra::Vector> u ) override;

    /**
     * Provide the initial guess for the solver. This is a pure virtual function that the derived
     * classes
     * need to provide an implementation of.
     * @param [in] initialGuess: shared pointer to the initial guess vector.
     */
    void setInitialGuess( std::shared_ptr<AMP::LinearAlgebra::Vector> initialGuess ) override
    {
        d_initialGuess = initialGuess;
    }


protected:
    void initialize( std::shared_ptr<const SolverStrategyParameters> parameters ) override;

    AMP_MPI d_comm;

    AMP::LinearAlgebra::Vector::shared_ptr d_initialGuess;

    Teuchos::RCP<NOX::Solver::Generic> d_solver;

    Teuchos::RCP<TrilinosThyraModelEvaluator> d_thyraModel;
    Teuchos::RCP<Teuchos::ParameterList> d_nlParams;
    Teuchos::RCP<Thyra::PreconditionerBase<double>> d_precOp;
    Teuchos::RCP<::Thyra::LinearOpWithSolveFactoryBase<double>> d_lowsFactory;
    Teuchos::RCP<NOX::StatusTest::Combo> d_status;
};
} // namespace AMP::Solver

#endif
