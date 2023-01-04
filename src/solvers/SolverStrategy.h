#ifndef included_AMP_SolverStrategy
#define included_AMP_SolverStrategy

#include "AMP/operators/Operator.h"
#include "AMP/vectors/Vector.h"
#include "SolverStrategyParameters.h"

#include <memory>


// Declare some classes
namespace AMP::IO {
class Writer;
} // namespace AMP::IO


namespace AMP::Solver {

/**
 * Class SolverStrategy is a base class for methods to solve
 * equations of the form \f$A(u) = f\f$. $A$ may be a nonlinear
 * or linear operator.
 */

class SolverStrategy
{
public:
    typedef std::shared_ptr<AMP::Solver::SolverStrategy> shared_ptr;

    /**
     * Default constructor
     */
    SolverStrategy();

    /**
     *  Main constructor for the base class.
     *  @param[in] parameters   The parameters object contains a database object which must contain
     * the
     *                          following fields:
     *                          1. type: integer, name: max_iterations (required)
     *                             acceptable values (non-negative integer values)
     *                          2. type: double, name: max_error, (required)
     *                             acceptable values (non-negative real values)
     *                          3. type: integer, name: print_info_level, default value: 0,
     *                             acceptable values (non negative integer values, the higher
     *                             the value the more verbose the debugging information provided)
     *                          4. type: bool, name: zero_initial_guess, default value: false,
     *                             acceptable values (TRUE, FALSE)
     */
    explicit SolverStrategy( std::shared_ptr<const SolverStrategyParameters> parameters );

    /**
     * Default destructor. Currently does not do anything.
     */
    virtual ~SolverStrategy();

    //! Return the name of the solver
    virtual std::string type() const = 0;

    enum class SolverStatus {
        ConvergedOnAbsTol,
        ConvergedOnRelTol,
        ConvergedIterations,
        ConvergedUserCondition,
        DivergedMaxIterations,
        DivergedLineSearch,
        DivergedStepSize,
        DivergedFunctionCount,
        DivergedOnNan,
        DivergedNestedSolver,
        DivergedOther
    };

    /**
     * Solve the system \f$A(u) = f\f$.  This is a pure virtual function that the derived classes
     * need to provide an implementation of.
     * @param[in]  f    shared pointer to right hand side vector
     * @param[out] u    shared pointer to approximate computed solution
     */
    virtual void apply( std::shared_ptr<const AMP::LinearAlgebra::Vector> f,
                        std::shared_ptr<AMP::LinearAlgebra::Vector> u ) = 0;

    /**
     * Initialize the solution vector and potentially create internal vectors needed for solution
     * @param[in] parameters    The parameters object contains a database object.
     *                          Currently there are no required fields for the database object.
     */
    virtual void initialize( std::shared_ptr<const SolverStrategyParameters> parameters );

    /**
     * Provide the initial guess for the solver. This is a pure virtual function that the derived
     * classes
     * need to provide an implementation of.
     * @param[in] initialGuess: shared pointer to the initial guess vector.
     */
    virtual void setInitialGuess( std::shared_ptr<AMP::LinearAlgebra::Vector> initialGuess );

    /**
     * Specify level of diagnostic information printed during iterations.
     * @param[in] print_level    integer level value with permissible values 0 and higher. Setting
     *                           to zero should provide minimial debugging information with higher
     *                           values resulting in increasingly verbose information being printed
     * out.
     */
    virtual void setDebugPrintInfoLevel( int print_level ) { d_iDebugPrintInfoLevel = print_level; }

    /**
     * Get level of diagnostic information printed during iterations.
     */
    int getDebugPrintInfoLevel( void ) { return d_iDebugPrintInfoLevel; }

    /**
     * Return the number of iterations taken by the solver to converge.
     */
    virtual int getIterations( void ) const { return ( d_iNumberIterations ); }

    /**
     * Tells the solver to use an initial guess of zero and not try to
     * copy an initial guess into the solution vector
     * @param[in] use_zero_guess    boolean to specify whether zero initial guess should be used or
     * not.
     */
    virtual void setZeroInitialGuess( bool use_zero_guess )
    {
        d_bUseZeroInitialGuess = use_zero_guess;
    }

    /**
     * Set a nested solver, eg, Krylov for Newton, preconditioner for Krylov etc. Null op in base
     * class
     */
    virtual void setNestedSolver( std::shared_ptr<SolverStrategy> ) {}

    /**
     * Return a nested solver (eg preconditioner) if it exists. By default return a nullptr
     */

    virtual std::shared_ptr<SolverStrategy> getNestedSolver( void ) { return nullptr; }

    /**
     * Register the operator that the solver will use during solves
     * @param [in] op shared pointer to operator \f$A()\f$ for equation \f$A(u) = f\f$
     */
    virtual void registerOperator( std::shared_ptr<AMP::Operator::Operator> op )
    {
        d_pOperator = op;
    }

    /**
     * \brief  Registers a writer with the solver
     * \details  This function will register a writer with the solver.  The solver
     *  may then register any vector components it "owns" with the writer.
     * \param writer   The writer to register
     */
    virtual void registerWriter( std::shared_ptr<AMP::IO::Writer> writer ) { d_writer = writer; }

    /**
     * Resets the operator registered with the solver with new parameters if necessary
     * @param parameters
     *        OperatorParameters object that is NULL by default
     */
    virtual void
    resetOperator( std::shared_ptr<const AMP::Operator::OperatorParameters> parameters );

    /**
     * Resets the solver internally with new parameters if necessary
     * @param parameters
     *        SolverStrategyParameters object that is NULL by default
     */
    virtual void reset( std::shared_ptr<SolverStrategyParameters> parameters );

    /**
     * Return a shared pointer to the operator registered with the solver.
     */
    virtual std::shared_ptr<AMP::Operator::Operator> getOperator( void ) { return d_pOperator; }

    /*!
     *  Get absolute tolerance for solver.
     */
    double getAbsoluteTolerance() const { return ( d_dAbsoluteTolerance ); }

    /*!
     *  Set absolute tolerance for nonlinear solver.
     */
    virtual void setAbsoluteTolerance( double abs_tol ) { d_dAbsoluteTolerance = abs_tol; }

    double getRelativeTolerance() const { return ( d_dRelativeTolerance ); }

    virtual void setRelativeTolerance( double rel_tol ) { d_dRelativeTolerance = rel_tol; }

    /**
     * Set the maximum number of iterations for the solver
     */
    virtual void setMaxIterations( const int max_iterations ) { d_iMaxIterations = max_iterations; }

    int getMaxIterations( void ) const { return d_iMaxIterations; }

    virtual void printStatistics( std::ostream &os = AMP::pout )
    {
        os << "Not implemented for this solver!" << std::endl;
    }

    int getConvergenceStatus( void ) const
    {
        return d_ConvergenceStatus <= SolverStatus::ConvergedUserCondition ? 1 : 0;
    }

    virtual void print( std::ostream &os ) { NULL_USE( os ); }

    /**
     * Return the residual norm.
     */
    virtual double getResidualNorm( void ) const { return d_dResidualNorm; }

    /**
     * returns whether the solver has converged or not
     */
    virtual bool checkConvergence( std::shared_ptr<const AMP::LinearAlgebra::Vector> residual );

    virtual const std::vector<int> &getIterationHistory( void ) { return d_iterationHistory; }

    int getTotalNumberOfIterations( void );

    virtual void residual( std::shared_ptr<const AMP::LinearAlgebra::Vector> f,
                           std::shared_ptr<const AMP::LinearAlgebra::Vector> u,
                           std::shared_ptr<AMP::LinearAlgebra::Vector> r );

    virtual void printConvergenceStatus( SolverStrategy::SolverStatus status,
                                         std::ostream &os = AMP::pout ) const
    {
    }

protected:
    void getFromInput( std::shared_ptr<AMP::Database> db );

    SolverStatus d_ConvergenceStatus = SolverStatus::DivergedOther;

    std::string d_sName;

    int d_iNumberIterations = 0; // iterations in solver

    int d_iMaxIterations = 100;

    double d_dResidualNorm    = 0.0;
    double d_dInitialResidual = 0.0;

    double d_dAbsoluteTolerance = 1.0e-14;
    double d_dRelativeTolerance = 1.0e-09;

    int d_iDebugPrintInfoLevel = 0;

    bool d_bUseZeroInitialGuess;

    int d_iObjectId;

    static int d_iInstanceId; // used to differentiate between different instances of the class

    //! keeps track of iteration statistics over solver lifetime
    std::vector<int> d_iterationHistory;

    std::shared_ptr<AMP::Database> d_db = nullptr;

    /** Pointer to global database
     *  This is temporary fix and eventually either d_global_db or d_db should go away
     *  This is introduced to allow for solver factories to access databases in the global
     *  database for the construction of nested solvers
     */
    std::shared_ptr<AMP::Database> d_global_db;

    std::shared_ptr<AMP::Operator::Operator> d_pOperator = nullptr;

    std::shared_ptr<AMP::IO::Writer> d_writer = nullptr;


private:
};

} // namespace AMP::Solver

#endif
