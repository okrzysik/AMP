#ifndef included_AMP_QMRCGSTABSolver
#define included_AMP_QMRCGSTABSolver

#include "solvers/SolverStrategy.h"
#include "utils/AMP_MPI.h"

namespace AMP {
namespace Solver {

/**
 * The QMRCGSTABSolver class implements the QMRCGSTAB method for non-symmetric linear systems
 * introduced by Chan et. al.
 * 
 * The implementation here is mostly based on the MATLAB code at
 * https://link.springer.com/content/pdf/bbm%3A978-3-8348-8100-7%2F1.pdf
 * Currently no preconditioning
 */

class QMRCGSTABSolver : public SolverStrategy
{
public:
    /**
     * default constructor
     */
    QMRCGSTABSolver();

    /**
     * main constructor
     @param [in] parameters The parameters object
     contains a database objects containing the following fields:

     1. type: double, name : relative_tolerance, default value of $1.0e-9$, relative tolerance for
    CG solver
    acceptable values (non-negative real values)

     2. type: bool, name : uses_preconditioner, default value false
        acceptable values (false, true),
        side effect: if false sets string pc_type to "none"

     3. type: string, name : pc_side, default value "RIGHT",
     acceptable values ("RIGHT", "LEFT", "SYMMETRIC" )
         active only when uses_preconditioner set to true
     */
    explicit QMRCGSTABSolver( AMP::shared_ptr<SolverStrategyParameters> parameters );

    /** 
     * static create routine that is used by SolverFactory 
     @param [in] parameters The parameters object
     contains a database objects with the fields listed for the constructor above
     */
    static AMP::shared_ptr<SolverStrategy> createSolver( AMP::shared_ptr<SolverStrategyParameters> solverStrategyParameters )
    {
      return AMP::make_shared<QMRCGSTABSolver> ( solverStrategyParameters );
    }

    /**
     * Default destructor
     */
    virtual ~QMRCGSTABSolver();

    /**
     * Solve the system \f$Au = 0\f$.
     * @param [in] f : shared pointer to right hand side vector
     * @param [out] u : shared pointer to approximate computed solution
     */
    void solve( AMP::shared_ptr<const AMP::LinearAlgebra::Vector> f,
                AMP::shared_ptr<AMP::LinearAlgebra::Vector>
                    u ) override;

    /**
     * Initialize the QMRCGSTABSolver. Should not be necessary for the user to call in general.
     * @param parameters
     */
    void initialize( AMP::shared_ptr<SolverStrategyParameters> const parameters ) override;

    /**
     * returns a shared pointer to a preconditioner object. The preconditioner is derived from
     * a SolverStrategy class
     */
    inline AMP::shared_ptr<AMP::Solver::SolverStrategy> getPreconditioner( void )
    {
        return d_pPreconditioner;
    }

    /**
     * sets a shared pointer to a preconditioner object. The preconditioner is derived from
     * a SolverStrategy class
     * @param pc shared pointer to preconditioner
     */
    inline void setPreconditioner( AMP::shared_ptr<AMP::Solver::SolverStrategy> pc )
    {
        d_pPreconditioner = pc;
    }

    /**
     * Register the operator that the solver will use during solves
     * @param [in] op shared pointer to operator $A()$ for equation \f$A(u) = f\f$
     */
    void registerOperator( const AMP::shared_ptr<AMP::Operator::Operator> op ) override;

    /**
     * Resets the registered operator internally with new parameters if necessary
     * @param parameters    OperatorParameters object that is NULL by default
     */
    void resetOperator( const AMP::shared_ptr<AMP::Operator::OperatorParameters> parameters ) override;

protected:
    void getFromInput( const AMP::shared_ptr<AMP::Database> &db );

private:
    AMP_MPI d_comm;

    double d_dRelativeTolerance = 1.0e-10; //! relative tolerance to converge to

    int d_restarts = 0; //! number of times the solver is restarted

    bool d_bUsesPreconditioner = false;

    std::string d_preconditioner_side;
    
    AMP::shared_ptr<AMP::Solver::SolverStrategy> d_pPreconditioner;
};
}
}

#endif
