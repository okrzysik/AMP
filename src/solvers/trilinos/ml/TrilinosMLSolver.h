
#ifndef included_AMP_TrilinosMLSolver
#define included_AMP_TrilinosMLSolver

#include "AMP/matrices/trilinos/EpetraMatrix.h"
#include "AMP/solvers/SolverStrategy.h"
#include "AMP/solvers/SolverStrategyParameters.h"
#include "AMP/solvers/trilinos/ml/MLoptions.h"

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-parameter"
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wextra"
#include "ml_MultiLevelPreconditioner.h"
#include "ml_include.h"
#pragma GCC diagnostic pop
#pragma GCC diagnostic pop


namespace AMP {
namespace Solver {


typedef SolverStrategyParameters TrilinosMLSolverParameters;


/**
 * The TrilinosMLSolver is a wrapper to the Trilinos ML solver. ML provides implementations of
 * various algebraic multigrid methods. The wrapper at present simply provides an adaptor
 * to enable AMP users to use the black box ML preconditioner.
 */

class TrilinosMLSolver : public SolverStrategy
{

public:
    /**
     * Default constructor
     */
    TrilinosMLSolver();

    /**
     * Main constructor.
     @param [in] parameters The parameters object contains a database object which must contain the
     following fields in addition to the fields expected by the base class SolverStrategy class:

     1. name:  print_info_level, type: integer, (optional), default value: 0
     acceptable values (non-negative integer values)

     2. name:  prec_type, type: string, (optional), default value: "MGV"
     acceptable values (see ML manual)

     3. name:  max_levels, type: integer, (optional), default value: 5
     acceptable values (see ML manual)

     4. name:  increasingordecreasing, type: string, (optional), default value: "increasing"
     acceptable values (see ML manual)

     5. name:  aggregation_dampingfactor, type: double, (optional), default value: 4.0/3/0
     acceptable values (see ML manual)

     6. name:  aggregationthreshold, type: double, (optional), default value: 0.0
     acceptable values (see ML manual)

     7. name:  eigen-analysis_type, type: string, (optional), default value: "cg"
     acceptable values (see ML manual)

     8. name:  eigen-analysis_iterations, type: integer, (optional), default value: 10
     acceptable values (see ML manual)

     9. name:  smoother_sweeps, type: integer, (optional), default value: 2
     acceptable values (see ML manual)

     10. name:  smoother_dampingfactor, type: double, (optional), default value: 1.0
     acceptable values (see ML manual)

     11. name:  aggregation_nodes_per_aggregate, type: integer, (optional), default value: none
     acceptable values (see ML manual)

     12. name:  aggregation_nextlevel_aggregates_per_process, type: integer, (optional), default
    value: none
     acceptable values (see ML manual)

     13. name:  aggregation_damping_factor, type: , (optional), default value: none
     acceptable values (see ML manual)

     14. name:  energy_minimization_enable, type: , (optional), default value: none
     acceptable values (see ML manual)

     15. name:  smoother_preorpost, type: string, (optional), default value: "both"
     acceptable values (see ML manual)

     16. name:  smoothertype, type: string, (optional), default value: "symmetric Gauss-Seidel"
     acceptable values (see ML manual)

     17. name:  coarse_type, type: string, (optional), default value: "Amesos-KLU"
     acceptable values (see ML manual)

     18. name:  PDE_equations, type: integer, (optional), default value: 1
     acceptable values (see ML manual)

     19. name:  coarse_maxsize, type: integer, (optional), default value: 128
     acceptable values (see ML manual)

     20. name: USE_EPETRA, type: bool, (optional), default value: true

     21. name: problem_type, type: string, (optional), default value: "SA"
     acceptable values "SA" for symmetric and "NSSA" for unsymmetric problems

     22. name: aggregation_aux_enable, type: bool, (optional), default_value: false

     23. name: aggregation_aux_threshold, type: double, (optional), default_values: 0.0

    24. name: null_space_type, type: string, (optional), default value: none

      25. name: null_space_dimension, type: integer, (optional), default value: 3

      26. name: null_space_add_default_vectors, type: bool, (optional), default value: true
      */
    explicit TrilinosMLSolver( std::shared_ptr<TrilinosMLSolverParameters> parameters );

    /**
     * Default destructor
     */
    virtual ~TrilinosMLSolver();

    /**
     * Solve the system \f$Au = f\f$.
     @param [in] f : shared pointer to right hand side vector
     @param [out] u : shared pointer to approximate computed solution
     */
    void solve( std::shared_ptr<const AMP::LinearAlgebra::Vector> f,
                std::shared_ptr<AMP::LinearAlgebra::Vector> u ) override;

    /**
     * Return a shared pointer to the ML_Epetra::MultiLevelPreconditioner object
     */
    inline const std::shared_ptr<ML_Epetra::MultiLevelPreconditioner> getMLSolver( void )
    {
        return d_mlSolver;
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
    void initialize( std::shared_ptr<SolverStrategyParameters> const parameters ) override;

    /**
     * Register the operator that the solver will use during solves
     @param [in] op shared pointer to the linear operator $A$ for equation \f$A u = f\f$
     */
    void registerOperator( const std::shared_ptr<AMP::Operator::Operator> op ) override;

    /**
     * Resets the associated operator internally with new parameters if necessary
     * @param [in] params
     *        OperatorParameters object that is NULL by default
     */
    void resetOperator( const std::shared_ptr<AMP::Operator::OperatorParameters> params ) override;

    /**
     * Resets the solver internally with new parameters if necessary
     * @param [in] params
     *        SolverStrategyParameters object that is NULL by default
     * Currently every call to reset destroys the ML preconditioner object
     * and recreates it based on the parameters object. See constructor for
     * fields required for parameter object.
     */
    void reset( std::shared_ptr<SolverStrategyParameters> params ) override;

protected:
    void reSolveWithLU( std::shared_ptr<const AMP::LinearAlgebra::Vector> f,
                        std::shared_ptr<AMP::LinearAlgebra::Vector> u );

    void getFromInput( std::shared_ptr<AMP::Database> db );

    void convertMLoptionsToTeuchosParameterList();

    void buildML();

    void computeCoordinates( const std::shared_ptr<AMP::Operator::Operator> op );
    void computeNullSpace( const std::shared_ptr<AMP::Operator::Operator> op );

private:
    bool d_bUseEpetra;
    ML *d_ml;
    ML_Aggregate *d_mlAggregate;

    AMP_MPI d_comm;

    bool d_bCreationPhase; /**< set to true if the PC is not ready and false otherwise. */
    bool d_bRobustMode;

    std::shared_ptr<MLoptions> d_mlOptions;

    std::shared_ptr<ML_Epetra::MultiLevelPreconditioner> d_mlSolver;
    std::shared_ptr<AMP::LinearAlgebra::EpetraMatrix> d_matrix;
    Teuchos::ParameterList d_MLParameterList;

    std::vector<double> d_x_values;
    std::vector<double> d_y_values;
    std::vector<double> d_z_values;
    std::vector<double> d_null_space;
};
} // namespace Solver
} // namespace AMP

#endif
