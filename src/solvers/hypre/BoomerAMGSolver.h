#ifndef included_AMP_BoomerAMGSolver
#define included_AMP_BoomerAMGSolver

#include "solvers/SolverStrategy.h"
#include "solvers/SolverStrategyParameters.h"

extern "C" {
#include "HYPRE.h"
#include "HYPRE_IJ_mv.h"
#include "HYPRE_parcsr_ls.h"
}

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
                AMP::shared_ptr<AMP::LinearAlgebra::Vector> u ) override;

    /**
     * Initialize the solution vector and potentially create internal vectors needed for solution
     @param [in] parameters The parameters object
     contains a database object. Refer to the documentation for the constructor to see what fields
     are required.
     This routine assumes that a non-NULL operator of type LinearOperator has been registered with
     the solver.
     The LinearOperator currently is assumed to contain a pointer to an EpetraMatrix object.
     */
    void initialize( AMP::shared_ptr<SolverStrategyParameters> const parameters ) override;

    /**
     * Register the operator that the solver will use during solves
     @param [in] op shared pointer to the linear operator $A$ for equation \f$A u = f\f$
     */
    void registerOperator( const AMP::shared_ptr<AMP::Operator::Operator> op ) override;

    /**
     * Resets the associated operator internally with new parameters if necessary
     * @param [in] params
     *        OperatorParameters object that is NULL by default
     */
    void resetOperator( const AMP::shared_ptr<AMP::Operator::OperatorParameters> params ) override;

    /**
     * Resets the solver internally with new parameters if necessary
     * @param [in] params
     *        SolverStrategyParameters object that is NULL by default
     * Currently every call to reset destroys the BoomerAMG preconditioner object
     * and recreates it based on the parameters object. See constructor for
     * fields required for parameter object.
     */
    void reset( AMP::shared_ptr<SolverStrategyParameters> params ) override;

    void getFromInput( const AMP::shared_ptr<AMP::Database> &db );

private:

    /**
     * create the internal HYPRE_IJMatrix based on the AMP matrix
     */
    void createHYPREMatrix( const AMP::shared_ptr<AMP::LinearAlgebra::Matrix> matrix );

    /** 
     * create and initialize the internal hypre vectors for rhs and solution
     */
    void createHYPREVectors( );

    /**
     *  copy values from amp vector to hypre vector
     */
    void copyToHypre( AMP::shared_ptr<const AMP::LinearAlgebra::Vector> amp_v, 
                      HYPRE_IJVector hypre_v );

    /**
     *  copy values from hypre vector to amp vector
     */
    void copyFromHypre( HYPRE_IJVector hypre_v, 
                        AMP::shared_ptr<AMP::LinearAlgebra::Vector> amp_v );


    void setParameters(void); //! set BoomerAMG parameters based on internally set variables

    AMP_MPI d_comm;

    HYPRE_IJMatrix d_ijMatrix;  //! pointer to HYPRE matrix struct

    HYPRE_IJVector d_hypre_rhs;       //! pointer to HYPRE representation of rhs 
    HYPRE_IJVector d_hypre_sol;       //! pointer to HYPRE representation of solution

    HYPRE_Solver d_solver;      //! pointer to HYPRE BoomerAMG solver

    int d_num_functions              = 1;
    int d_min_iterations             = 0;
    int d_max_coarse_size            = 800;
    int d_min_coarse_size            = 100;
    int d_max_levels                 = 10;
    int d_coarsen_type;
    int d_measure_type;
    int d_agg_num_levels;
    int d_num_paths;
    int d_cgc_iterations;
    int d_nodal;
    int d_nodal_diag;
    int d_interp_type;
    int d_P_max_elements;
    int d_separate_weights;
    int d_agg_interp_type;
    int d_agg_P_max_elements;
    int d_agg_P12_max_elements;
    int d_number_samples; 
    int d_cycle_type;
    int d_additive_level;
    int d_mult_additive_level;
    int d_simple_level;
    int d_add_P_max_elmts;
    int d_number_sweeps;
    int d_relax_type;
    int d_relax_order;
    int d_chebyshev_order;
    int d_smooth_type;
    int d_smooth_number_levels;
    int d_smooth_number_sweeps;
    int d_schwarz_variant;
    int d_schwarz_overlap;
    int d_schwarz_domain_type;
    int d_schwarz_nonsymmetric;
    int d_logging;
    int d_debug_flag;
    int d_rap2;
    int d_keep_transpose;

    double d_strong_threshold;
    double d_max_row_sum;
    double d_non_galerkin_tol;
    double d_trunc_factor;
    double d_agg_trunc_factor;
    double d_agg_P12_trunc_factor;
    double d_additive_trunc_factor;
    double d_relax_weight;
    double d_outer_weight;
    double d_chebyshev_fraction;
    double d_schwarz_weight;

    bool d_bCreationPhase; /**< set to true if the PC is not ready and false otherwise. */
};
}
}

#endif
