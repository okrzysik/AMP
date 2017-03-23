#include "solvers/hypre/BoomerAMGSolver.h"

#include "ProfilerApp.h"
#include "matrices/Matrix.h"
#include "operators/LinearOperator.h"
#include "utils/Utilities.h"
#include "vectors/DataChangeFirer.h"

#include <iomanip>

namespace AMP {
namespace Solver {


/****************************************************************
* Constructors / Destructor                                     *
****************************************************************/
BoomerAMGSolver::BoomerAMGSolver()
{
    d_bCreationPhase = true;
}
BoomerAMGSolver::BoomerAMGSolver( AMP::shared_ptr<SolverStrategyParameters> parameters )
    : SolverStrategy( parameters )
{
    AMP_ERROR("Not implemented");
    AMP_ASSERT( parameters.get() != nullptr );
    initialize( parameters );
}
BoomerAMGSolver::~BoomerAMGSolver()
{
    AMP_ERROR("Not implemented");
}

void BoomerAMGSolver::initialize( AMP::shared_ptr<SolverStrategyParameters> const parameters )
{
    AMP_ERROR("Not implemented");
    getFromInput( parameters->d_db );
    if ( d_pOperator.get() != nullptr ) {
        registerOperator( d_pOperator );
    }
}

void BoomerAMGSolver::getFromInput( const AMP::shared_ptr<AMP::Database> &db )
{
    d_num_functions    = db->getIntegerWithDefault( "num_functions", 1);
    d_max_coarse_size  = db->getIntegerWithDefault( "max_coarse_size", 800 );
    d_min_coarse_size  = db->getIntegerWithDefault( "min_coarse_size", 100 );
    d_max_levels       = db->getIntegerWithDefault( "max_levels", 10);

    // 6.2.18 in hypre 11.2 manual 
    if( db->keyExists( "strong_threshold" ) )
        d_strong_threshold = db->getDouble( "strong_threshold" );

    // 6.2.20 in hypre 11.2 manual 
    if( db->keyExists( "max_row_sum" ) )
        d_max_row_sum = db->getDouble( "max_row_sum" );

    // 6.2.21 in hypre 11.2 manual 
    if( db->keyExists( "coarsen_type" ) )
        d_coarsen_type = db->getInteger( "coarsen_type" );

    // 6.2.23 in hypre 11.2 manual 
    if( db->keyExists( "non_galerkin_tol" ) )
        d_non_galerkin_tol = db->getDouble( "non_galerkin_tol" );

    // 6.2.24 in hypre 11.2 manual 
    if( db->keyExists( "measure_type" ) )
        d_measure_type = db->getInteger( "measure_type" );

    // 6.2.25 in hypre 11.2 manual 
    if( db->keyExists( "agg_num_levels" ) )
        d_agg_num_levels = db->getInteger( "agg_num_levels" );

    // 6.2.26 in hypre 11.2 manual 
    if( db->keyExists( "num_paths" ) )
        d_num_paths = db->getInteger( "num_paths" );

    // 6.2.27 in hypre 11.2 manual 
    if( db->keyExists( "cgc_iterations" ) )
        d_cgc_iterations = db->getInteger( "cgc_iterations" );

    // 6.2.28 in hypre 11.2 manual 
    if( db->keyExists( "nodal" ) )
        d_nodal = db->getInteger( "nodal" );

    // 6.2.29 in hypre 11.2 manual 
    if( db->keyExists( "nodal_diag" ) )
        d_nodal_diag = db->getInteger( "nodal_diag" );

    // 6.2.30 in hypre 11.2 manual 
    if( db->keyExists( "interp_type" ) )
        d_interp_type = db->getInteger( "interp_type" );

    // 6.2.31 in hypre 11.2 manual 
    if( db->keyExists( "trunc_factor" ) )
        d_trunc_factor = db->getDouble( "trunc_factor" );

    // 6.2.32 in hypre 11.2 manual 
    if( db->keyExists( "P_max_elements" ) )
        d_P_max_elements = db->getInteger( "P_max_elements" );

    // 6.2.33 in hypre 11.2 manual 
    if( db->keyExists( "separate_weights" ) )
        d_separate_weights = db->getInteger( "separate_weights" );

    // 6.2.34 in hypre 11.2 manual 
    if( db->keyExists( "agg_interp_type" ) )
        d_agg_interp_type = db->getInteger( "agg_interp_type" );

    // 6.2.35 in hypre 11.2 manual 
    if( db->keyExists( "agg_trunc_factor" ) )
        d_agg_trunc_factor = db->getDouble( "agg_trunc_factor" );

    // 6.2.36 in hypre 11.2 manual 
    if( db->keyExists( "agg_P12_trunc_factor" ) )
        d_agg_P12_trunc_factor = db->getDouble( "agg_P12_trunc_factor" );

    // 6.2.37 in hypre 11.2 manual 
    if( db->keyExists( "agg_P_max_elements" ) )
        d_agg_P_max_elements = db->getInteger( "agg_P_max_elements" );

    // 6.2.38 in hypre 11.2 manual 
    if( db->keyExists( "agg_P12_max_elements" ) )
        d_agg_P12_max_elements = db->getInteger( "agg_P12_max_elements" );

    // 6.2.45 in hypre 11.2 manual 
    if( db->keyExists( "cycle_type" ) )
        d_cycle_type = db->getInteger( "cycle_type" );

    // 6.2.46 in hypre 11.2 manual 
    if( db->keyExists( "additive_level" ) )
        d_additive_level = db->getInteger( "additive_level" );

    // 6.2.47 in hypre 11.2 manual 
    if( db->keyExists( "mult_additive_level" ) )
        d_mult_additive_level = db->getInteger( "mult_additive_level" );

    // 6.2.48 in hypre 11.2 manual 
    if( db->keyExists( "simple_level" ) )
        d_simple_level = db->getInteger( "simple_level" );

    // 6.2. in hypre 11.2 manual 
    if( db->keyExists( "additive_trunc_factor" ) )
        d_additive_trunc_factor = db->getDouble( "additive_trunc_factor" );

    // 6.2. in hypre 11.2 manual 
    if( db->keyExists( "number_sweeps" ) )
        d_number_sweeps = db->getInteger( "number_sweeps" );

    // 6.2. in hypre 11.2 manual 
    if( db->keyExists( "relax_type" ) )
        d_relax_type = db->getInteger( "relax_type" );

    // 6.2. in hypre 11.2 manual 
    if( db->keyExists( "relax_order" ) )
        d_relax_order = db->getInteger( "relax_order" );

    // 6.2. in hypre 11.2 manual 
    if( db->keyExists( "relax_weight" ) )
        d_relax_weight = db->getDouble( "relax_weight" );

    // 6.2. in hypre 11.2 manual 
    if( db->keyExists( "chebyshev_order" ) )
        d_chebyshev_order = db->getInteger( "chebyshev_order" );

    // 6.2. in hypre 11.2 manual 
    if( db->keyExists( "chebyshev_fraction" ) )
        d_chebyshev_fraction = db->getDouble( "chebyshev_fraction" );

    // 6.2. in hypre 11.2 manual 
    if( db->keyExists( "smooth_type" ) )
        d_smooth_type = db->getInteger( "smooth_type" );

    // 6.2. in hypre 11.2 manual 
    if( db->keyExists( "smooth_number_levels" ) )
        d_smooth_number_levels = db->getInteger( "smooth_number_levels" );

    // 6.2. in hypre 11.2 manual 
    if( db->keyExists( "smooth_number_sweeps" ) )
        d_smooth_number_sweeps = db->getInteger( "smooth_number_sweeps" );

    // 6.2. in hypre 11.2 manual 
    if( db->keyExists( "schwarz_variant" ) )
        d_schwarz_variant = db->getInteger( "schwarz_variant" );

    // 6.2. in hypre 11.2 manual 
    if( db->keyExists( "schwarz_overlap" ) )
        d_schwarz_overlap = db->getInteger( "schwarz_overlap" );

    // 6.2. in hypre 11.2 manual 
    if( db->keyExists( "schwarz_domain_type" ) )
        d_schwarz_domain_type = db->getInteger( "schwarz_domain_type" );

    // 6.2. in hypre 11.2 manual 
    if( db->keyExists( "logging" ) )
        d_logging = db->getInteger( "logging" );

    // 6.2. in hypre 11.2 manual 
    if( db->keyExists( "debug_flag" ) )
        d_debug_flag = db->getInteger( "debug_flag" );

    // 6.2. in hypre 11.2 manual 
    d_rap2 = db->getIntegerWithDefault( "rap2", 0 );

    // 6.2. in hypre 11.2 manual 
    if( db->keyExists( "keep_transpose" ) )
        d_keep_transpose = db->getInteger( "keep_transpose" );


    AMP_ERROR("Not implemented");
}

void BoomerAMGSolver::createHYPREMatrix( const AMP::shared_ptr<AMP::LinearAlgebra::Matrix> matrix )
{
    AMP_ERROR("Not implemented");
    int ierr;

    const auto myFirstRow = matrix->getLeftDOFManager()->beginDOF();
    const auto myEndRow   = matrix->getLeftDOFManager()->endDOF(); // check whether endDOF is truly the last -1 

    ierr = HYPRE_IJMatrixCreate( d_comm.getCommunicator(), myFirstRow, myEndRow-1, myFirstRow, myEndRow-1, &d_ijMatrix );
    ierr = HYPRE_IJMatrixSetObjectType( d_ijMatrix, HYPRE_PARCSR );
    ierr = HYPRE_IJMatrixInitialize( d_ijMatrix );

    std::vector<unsigned int> cols;
    std::vector<double> values;

    // iterate over all rows
    for(auto i=myFirstRow; i!=myEndRow; ++i) {
        matrix->getRowByGlobalID(i, cols, values);
        const int nrows = 1;
        const auto irow = i;
        const auto ncols = cols.size();
        ierr = HYPRE_IJMatrixSetValues( d_ijMatrix, 
                                        nrows, 
                                        (HYPRE_Int *)&ncols, 
                                        (HYPRE_Int *)&irow, 
                                        (HYPRE_Int *) &cols[0], 
                                        (const double *) &values[0] );
    }

    ierr = HYPRE_IJMatrixAssemble( d_ijMatrix );
}

void BoomerAMGSolver::registerOperator( const AMP::shared_ptr<AMP::Operator::Operator> op )
{
    AMP_ERROR("Not implemented");
    d_pOperator = op;
    AMP_INSIST( d_pOperator.get() != nullptr,
                "ERROR: BoomerAMGSolver::registerOperator() operator cannot be NULL" );

    auto linearOperator =
        AMP::dynamic_pointer_cast<AMP::Operator::LinearOperator>( d_pOperator );
    AMP_INSIST( linearOperator.get() != nullptr, "linearOperator cannot be NULL" );

    auto matrix = linearOperator->getMatrix();
    AMP_INSIST( matrix.get() != nullptr, "matrix cannot be NULL" );

    createHYPREMatrix( matrix );

    // the next section of code should initialize a hypre IJ matrix based on the AMP matrix
    d_bCreationPhase = false;
}


void BoomerAMGSolver::resetOperator(
    const AMP::shared_ptr<AMP::Operator::OperatorParameters> params )
{
    AMP_ERROR("Not implemented");
    PROFILE_START( "resetOperator" );
    AMP_INSIST( ( d_pOperator.get() != nullptr ),
                "ERROR: BoomerAMGSolver::resetOperator() operator cannot be NULL" );
    d_pOperator->reset( params );
    reset( AMP::shared_ptr<SolverStrategyParameters>() );
    PROFILE_STOP( "resetOperator" );
}


void BoomerAMGSolver::reset( AMP::shared_ptr<SolverStrategyParameters> )
{
    AMP_ERROR("Not implemented");
    PROFILE_START( "reset" );
    registerOperator( d_pOperator );
    PROFILE_STOP( "reset" );
}


void BoomerAMGSolver::solve( AMP::shared_ptr<const AMP::LinearAlgebra::Vector> f,
                                 AMP::shared_ptr<AMP::LinearAlgebra::Vector> u )
{
    AMP_ERROR("Not implemented");
    PROFILE_START( "solve" );
    // in this case we make the assumption we can access a EpetraMat for now
    AMP_INSIST( d_pOperator.get() != nullptr,
                "ERROR: BoomerAMGSolver::solve() operator cannot be NULL" );

    if ( d_bUseZeroInitialGuess ) {
        u->zero();
    }

    if ( d_bCreationPhase ) {
        d_bCreationPhase = false;
    }

    AMP::shared_ptr<AMP::LinearAlgebra::Vector> r;

    bool computeResidual = false;

    double initialResNorm = 0., finalResNorm = 0.;

    if ( computeResidual ) {
        r = f->cloneVector();
        d_pOperator->residual( f, u, r );
        initialResNorm = r->L2Norm();

        if ( d_iDebugPrintInfoLevel > 1 ) {
            AMP::pout << "BoomerAMGSolver::solve(), L2 norm of residual before solve "
                      << std::setprecision( 15 ) << initialResNorm << std::endl;
        }
    }

    if ( d_iDebugPrintInfoLevel > 2 ) {
        double solution_norm = u->L2Norm();
        AMP::pout << "BoomerAMGSolver : before solve solution norm: " << std::setprecision( 15 )
                  << solution_norm << std::endl;
    }

    // add in code for solve here

    // Check for NaNs in the solution (no communication necessary)
    double localNorm = u->localL2Norm();
    AMP_INSIST( localNorm == localNorm, "NaNs detected in solution" );

    // we are forced to update the state of u here
    // as Hypre is not going to change the state of a managed vector
    // an example where this will and has caused problems is when the
    // vector is a petsc managed vector being passed back to PETSc
    if ( u->isA<AMP::LinearAlgebra::DataChangeFirer>() ) {
        u->castTo<AMP::LinearAlgebra::DataChangeFirer>().fireDataChange();
    }

    if ( d_iDebugPrintInfoLevel > 2 ) {
        double solution_norm = u->L2Norm();
        AMP::pout << "BoomerAMGSolver : after solve solution norm: " << std::setprecision( 15 )
                  << solution_norm << std::endl;
    }

    if ( computeResidual ) {
        d_pOperator->residual( f, u, r );
        finalResNorm = r->L2Norm();

        if ( d_iDebugPrintInfoLevel > 1 ) {
            AMP::pout << "BoomerAMGSolver::solve(), L2 norm of residual after solve "
                      << std::setprecision( 15 ) << finalResNorm << std::endl;
        }
    }

    PROFILE_STOP( "solve" );
}

} // Solver
} // AMP
