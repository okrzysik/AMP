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
    AMP_ERROR("Not implemented");
}

void BoomerAMGSolver::createHYPREMatrix( const AMP::shared_ptr<AMP::LinearAlgebra::Matrix> matrix )
{
    AMP_ERROR("Not implemented");
    int ierr;

    const auto myFirstRow = matrix->getLeftDOFManager()->beginDOF();
    const auto myEndRow   = matrix->getLeftDOFManager()->endDOF(); // check whether endDOF is truly the last -1 

    ierr = HYPRE_IJMatrixCreate( d_comm.getCommunicator(), myFirstRow, myEndRow-1, myFirstRow, myEndRow-1, &d_ijMatrix );
    ierr = HYPRE_IJMatrixSetObjectType( d_ijMatrix, HYPRE_PARCSR);
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
