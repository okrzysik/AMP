#include "AMP/operators/LinearOperator.h"
#include "AMP/solvers/DiagonalSolver.h"
#include "AMP/solvers/SolverFactory.h"
#include "ProfilerApp.h"

namespace AMP::Solver {

/****************************************************************
 *  Constructors                                                 *
 ****************************************************************/

template<typename T>
DiagonalSolver<T>::DiagonalSolver( std::shared_ptr<AMP::Solver::SolverStrategyParameters> parameters )
    : SolverStrategy( parameters )
{
    AMP_ASSERT( parameters );

    // Initialize
    initialize( parameters );
}

/****************************************************************
 *  Initialize                                                   *
 ****************************************************************/
template<typename T>
void DiagonalSolver<T>::initialize(
    std::shared_ptr<const AMP::Solver::SolverStrategyParameters> parameters )
{
    AMP_ASSERT( parameters );

    auto db = parameters->d_db;
    getFromInput( db );

    if ( parameters->d_pNestedSolver ) {
        d_pNestedSolver = parameters->d_pNestedSolver;
    } else {
        if ( d_bUsesNestedSolver ) {
            auto nestedName =
                db->getWithDefault<std::string>( "nested_solver_name", "NestedSolver" );
            auto outerDB = db->keyExists( nestedName ) ? db : parameters->d_global_db;
            if ( outerDB ) {
                auto nestedDB = outerDB->getDatabase( nestedName );
                auto innerParameters =
                    std::make_shared<AMP::Solver::SolverStrategyParameters>( nestedDB );
                innerParameters->d_global_db = parameters->d_global_db;
                innerParameters->d_pOperator = d_pOperator;
                d_pNestedSolver = AMP::Solver::SolverFactory::create( innerParameters );
                AMP_ASSERT( d_pNestedSolver );
            }
        }
    }
}

// Function to get values from input
template<typename T>
void DiagonalSolver<T>::getFromInput( std::shared_ptr<AMP::Database> db )
{
    d_dDivergenceTolerance = db->getWithDefault<T>( "divergence_tolerance", 1.0e+03 );
    d_bUsesNestedSolver    = db->getWithDefault<bool>( "uses_nested_solver", false );
}

/****************************************************************
 *  Solve                                                        *
 ****************************************************************/
template<typename T>
void DiagonalSolver<T>::apply( std::shared_ptr<const AMP::LinearAlgebra::Vector> f,
                               std::shared_ptr<AMP::LinearAlgebra::Vector> u )
{
    PROFILE( "DiagonalSolver<T>::apply" );

    // Always zero before checking stopping criteria for any reason
    d_iNumberIterations = 1;

    // Check input vector states
    AMP_ASSERT( ( u->getUpdateStatus() == AMP::LinearAlgebra::UpdateState::UNCHANGED ) ||
                ( u->getUpdateStatus() == AMP::LinearAlgebra::UpdateState::LOCAL_CHANGED ) );

    {
        PROFILE( "DiagonalSolver<T>:: u->makeConsistent" );
        u->makeConsistent( AMP::LinearAlgebra::ScatterType::CONSISTENT_SET );
    }
    {
        PROFILE( "DiagonalSolver<T>:: r = f-Au (final)" );
        d_pOperator->residual( f, u, r );
    }
    {
        PROFILE( "DiagonalSolver<T>:: r->L2Norm (final)" );
        d_dResidualNorm = static_cast<T>( r->L2Norm() );
    }
    // final check updates flags if needed
    checkStoppingCriteria( d_dResidualNorm );

    if ( d_iDebugPrintInfoLevel > 0 ) {
        AMP::pout << "CG: final residual: " << d_dResidualNorm
                  << ", iterations: " << d_iNumberIterations << ", convergence reason: "
                  << SolverStrategy::statusToString( d_ConvergenceStatus ) << std::endl;
    }
    if ( d_iDebugPrintInfoLevel > 1 ) {
        AMP::pout << "CG: final L2Norm of solution: " << u->L2Norm() << std::endl;
    }
}

template<typename T>
void DiagonalSolver<T>::resetOperator(
    std::shared_ptr<const AMP::Operator::OperatorParameters> params )
{
    if ( d_pOperator ) {
        d_pOperator->reset( params );
    }

    // should add a mechanism for the linear operator to provide updated parameters for the
    // preconditioner operator
    // though it's unclear where this might be necessary
    if ( d_pNestedSolver ) {
        d_pNestedSolver->resetOperator( params );
    }
}
} // namespace AMP::Solver
