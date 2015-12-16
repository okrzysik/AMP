
#include "solvers/PelletStackMechanicsSolver.h"

namespace AMP {
namespace Solver {


PelletStackMechanicsSolver::PelletStackMechanicsSolver(
    AMP::shared_ptr<PelletStackMechanicsSolverParameters> params )
    : SolverStrategy( params )
{
    d_columnSolver  = params->d_columnSolver;
    d_pelletStackOp = AMP::dynamic_pointer_cast<AMP::Operator::PelletStackOperator>( d_pOperator );
}


void PelletStackMechanicsSolver::resetOperator(
    const AMP::shared_ptr<AMP::Operator::OperatorParameters> params )
{
    d_columnSolver->resetOperator( params );
}


void PelletStackMechanicsSolver::solve( AMP::shared_ptr<const AMP::LinearAlgebra::Vector> f,
                                        AMP::shared_ptr<AMP::LinearAlgebra::Vector>
                                            u )
{
    AMP_ASSERT( ( f->getUpdateStatus() == AMP::LinearAlgebra::Vector::UNCHANGED ) ||
                ( f->getUpdateStatus() == AMP::LinearAlgebra::Vector::LOCAL_CHANGED ) );
    AMP::shared_ptr<const AMP::LinearAlgebra::Vector> fInternal = f;
    if ( d_pelletStackOp->useScaling() ) {
        if ( d_fbuffer1 == NULL ) {
            d_fbuffer1 = f->cloneVector();
        }
        d_fbuffer1->copyVector( f );
        d_pelletStackOp->applyUnscaling( d_fbuffer1 );
        fInternal = d_fbuffer1;
    }
    if ( d_pelletStackOp->useSerial() ) {
        solveSerial( fInternal, u );
    } else {
        solveScan( fInternal, u );
    }
}

void PelletStackMechanicsSolver::solveSerial( AMP::shared_ptr<const AMP::LinearAlgebra::Vector> f,
                                              AMP::shared_ptr<AMP::LinearAlgebra::Vector>
                                                  u )
{
    if ( d_fbuffer2 == NULL ) {
        d_fbuffer2 = f->cloneVector();
    }

    unsigned int totalNumberOfPellets = d_pelletStackOp->getTotalNumberOfPellets();

    int locPellIdx = d_pelletStackOp->getLocalIndexForPellet( 0 );

    if ( locPellIdx != -1 ) {
        AMP::shared_ptr<AMP::Solver::SolverStrategy> currSolver =
            d_columnSolver->getSolver( locPellIdx );
        AMP::shared_ptr<AMP::Operator::Operator> currOp = currSolver->getOperator();
        AMP::LinearAlgebra::Vector::shared_ptr subUvec  = currOp->subsetInputVector( u );
        AMP_ASSERT( subUvec != NULL );
        AMP::LinearAlgebra::Vector::const_shared_ptr subFvec = currOp->subsetOutputVector( f );
        AMP_ASSERT( subFvec != NULL );
        currSolver->solve( subFvec, subUvec );
    }

    for ( unsigned int pellId = 1; pellId < totalNumberOfPellets; pellId++ ) {
        AMP::shared_ptr<AMP::Database> emptyDb;
        AMP::shared_ptr<AMP::Operator::PelletStackOperatorParameters> pelletStackOpParams(
            new AMP::Operator::PelletStackOperatorParameters( emptyDb ) );
        pelletStackOpParams->d_currentPellet = pellId;
        d_pelletStackOp->reset( pelletStackOpParams );
        d_pelletStackOp->residual( f, u, d_fbuffer2 );
        locPellIdx = d_pelletStackOp->getLocalIndexForPellet( pellId );
        if ( locPellIdx != -1 ) {
            AMP::shared_ptr<AMP::Solver::SolverStrategy> currSolver =
                d_columnSolver->getSolver( locPellIdx );
            AMP::shared_ptr<AMP::Operator::Operator> currOp = currSolver->getOperator();
            AMP::LinearAlgebra::Vector::shared_ptr subUvec  = currOp->subsetInputVector( u );
            AMP_ASSERT( subUvec != NULL );
            AMP::LinearAlgebra::Vector::shared_ptr subFvec =
                currOp->subsetOutputVector( d_fbuffer2 );
            AMP_ASSERT( subFvec != NULL );
            currSolver->solve( subFvec, subUvec );
        }
    } // end for pellId
}

void PelletStackMechanicsSolver::solveScan( AMP::shared_ptr<const AMP::LinearAlgebra::Vector> f,
                                            AMP::shared_ptr<AMP::LinearAlgebra::Vector>
                                                u )
{
    d_columnSolver->solve( f, u );
    if ( d_pelletStackOp->onlyZcorrection() ) {
        AMP::LinearAlgebra::Vector::shared_ptr nullVec;
        d_pelletStackOp->apply( nullVec, u );
    } else {
        if ( d_fbuffer2 == NULL ) {
            d_fbuffer2 = f->cloneVector();
        }
        d_pelletStackOp->residual( f, u, d_fbuffer2 );
        d_columnSolver->solve( d_fbuffer2, u );
    }
}


} // Solver
} // AMP
