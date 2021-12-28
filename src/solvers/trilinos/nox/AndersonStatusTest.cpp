#include "AMP/solvers/trilinos/nox/AndersonStatusTest.h"

#include "AMP/solvers/trilinos/thyra/TrilinosThyraModelEvaluator.h"
#include "AMP/vectors/trilinos/thyra/ThyraVector.h"
#include "AMP/vectors/trilinos/thyra/ThyraVectorWrapper.h"

// Trilinos includes
DISABLE_WARNINGS
#include "NOX_Thyra.H"
#include "NOX_Thyra_Group.H"
ENABLE_WARNINGS


namespace AMP::Solver {


/****************************************************************
 *  Constructors                                                 *
 ****************************************************************/
AndersonStatusTest::AndersonStatusTest( std::shared_ptr<AMP::Database> db )
{
    if ( db->keyExists( "AndersonConvergenceVariables" ) ) {
        d_variableNames = db->getVector<std::string>( "AndersonConvergenceVariables" );
        if ( db->keyExists( "AndersonConvergenceTolerances" ) ) {
            d_tolerances = db->getVector<double>( "AndersonConvergenceTolerances" );
        }
    }
    AMP_ASSERT( d_variableNames.size() == d_tolerances.size() );
    d_relativeResiduals.resize( d_variableNames.size(), 1e12 );
}
AndersonStatusTest::~AndersonStatusTest() = default;

NOX::StatusTest::StatusType AndersonStatusTest::checkStatus( const NOX::Solver::Generic &solver,
                                                             NOX::StatusTest::CheckType checkType )
{
    // Check for early exit
    if ( checkType == NOX::StatusTest::None )
        return NOX::StatusTest::Unevaluated;

    // Get the current and previous solutions from solver
    auto curSolVec  = solver.getSolutionGroup().getXPtr();
    auto prevSolVec = solver.getPreviousSolutionGroup().getXPtr();
    AMP_ASSERT( curSolVec != Teuchos::null );
    AMP_ASSERT( prevSolVec != Teuchos::null );

    // Cast vectors to Thyra::Vector base class
    auto curSolThyraVec  = Teuchos::rcp_dynamic_cast<const NOX::Thyra::Vector>( curSolVec );
    auto prevSolThyraVec = Teuchos::rcp_dynamic_cast<const NOX::Thyra::Vector>( prevSolVec );
    AMP_ASSERT( curSolThyraVec != Teuchos::null );
    AMP_ASSERT( prevSolThyraVec != Teuchos::null );

    // Cast to Thyra wrappers
    auto curSolWrappedVec = dynamic_cast<const AMP::LinearAlgebra::ThyraVectorWrapper *>(
        curSolThyraVec->getThyraRCPVector().getRawPtr() );
    auto prevSolWrappedVec = dynamic_cast<const AMP::LinearAlgebra::ThyraVectorWrapper *>(
        prevSolThyraVec->getThyraRCPVector().getRawPtr() );
    AMP_ASSERT( curSolWrappedVec );
    AMP_ASSERT( prevSolWrappedVec );

    // Finally get the AMP vectors from wrapper
    AMP_ASSERT( curSolWrappedVec->numVecs() == 1 );
    AMP_ASSERT( prevSolWrappedVec->numVecs() == 1 );
    AMP::LinearAlgebra::Vector::const_shared_ptr curSolAmpVec  = curSolWrappedVec->getVec( 0 );
    AMP::LinearAlgebra::Vector::const_shared_ptr prevSolAmpVec = prevSolWrappedVec->getVec( 0 );
    AMP_ASSERT( curSolAmpVec );
    AMP_ASSERT( prevSolAmpVec );

    // Process each variable in list
    bool converged = true;
    for ( size_t i = 0; i < d_variableNames.size(); ++i ) {
        auto thisVar     = std::make_shared<AMP::LinearAlgebra::Variable>( d_variableNames[i] );
        auto thisCurVec  = curSolAmpVec->subsetVectorForVariable( thisVar );
        auto thisPrevVec = prevSolAmpVec->subsetVectorForVariable( thisVar );
        if ( thisCurVec ) {
            AMP_ASSERT( thisPrevVec );
            auto thisDiffVec = thisCurVec->cloneVector();
            thisDiffVec->subtract( *thisCurVec, *thisPrevVec );
            d_relativeResiduals[i] =
                double( thisDiffVec->L2Norm() ) / double( thisCurVec->L2Norm() );
            if ( d_relativeResiduals[i] > d_tolerances[i] )
                converged = false;
        }
    }

    // Perform reduction on convergence
    converged = curSolAmpVec->getComm().allReduce( converged );
    d_status  = NOX::StatusTest::Unconverged;
    if ( converged )
        d_status = NOX::StatusTest::Converged;

    return d_status;
}

NOX::StatusTest::StatusType AndersonStatusTest::getStatus() const { return d_status; }

std::ostream &AndersonStatusTest::print( std::ostream &stream, int indent ) const
{
    AMP_ASSERT( indent >= 0 );

    for ( size_t i = 0; i < d_variableNames.size(); ++i ) {
        for ( int j = 0; j < indent; j++ )
            stream << ' ';
        stream << "Relative update on " << d_variableNames[i] << " = ";
        if ( d_status == NOX::StatusTest::Unevaluated ) {
            stream << "?";
        } else {
            stream << std::scientific << d_relativeResiduals[i];
        }
        stream << " < " << std::scientific << d_tolerances[i] << std::endl;
    }

    return stream;
}
} // namespace AMP::Solver
