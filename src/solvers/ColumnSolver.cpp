
#include "ColumnSolver.h"
#include "AMP/operators/ColumnOperatorParameters.h"

namespace AMP {
namespace Solver {


ColumnSolver::ColumnSolver( std::shared_ptr<const SolverStrategyParameters> parameters )
    : SolverStrategy( parameters )
{
    AMP_ASSERT( parameters );
    std::shared_ptr<AMP::Database> db = parameters->d_db;
    d_IterationType       = db->getWithDefault<std::string>( "IterationType", "GaussSeidel" );
    d_resetColumnOperator = db->getWithDefault( "ResetColumnOperator", false );
}

void ColumnSolver::apply( std::shared_ptr<const AMP::LinearAlgebra::Vector> f,
                          std::shared_ptr<AMP::LinearAlgebra::Vector> u )
{
    // u->zero();

    if ( d_IterationType == "GaussSeidel" ) {
        GaussSeidel( f, u );
    } else if ( d_IterationType == "SymmetricGaussSeidel" ) {
        SymmetricGaussSeidel( f, u );
    } else {
        AMP::pout << "ERROR: Invalid iteration type specified " << std::endl;
    }
}

void ColumnSolver::GaussSeidel( std::shared_ptr<const AMP::LinearAlgebra::Vector> &f,
                                std::shared_ptr<AMP::LinearAlgebra::Vector> &u )
{
    for ( int it = 0; it < d_iMaxIterations; it++ ) {
        for ( auto &elem : d_Solvers ) {
            std::shared_ptr<AMP::Operator::Operator> op = elem->getOperator();
            AMP_INSIST( op, "EROR: NULL Operator returned by SolverStrategy::getOperator" );

            std::shared_ptr<const AMP::LinearAlgebra::Vector> sf = op->subsetOutputVector( f );
            AMP_INSIST( sf, "ERROR: subset on rhs f yields NULL vector in ColumnSolver::solve" );
            std::shared_ptr<AMP::LinearAlgebra::Vector> su = op->subsetInputVector( u );
            AMP_INSIST( su,
                        "ERROR: subset on solution u yields NULL vector in ColumnSolver::solve" );

            elem->apply( sf, su );
        }
    }
}

void ColumnSolver::SymmetricGaussSeidel( std::shared_ptr<const AMP::LinearAlgebra::Vector> &f,
                                         std::shared_ptr<AMP::LinearAlgebra::Vector> &u )
{
    for ( int it = 0; it < d_iMaxIterations; it++ ) {
        for ( auto &elem : d_Solvers ) {
            std::shared_ptr<AMP::Operator::Operator> op = elem->getOperator();
            AMP_INSIST( op, "EROR: NULL Operator returned by SolverStrategy::getOperator" );

            std::shared_ptr<const AMP::LinearAlgebra::Vector> sf = op->subsetOutputVector( f );
            AMP_INSIST( sf, "ERROR: subset on rhs f yields NULL vector in ColumnSolver::solve" );
            std::shared_ptr<AMP::LinearAlgebra::Vector> su = op->subsetInputVector( u );
            AMP_INSIST( su,
                        "ERROR: subset on solution u yields NULL vector in ColumnSolver::solve" );

            elem->apply( sf, su );
        }

        for ( int i = (int) d_Solvers.size() - 1; i >= 0; i-- ) {
            std::shared_ptr<AMP::Operator::Operator> op = d_Solvers[i]->getOperator();
            AMP_INSIST( op, "EROR: NULL Operator returned by SolverStrategy::getOperator" );

            std::shared_ptr<const AMP::LinearAlgebra::Vector> sf = op->subsetOutputVector( f );
            AMP_INSIST( sf, "ERROR: subset on rhs f yields NULL vector in ColumnSolver::solve" );
            std::shared_ptr<AMP::LinearAlgebra::Vector> su = op->subsetInputVector( u );
            AMP_INSIST( su,
                        "ERROR: subset on solution u yields NULL vector in ColumnSolver::solve" );

            d_Solvers[i]->apply( sf, su );
        }
    }
}

void ColumnSolver::setInitialGuess( std::shared_ptr<AMP::LinearAlgebra::Vector> initialGuess )
{
    for ( auto &elem : d_Solvers ) {
        elem->setInitialGuess( initialGuess );
    }
}

void ColumnSolver::append( std::shared_ptr<AMP::Solver::SolverStrategy> solver )
{
    AMP_INSIST( ( solver ), "AMP::Solver::ColumnSolver::append input argument is a NULL solver" );
    d_Solvers.push_back( solver );
}

void ColumnSolver::resetOperator( std::shared_ptr<const AMP::Operator::OperatorParameters> params )
{
    if ( d_resetColumnOperator ) {
        d_pOperator->reset( params );

        std::shared_ptr<SolverStrategyParameters> solverParams;

        for ( auto &elem : d_Solvers ) {
            elem->reset( solverParams );
        }
    } else {
        auto columnParams =
            std::dynamic_pointer_cast<const AMP::Operator::ColumnOperatorParameters>( params );
        AMP_INSIST( columnParams, "Dynamic cast failed!" );

        for ( unsigned int i = 0; i < d_Solvers.size(); i++ ) {
            d_Solvers[i]->resetOperator( columnParams->d_OperatorParameters[i] );
        }
    }
}
} // namespace Solver
} // namespace AMP
