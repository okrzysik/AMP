
#include <cstdio>
#include <cstring>
#include <iostream>
#include <vector>

#include "AMP/utils/AMPManager.h"
#include "AMP/utils/Database.h"
#include "AMP/utils/PIO.h"
#include "AMP/utils/UnitTest.h"

#include "AMP/vectors/MultiVector.h"
#include "AMP/vectors/SimpleVector.h"
#include "AMP/vectors/newFrozenVectorDesign/newFrozenVectorDesignHelpers.h"

#include "AMP/operators/ColumnOperator.h"
#include "AMP/operators/newFrozenVectorDesign/FirstOperator.h"
#include "AMP/operators/newFrozenVectorDesign/SecondOperator.h"

#include "AMP/solvers/ColumnSolver.h"
#include "AMP/solvers/OnePointSolver.h"

int main( int argc, char *argv[] )
{
    AMP::AMPManager::startup( argc, argv );
    AMP::UnitTest ut;

    std::string exeName    = "testNewFrozenVectorDesign";
    std::string input_file = "input_" + exeName;
    std::string log_file   = "output_" + exeName;
    AMP::PIO::logOnlyNodeZero( log_file );


    auto input_db = AMP::Database::parseInputFile( input_file );

    std::shared_ptr<AMP::Database> firstOp_db = input_db->getDatabase( "FirstOperator" );
    std::shared_ptr<AMP::Operator::OperatorParameters> firstOpParams(
        new AMP::Operator::OperatorParameters( firstOp_db ) );
    std::shared_ptr<AMP::Operator::FirstOperator> firstOp(
        new AMP::Operator::FirstOperator( firstOpParams ) );

    std::shared_ptr<AMP::Database> secondOp_db = input_db->getDatabase( "SecondOperator" );
    std::shared_ptr<AMP::Operator::OperatorParameters> secondOpParams(
        new AMP::Operator::OperatorParameters( secondOp_db ) );
    std::shared_ptr<AMP::Operator::SecondOperator> secondOp(
        new AMP::Operator::SecondOperator( secondOpParams ) );

    std::shared_ptr<AMP::Operator::ColumnOperator> columnOp( new AMP::Operator::ColumnOperator );
    columnOp->append( firstOp );
    columnOp->append( secondOp );

    std::cout << "Constructed Operators." << std::endl;

    std::shared_ptr<AMP::Database> solver_db = input_db->getDatabase( "Solver" );

    std::shared_ptr<AMP::Solver::SolverStrategyParameters> firstSolverParams(
        new AMP::Solver::SolverStrategyParameters( solver_db ) );
    firstSolverParams->d_pOperator = firstOp;
    std::shared_ptr<AMP::Solver::OnePointSolver> firstSolver(
        new AMP::Solver::OnePointSolver( firstSolverParams ) );

    std::shared_ptr<AMP::Solver::SolverStrategyParameters> secondSolverParams(
        new AMP::Solver::SolverStrategyParameters( solver_db ) );
    secondSolverParams->d_pOperator = secondOp;
    std::shared_ptr<AMP::Solver::OnePointSolver> secondSolver(
        new AMP::Solver::OnePointSolver( secondSolverParams ) );

    std::shared_ptr<AMP::Solver::SolverStrategyParameters> columnSolverParams(
        new AMP::Solver::SolverStrategyParameters( solver_db ) );
    columnSolverParams->d_pOperator = columnOp;
    std::shared_ptr<AMP::Solver::ColumnSolver> columnSolver(
        new AMP::Solver::ColumnSolver( columnSolverParams ) );
    columnSolver->append( firstSolver );
    columnSolver->append( secondSolver );

    std::cout << "Constructed Solvers." << std::endl;

    AMP::LinearAlgebra::Variable::shared_ptr firstVar  = firstOp->getOutputVariable();
    AMP::LinearAlgebra::Variable::shared_ptr secondVar = secondOp->getOutputVariable();

    AMP::LinearAlgebra::Vector::shared_ptr firstVec =
        AMP::LinearAlgebra::SimpleVector<double>::create( 1, firstVar );
    AMP::LinearAlgebra::Vector::shared_ptr secondVec =
        AMP::LinearAlgebra::SimpleVector<double>::create( 1, secondVar );

    AMP::LinearAlgebra::CommunicationList::shared_ptr commList =
        AMP::LinearAlgebra::CommunicationList::createEmpty( 1, AMP_COMM_SELF );

    firstVec->setCommunicationList( commList );
    secondVec->setCommunicationList( commList );

    AMP::LinearAlgebra::Vector::shared_ptr rhsVec =
        AMP::LinearAlgebra::joinVectors( firstVec, secondVec );
    AMP::LinearAlgebra::Vector::shared_ptr solVec = rhsVec->cloneVector();

    std::cout << "Constructed Vectors." << std::endl;

    firstVec->setToScalar( 3.0 );
    secondVec->setToScalar( 5.0 );

    std::cout << "Formed RHS." << std::endl;

    columnSolver->solve( rhsVec, solVec );

    std::cout << "Completed Solve." << std::endl;

    auto firstSol = std::dynamic_pointer_cast<AMP::LinearAlgebra::SimpleVector<double>>(
        solVec->subsetVectorForVariable( firstVar ) );

    auto secondSol = std::dynamic_pointer_cast<AMP::LinearAlgebra::SimpleVector<double>>(
        solVec->subsetVectorForVariable( secondVar ) );

    std::cout << "First Solution = " << ( ( *firstSol )[0] ) << std::endl;
    std::cout << "Second Solution = " << ( ( *secondSol )[0] ) << std::endl;

    ut.passes( exeName );

    ut.report();
    int num_failed = ut.NumFailGlobal();

    AMP::AMPManager::shutdown();

    return num_failed;
}
