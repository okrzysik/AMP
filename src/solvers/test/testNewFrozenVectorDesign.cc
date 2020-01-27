#include "AMP/operators/ColumnOperator.h"
#include "AMP/operators/newFrozenVectorDesign/FirstOperator.h"
#include "AMP/operators/newFrozenVectorDesign/SecondOperator.h"
#include "AMP/solvers/ColumnSolver.h"
#include "AMP/solvers/OnePointSolver.h"
#include "AMP/utils/AMPManager.h"
#include "AMP/utils/Database.h"
#include "AMP/utils/PIO.h"
#include "AMP/utils/UnitTest.h"
#include "AMP/vectors/MultiVector.h"
#include "AMP/vectors/SimpleVector.h"
#include "AMP/vectors/newFrozenVectorDesign/newFrozenVectorDesignHelpers.h"

#include <cstdio>
#include <cstring>
#include <iostream>
#include <vector>


int main( int argc, char *argv[] )
{
    AMP::AMPManager::startup( argc, argv );
    AMP::UnitTest ut;

    std::string exeName    = "testNewFrozenVectorDesign";
    std::string input_file = "input_" + exeName;
    std::string log_file   = "output_" + exeName;
    AMP::PIO::logOnlyNodeZero( log_file );


    auto input_db = AMP::Database::parseInputFile( input_file );

    auto firstOp_db    = input_db->getDatabase( "FirstOperator" );
    auto firstOpParams = std::make_shared<AMP::Operator::OperatorParameters>( firstOp_db );
    auto firstOp       = std::make_shared<AMP::Operator::FirstOperator>( firstOpParams );

    auto secondOp_db    = input_db->getDatabase( "SecondOperator" );
    auto secondOpParams = std::make_shared<AMP::Operator::OperatorParameters>( secondOp_db );
    auto secondOp       = std::make_shared<AMP::Operator::SecondOperator>( secondOpParams );

    auto columnOp = std::make_shared<AMP::Operator::ColumnOperator>();
    columnOp->append( firstOp );
    columnOp->append( secondOp );

    std::cout << "Constructed Operators." << std::endl;

    std::shared_ptr<AMP::Database> solver_db = input_db->getDatabase( "Solver" );

    auto firstSolverParams = std::make_shared<AMP::Solver::SolverStrategyParameters>( solver_db );
    firstSolverParams->d_pOperator = firstOp;
    auto firstSolver = std::make_shared<AMP::Solver::OnePointSolver>( firstSolverParams );

    auto secondSolverParams = std::make_shared<AMP::Solver::SolverStrategyParameters>( solver_db );
    secondSolverParams->d_pOperator = secondOp;
    auto secondSolver = std::make_shared<AMP::Solver::OnePointSolver>( secondSolverParams );

    auto columnSolverParams = std::make_shared<AMP::Solver::SolverStrategyParameters>( solver_db );
    columnSolverParams->d_pOperator = columnOp;
    auto columnSolver = std::make_shared<AMP::Solver::ColumnSolver>( columnSolverParams );
    columnSolver->append( firstSolver );
    columnSolver->append( secondSolver );

    std::cout << "Constructed Solvers." << std::endl;

    auto firstVar  = firstOp->getOutputVariable();
    auto secondVar = secondOp->getOutputVariable();

    auto firstVec  = AMP::LinearAlgebra::SimpleVector<double>::create( 1, firstVar );
    auto secondVec = AMP::LinearAlgebra::SimpleVector<double>::create( 1, secondVar );

    auto commList = AMP::LinearAlgebra::CommunicationList::createEmpty( 1, AMP_COMM_SELF );

    firstVec->setCommunicationList( commList );
    secondVec->setCommunicationList( commList );

    auto rhsVec = AMP::LinearAlgebra::joinVectors( firstVec, secondVec );
    auto solVec = rhsVec->cloneVector();

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
