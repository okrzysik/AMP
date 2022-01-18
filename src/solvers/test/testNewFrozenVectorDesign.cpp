#include "AMP/IO/PIO.h"
#include "AMP/operators/ColumnOperator.h"
#include "AMP/operators/newFrozenVectorDesign/FirstOperator.h"
#include "AMP/operators/newFrozenVectorDesign/SecondOperator.h"
#include "AMP/solvers/ColumnSolver.h"
#include "AMP/solvers/OnePointSolver.h"
#include "AMP/utils/AMPManager.h"
#include "AMP/utils/Database.h"
#include "AMP/utils/UnitTest.h"
#include "AMP/vectors/MultiVector.h"
#include "AMP/vectors/VectorBuilder.h"
#include "AMP/vectors/newFrozenVectorDesign/newFrozenVectorDesignHelpers.h"

#include <cstdio>
#include <cstring>
#include <iostream>
#include <vector>


void run( AMP::UnitTest &ut )
{
    std::string exeName    = "testNewFrozenVectorDesign";
    std::string input_file = "input_" + exeName;
    std::string log_file   = "output_" + exeName;
    AMP::logOnlyNodeZero( log_file );


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

    auto solver_db = input_db->getDatabase( "Solver" );

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

    auto firstVec  = AMP::LinearAlgebra::createSimpleVector<double>( 1, firstVar );
    auto secondVec = AMP::LinearAlgebra::createSimpleVector<double>( 1, secondVar );

    auto commList = AMP::LinearAlgebra::CommunicationList::createEmpty( 1, AMP_COMM_SELF );

    firstVec->setCommunicationList( commList );
    secondVec->setCommunicationList( commList );

    auto rhsVec = AMP::LinearAlgebra::joinVectors( firstVec, secondVec );
    auto solVec = rhsVec->cloneVector();

    std::cout << "Constructed Vectors." << std::endl;

    firstVec->setToScalar( 3.0 );
    secondVec->setToScalar( 5.0 );

    std::cout << "Formed RHS." << std::endl;

    columnSolver->apply( rhsVec, solVec );

    std::cout << "Completed Solve." << std::endl;

    auto firstSol = solVec->subsetVectorForVariable( firstVar );

    auto secondSol = solVec->subsetVectorForVariable( secondVar );

    auto firstSolData  = firstSol->getVectorData();
    auto secondSolData = secondSol->getVectorData();
    // std::cout << "First Solution = " << ( ( *firstSolData )[0] ) << std::endl;
    // std::cout << "Second Solution = " << ( ( *secondSolData )[0] ) << std::endl;

    ut.passes( exeName );
}


int main( int argc, char *argv[] )
{
    AMP::AMPManager::startup( argc, argv );
    AMP::UnitTest ut;

    run( ut );

    int N_failed = ut.NumFailGlobal();
    ut.report();
    ut.reset();
    AMP::AMPManager::shutdown();
    return N_failed;
}
