// This tests checks the creation of a TrilinosNOXSolver
#include "AMP/IO/PIO.h"
#include "AMP/mesh/Mesh.h"
#include "AMP/mesh/MeshFactory.h"
#include "AMP/operators/IdentityOperator.h"
#include "AMP/operators/NullOperator.h"
#include "AMP/solvers/trilinos/nox/TrilinosNOXSolver.h"
#include "AMP/utils/AMPManager.h"
#include "AMP/utils/AMP_MPI.h"
#include "AMP/utils/Database.h"
#include "AMP/utils/UnitTest.h"
#include "AMP/vectors/MultiVector.h"
#include "AMP/vectors/VectorBuilder.h"

#include <iostream>
#include <string>


static void myTest( AMP::UnitTest *ut, const std::string &exeName )
{
    std::string input_file = "input_" + exeName;
    std::string log_file   = "output_" + exeName;

    AMP::logOnlyNodeZero( log_file );
    AMP::AMP_MPI globalComm( AMP_COMM_WORLD );
    auto solverComm = globalComm.dup(); // Create a unique solver comm to test proper cleanup


    auto input_db = AMP::Database::parseInputFile( input_file );
    input_db->print( AMP::plog );

    // Create the solution and function variables
    auto var   = std::make_shared<AMP::LinearAlgebra::Variable>( "x" );
    auto u     = AMP::LinearAlgebra::createSimpleVector<double>( 25, var, solverComm );
    auto f     = u->clone();
    auto icVec = u->clone();

    // Create the operator
    auto op = std::make_shared<AMP::Operator::IdentityOperator>();
    op->setInputVariable( var );
    op->setOutputVariable( var );

    // Get the databases for the nonlinear and linear solvers
    auto nonlinearSolver_db = input_db->getDatabase( "NonlinearSolver" );
    // auto linearSolver_db = nonlinearSolver_db->getDatabase("LinearSolver");

    // initialize the nonlinear solver parameters
    auto nonlinearSolverParams =
        std::make_shared<AMP::Solver::TrilinosNOXSolverParameters>( nonlinearSolver_db );
    nonlinearSolverParams->d_comm            = solverComm;
    nonlinearSolverParams->d_pInitialGuess   = icVec;
    nonlinearSolverParams->d_pOperator       = op;
    nonlinearSolverParams->d_pLinearOperator = op;

    // Create the nonlinear solver
    auto nonlinearSolver =
        std::make_shared<AMP::Solver::TrilinosNOXSolver>( nonlinearSolverParams );
    ut->passes( "TrilinosNOXSolver created" );

    // Call solve with a simple vector
    u->setRandomValues();
    f->setRandomValues();
    nonlinearSolver->apply( f, u );
    ut->passes( "TrilinosNOXSolver solve called with simple vector" );
    auto x = u->clone();
    x->subtract( *u, *f );
    auto x_norm  = static_cast<double>( x->L2Norm() );
    auto f_norm  = static_cast<double>( f->L2Norm() );
    double error = x_norm / std::max( f_norm, 1.0 );
    if ( fabs( error ) < 1e-8 )
        ut->passes( "Solve with simple vector passed" );
    else
        ut->failure( "Solve with simple vector failed" );


    // Call solve with a multivector
    // (there can be bugs when solve is called with a single vector and then a multivector)
    auto mu = AMP::LinearAlgebra::MultiVector::create( "multivector", solverComm );
    auto mf = AMP::LinearAlgebra::MultiVector::create( "multivector", solverComm );
    mu->addVector( u );
    mf->addVector( f );
    mu->setRandomValues();
    mf->zero();
    nonlinearSolver->apply( mf, mu );
    ut->passes( "TrilinosNOXSolver solve called with multivector" );
}


int testTrilinosNOXSolver( int argc, char *argv[] )
{
    AMP::AMPManager::startup( argc, argv );
    AMP::UnitTest ut;

    myTest( &ut, "testTrilinosNOXSolver" );

    ut.report();

    int num_failed = ut.NumFailGlobal();
    AMP::AMPManager::shutdown();
    return num_failed;
}
