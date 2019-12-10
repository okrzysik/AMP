// This tests checks the creation of a PetscSNESSolver
// Note: the comm used should NOT be comm_world as there are cleanup issues for other comms when
// using the monitor
// option
#include "AMP/utils/AMPManager.h"
#include "AMP/utils/AMP_MPI.h"
#include "AMP/utils/Database.h"
#include "AMP/utils/PIO.h"
#include "AMP/utils/UnitTest.h"
#include "AMP/utils/Utilities.h"
#include <iostream>
#include <string>

#include "AMP/ampmesh/Mesh.h"
#include "AMP/operators/IdentityOperator.h"
#include "AMP/operators/NullOperator.h"
#include "AMP/solvers/petsc/PetscSNESSolver.h"
#include "AMP/vectors/MultiVector.h"
#include "AMP/vectors/NullVector.h"
#include "AMP/vectors/SimpleVector.h"


static void myTest( AMP::UnitTest *ut, const std::string &exeName )
{
    std::string input_file = "input_" + exeName;
    std::string log_file   = "output_" + exeName;

    AMP::PIO::logOnlyNodeZero( log_file );
    AMP::AMP_MPI globalComm( AMP_COMM_WORLD );
    auto solverComm = globalComm.dup(); // Create a unique solver comm to test proper cleanup

    auto input_db = AMP::Database::parseInputFile( input_file );
    input_db->print( AMP::plog );

    // Create a null vector for the initial guess
    auto nullVec = AMP::LinearAlgebra::NullVector<double>::create( "null" );

    // Create the solution and function variables
    auto var = AMP::make_shared<AMP::LinearAlgebra::Variable>( "x" );
    auto u   = AMP::LinearAlgebra::SimpleVector<double>::create( 10, var, solverComm );
    auto f   = u->cloneVector();

    // Create the operator
    auto op = AMP::make_shared<AMP::Operator::IdentityOperator>();
    op->setInputVariable( var );
    op->setOutputVariable( var );

    // Get the databases for the nonlinear and linear solvers
    auto nonlinearSolver_db = input_db->getDatabase( "NonlinearSolver" );
    // auto linearSolver_db = nonlinearSolver_db->getDatabase("LinearSolver");

    // initialize the nonlinear solver parameters
    auto nonlinearSolverParams =
        AMP::make_shared<AMP::Solver::PetscSNESSolverParameters>( nonlinearSolver_db );
    nonlinearSolverParams->d_comm          = solverComm;
    nonlinearSolverParams->d_pInitialGuess = nullVec;
    nonlinearSolverParams->d_pOperator     = op;


    // Create the nonlinear solver
    auto nonlinearSolver = AMP::make_shared<AMP::Solver::PetscSNESSolver>( nonlinearSolverParams );
    ut->passes( "PetscSNESSolver created" );

    // Call solve with a simple vector
    u->setRandomValues();
    f->setRandomValues();
    nonlinearSolver->solve( f, u );
    ut->passes( "PetscSNESSolver solve called with simple vector" );
    auto x = u->cloneVector();
    x->subtract( u, f );
    double error = x->L2Norm() / f->L2Norm();
    if ( fabs( error ) < 1e-8 )
        ut->passes( "Solve with simple vector passed" );
    else
        ut->failure( "Solve with simple vector failed" );

    // Call solve with a multivector (there can be bugs when solve is called with a single vector
    // and then a
    // multivector)
    auto mu = AMP::LinearAlgebra::MultiVector::create( "multivector", solverComm );
    auto mf = AMP::LinearAlgebra::MultiVector::create( "multivector", solverComm );
    mu->addVector( u );
    mf->addVector( f );
    mu->setRandomValues();
    mf->zero();
    nonlinearSolver->solve( mf, mu );
    ut->passes( "PetscSNESSolver solve called with multivector" );
}


int testPetscSNESSolver( int argc, char *argv[] )
{
    AMP::AMPManager::startup( argc, argv );
    AMP::UnitTest ut;

    myTest( &ut, "testPetscSNESSolver" );

    ut.report();

    int num_failed = ut.NumFailGlobal();
    AMP::AMPManager::shutdown();
    return num_failed;
}
