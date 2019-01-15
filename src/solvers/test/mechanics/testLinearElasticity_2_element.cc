

#include <fstream>
#include <iostream>
#include <string>

/* AMP files */
#include "AMP/utils/AMPManager.h"
#include "AMP/utils/AMP_MPI.h"
#include "AMP/utils/Database.h"
#include "AMP/utils/InputDatabase.h"
#include "AMP/utils/InputManager.h"
#include "AMP/utils/PIO.h"
#include "AMP/utils/ReadTestMesh.h"
#include "AMP/utils/UnitTest.h"
#include "AMP/utils/Utilities.h"
#include "AMP/utils/WriteSolutionToFile.h"

#include "AMP/ampmesh/Mesh.h"
#include "AMP/ampmesh/libmesh/libMesh.h"

#include "AMP/discretization/simpleDOF_Manager.h"
#include "AMP/vectors/VectorBuilder.h"

#include "AMP/operators/LinearBVPOperator.h"
#include "AMP/operators/OperatorBuilder.h"
#include "AMP/operators/boundary/DirichletVectorCorrection.h"

#include "AMP/utils/Writer.h"
#include "AMP/vectors/Vector.h"

#include "AMP/solvers/trilinos/ml/TrilinosMLSolver.h"

static void linearElasticTest( AMP::UnitTest *ut, int reduced, std::string mesh_file )
{
    std::string exeName;
    std::string input_file;

    if ( reduced ) {
        exeName    = "testLinearElasticity-reduced-" + mesh_file;
        input_file = "input_testLinearElasticity-reduced-mesh2elem";
    } else {
        exeName    = "testLinearElasticity-normal-" + mesh_file;
        input_file = "input_testLinearElasticity-normal-mesh2elem";
    }

    std::string log_file = "output_" + exeName;

    AMP::PIO::logOnlyNodeZero( log_file );
    AMP::AMP_MPI globalComm = AMP::AMP_MPI( AMP_COMM_WORLD );

    if ( globalComm.getSize() == 1 ) {
        AMP::shared_ptr<AMP::InputDatabase> input_db( new AMP::InputDatabase( "input_db" ) );
        AMP::InputManager::getManager()->parseInputFile( input_file, input_db );
        input_db->printClassData( AMP::plog );

        AMP::shared_ptr<AMP::InputDatabase> mesh_file_db(
            new AMP::InputDatabase( "mesh_file_db" ) );
        AMP::InputManager::getManager()->parseInputFile( mesh_file, mesh_file_db );

        const unsigned int mesh_dim = 3;
        AMP::shared_ptr<::Mesh> mesh( new ::Mesh( mesh_dim ) );

        AMP::readTestMesh( mesh_file_db, mesh );
        mesh->prepare_for_use( false );

        AMP::Mesh::Mesh::shared_ptr meshAdapter =
            AMP::Mesh::Mesh::shared_ptr( new AMP::Mesh::libMesh( mesh, "TestMesh" ) );

        AMP::shared_ptr<AMP::Operator::ElementPhysicsModel> elementPhysicsModel;
        AMP::shared_ptr<AMP::Operator::LinearBVPOperator> bvpOperator =
            AMP::dynamic_pointer_cast<AMP::Operator::LinearBVPOperator>(
                AMP::Operator::OperatorBuilder::createOperator(
                    meshAdapter, "MechanicsBVPOperator", input_db, elementPhysicsModel ) );

        AMP::shared_ptr<AMP::Operator::ElementPhysicsModel> dummyModel;
        AMP::shared_ptr<AMP::Operator::DirichletVectorCorrection> dirichletVecOp =
            AMP::dynamic_pointer_cast<AMP::Operator::DirichletVectorCorrection>(
                AMP::Operator::OperatorBuilder::createOperator(
                    meshAdapter, "Load_Boundary", input_db, dummyModel ) );

        AMP::LinearAlgebra::Variable::shared_ptr var = bvpOperator->getOutputVariable();

        // This has an in-place apply. So, it has an empty input variable and
        // the output variable is the same as what it is operating on.
        dirichletVecOp->setVariable( var );

        AMP::Discretization::DOFManager::shared_ptr dofMap =
            AMP::Discretization::simpleDOFManager::create(
                meshAdapter, AMP::Mesh::GeomType::Vertex, 1, 3, true );

        AMP::LinearAlgebra::Vector::shared_ptr nullVec;
        AMP::LinearAlgebra::Vector::shared_ptr mechSolVec =
            AMP::LinearAlgebra::createVector( dofMap, var, true );
        AMP::LinearAlgebra::Vector::shared_ptr mechRhsVec = mechSolVec->cloneVector();
        AMP::LinearAlgebra::Vector::shared_ptr mechResVec = mechSolVec->cloneVector();

        mechSolVec->setToScalar( 0.5 );
        mechRhsVec->setToScalar( 0.0 );
        mechResVec->setToScalar( 0.0 );

        dirichletVecOp->apply( nullVec, mechRhsVec );

        double rhsNorm = mechRhsVec->L2Norm();

        AMP::pout << "RHS Norm: " << rhsNorm << std::endl;

        double initSolNorm = mechSolVec->L2Norm();

        AMP::pout << "Initial Solution Norm: " << initSolNorm << std::endl;

        bvpOperator->residual( mechRhsVec, mechSolVec, mechResVec );

        double initResidualNorm = mechResVec->L2Norm();

        AMP::pout << "Initial Residual Norm: " << initResidualNorm << std::endl;

        AMP::shared_ptr<AMP::Database> mlSolver_db = input_db->getDatabase( "LinearSolver" );

        AMP::shared_ptr<AMP::Solver::SolverStrategyParameters> mlSolverParams(
            new AMP::Solver::SolverStrategyParameters( mlSolver_db ) );

        mlSolverParams->d_pOperator = bvpOperator;

        // create the ML solver interface
        AMP::shared_ptr<AMP::Solver::TrilinosMLSolver> mlSolver(
            new AMP::Solver::TrilinosMLSolver( mlSolverParams ) );

        mlSolver->setZeroInitialGuess( false );

        mlSolver->solve( mechRhsVec, mechSolVec );

        double finalSolNorm = mechSolVec->L2Norm();

        AMP::pout << "Final Solution Norm: " << finalSolNorm << std::endl;

        bvpOperator->residual( mechRhsVec, mechSolVec, mechResVec );

        double finalResidualNorm = mechResVec->L2Norm();

        AMP::pout << "Final Residual Norm: " << finalResidualNorm << std::endl;

        printSolution( meshAdapter, mechSolVec, exeName );

        if ( finalResidualNorm > ( 1e-10 * initResidualNorm ) ) {
            ut->failure( exeName );
        } else {
            ut->passes( exeName );
        }
    } else {
        AMP::pout << "WARNING: This is a single processor test!" << std::endl;
        ut->passes( exeName );
    }
}

int testLinearElasticity_2_element( int argc, char *argv[] )
{
    AMP::AMPManager::startup( argc, argv );
    AMP::shared_ptr<AMP::Mesh::initializeLibMesh> libmeshInit(
        new AMP::Mesh::initializeLibMesh( AMP_COMM_WORLD ) );

    AMP::UnitTest ut;

    std::vector<std::string> mesh_files;

    if ( argc == 1 ) {
        for ( int i = 1; i <= 6; i++ ) {
            char name[100];
            sprintf( name, "mesh2elem-%d", i );
            mesh_files.emplace_back( name );
        } // end for i
    } else {
        for ( int i = 1; i < argc; i++ ) {
            char name[100];
            sprintf( name, "mesh2elem-%d", atoi( argv[i] ) );
            mesh_files.emplace_back( name );
        } // end for i
    }

    for ( auto &mesh_file : mesh_files ) {
        for ( int reduced = 0; reduced < 2; reduced++ ) {
            try {
                linearElasticTest( &ut, reduced, mesh_file );
            } catch ( std::exception &err ) {
                AMP::pout << "ERROR: " << err.what() << std::endl;
                ut.failure( "ERROR" );
            } catch ( ... ) {
                AMP::pout << "ERROR: "
                          << "An unknown exception was thrown." << std::endl;
                ut.failure( "ERROR" );
            }
        } // end for reduced
    }     // end for i

    ut.report();
    int num_failed = ut.NumFailGlobal();

    libmeshInit.reset();
    AMP::AMPManager::shutdown();
    return num_failed;
}
