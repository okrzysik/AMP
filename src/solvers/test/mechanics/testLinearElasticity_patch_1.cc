#include "AMP/ampmesh/Mesh.h"
#include "AMP/ampmesh/libmesh/libmeshMesh.h"
#include "AMP/discretization/simpleDOF_Manager.h"
#include "AMP/operators/LinearBVPOperator.h"
#include "AMP/operators/OperatorBuilder.h"
#include "AMP/operators/boundary/DirichletVectorCorrection.h"
#include "AMP/operators/mechanics/MechanicsLinearFEOperator.h"
#include "AMP/solvers/trilinos/ml/TrilinosMLSolver.h"
#include "AMP/utils/AMPManager.h"
#include "AMP/utils/AMP_MPI.h"
#include "AMP/utils/Database.h"
#include "AMP/utils/PIO.h"
#include "AMP/utils/ReadTestMesh.h"
#include "AMP/utils/UnitTest.h"
#include "AMP/utils/Utilities.h"
#include "AMP/utils/Writer.h"
#include "AMP/vectors/Vector.h"
#include "AMP/vectors/VectorBuilder.h"

#include <fstream>
#include <iostream>
#include <string>


static void linearElasticTest( AMP::UnitTest *ut, const std::string &exeName )
{
    std::string input_file = "input_" + exeName;
    std::string log_file   = "output_" + exeName;

    AMP::PIO::logOnlyNodeZero( log_file );
    AMP::AMP_MPI globalComm = AMP::AMP_MPI( AMP_COMM_WORLD );

    if ( globalComm.getSize() == 1 ) {

        auto input_db = AMP::Database::parseInputFile( input_file );
        input_db->print( AMP::plog );

        auto mesh_file_db = AMP::Database::parseInputFile( input_db->getString( "mesh_file" ) );

        const unsigned int mesh_dim = 3;
        auto mesh                   = std::make_shared<libMesh::Mesh>(libMesh::Parallel::Communicator(), mesh_dim );

        AMP::readTestMesh( mesh_file_db, mesh );
        mesh->prepare_for_use( false );

        auto meshAdapter = std::make_shared<AMP::Mesh::libmeshMesh>( mesh, "TestMesh" );

        std::shared_ptr<AMP::Operator::ElementPhysicsModel> elementPhysicsModel;
        auto bvpOperator = std::dynamic_pointer_cast<AMP::Operator::LinearBVPOperator>(
            AMP::Operator::OperatorBuilder::createOperator(
                meshAdapter, "MechanicsBVPOperator", input_db, elementPhysicsModel ) );

        std::shared_ptr<AMP::Operator::ElementPhysicsModel> dummyModel;
        auto dirichletVecOp = std::dynamic_pointer_cast<AMP::Operator::DirichletVectorCorrection>(
            AMP::Operator::OperatorBuilder::createOperator(
                meshAdapter, "Load_Boundary", input_db, dummyModel ) );

        auto var = bvpOperator->getOutputVariable();

        // This has an in-place apply. So, it has an empty input variable and
        // the output variable is the same as what it is operating on.
        dirichletVecOp->setVariable( var );

        auto dofMap = AMP::Discretization::simpleDOFManager::create(
            meshAdapter, AMP::Mesh::GeomType::Vertex, 1, 3, true );

        AMP::LinearAlgebra::Vector::shared_ptr nullVec;
        auto mechSolVec = AMP::LinearAlgebra::createVector( dofMap, var, true );
        auto mechRhsVec = mechSolVec->cloneVector();
        auto mechResVec = mechSolVec->cloneVector();

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

        auto mlSolver_db = input_db->getDatabase( "LinearSolver" );

        auto mlSolverParams =
            std::make_shared<AMP::Solver::SolverStrategyParameters>( mlSolver_db );

        mlSolverParams->d_pOperator = bvpOperator;

        // create the ML solver interface
        auto mlSolver = std::make_shared<AMP::Solver::TrilinosMLSolver>( mlSolverParams );

        mlSolver->setZeroInitialGuess( false );

        mlSolver->solve( mechRhsVec, mechSolVec );

        double finalSolNorm = mechSolVec->L2Norm();

        AMP::pout << "Final Solution Norm: " << finalSolNorm << std::endl;

        std::string fname = exeName + "-stressAndStrain.txt";

        ( std::dynamic_pointer_cast<AMP::Operator::MechanicsLinearFEOperator>(
              bvpOperator->getVolumeOperator() ) )
            ->printStressAndStrain( mechSolVec, fname );

        bvpOperator->residual( mechRhsVec, mechSolVec, mechResVec );

        double finalResidualNorm = mechResVec->L2Norm();

        AMP::pout << "Final Residual Norm: " << finalResidualNorm << std::endl;

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

int testLinearElasticity_patch_1( int argc, char *argv[] )
{
    AMP::AMPManager::startup( argc, argv );
    auto libmeshInit = std::make_shared<AMP::Mesh::initializeLibMesh>( AMP_COMM_WORLD );

    AMP::UnitTest ut;

    std::vector<std::string> exeNames;
    exeNames.emplace_back( "testLinearElasticity-patch-1-normal" );
    exeNames.emplace_back( "testLinearElasticity-patch-1-reduced" );

    for ( auto &exeName : exeNames )
        linearElasticTest( &ut, exeName );

    ut.report();
    int num_failed = ut.NumFailGlobal();

    libmeshInit.reset();
    AMP::AMPManager::shutdown();
    return num_failed;
}
