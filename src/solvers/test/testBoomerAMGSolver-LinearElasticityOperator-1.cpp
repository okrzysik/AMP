#include "AMP/IO/PIO.h"
#include "AMP/IO/Writer.h"
#include "AMP/discretization/simpleDOF_Manager.h"
#include "AMP/mesh/MeshParameters.h"
#include "AMP/operators/LinearBVPOperator.h"
#include "AMP/operators/OperatorBuilder.h"
#include "AMP/operators/boundary/DirichletVectorCorrection.h"
#include "AMP/solvers/hypre/BoomerAMGSolver.h"
#include "AMP/utils/AMPManager.h"
#include "AMP/utils/AMP_MPI.h"
#include "AMP/utils/Database.h"
#include "AMP/utils/UnitTest.h"
#include "AMP/utils/Utilities.h"
#include "AMP/vectors/VectorBuilder.h"


void linearElasticTest( AMP::UnitTest *ut )
{
    std::string exeName( "testBoomerAMGSolver-LinearElasticityOperator-1" );
    std::string input_file = "input_" + exeName;
    std::string log_file   = "output_" + exeName;

    AMP::logOnlyNodeZero( log_file );

    auto input_db = AMP::Database::parseInputFile( input_file );
    input_db->print( AMP::plog );

    AMP_INSIST( input_db->keyExists( "Mesh" ), "Key ''Mesh'' is missing!" );
    auto mesh_db    = input_db->getDatabase( "Mesh" );
    auto meshParams = std::make_shared<AMP::Mesh::MeshParameters>( mesh_db );
    meshParams->setComm( AMP::AMP_MPI( AMP_COMM_WORLD ) );
    auto meshAdapter = AMP::Mesh::Mesh::buildMesh( meshParams );

    std::shared_ptr<AMP::Operator::ElementPhysicsModel> elementPhysicsModel;
    auto bvpOperator = std::dynamic_pointer_cast<AMP::Operator::LinearBVPOperator>(
        AMP::Operator::OperatorBuilder::createOperator(
            meshAdapter, "MechanicsBVPOperator", input_db, elementPhysicsModel ) );

    auto var = bvpOperator->getOutputVariable();

    std::shared_ptr<AMP::Operator::ElementPhysicsModel> dummyModel;
    auto dirichletVecOp = std::dynamic_pointer_cast<AMP::Operator::DirichletVectorCorrection>(
        AMP::Operator::OperatorBuilder::createOperator(
            meshAdapter, "Load_Boundary", input_db, dummyModel ) );
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


    std::cout << "RHS Norm: " << mechRhsVec->L2Norm() << std::endl;
    std::cout << "Initial Solution Norm: " << mechSolVec->L2Norm() << std::endl;

    bvpOperator->residual( mechRhsVec, mechSolVec, mechResVec );

    double initResidualNorm = static_cast<double>( mechResVec->L2Norm() );
    std::cout << "Initial Residual Norm: " << initResidualNorm << std::endl;

    auto mlSolver_db = input_db->getDatabase( "LinearSolver" );

    auto mlSolverParams = std::make_shared<AMP::Solver::SolverStrategyParameters>( mlSolver_db );

    mlSolverParams->d_pOperator = bvpOperator;

    // create the ML solver interface
    auto mlSolver = std::make_shared<AMP::Solver::BoomerAMGSolver>( mlSolverParams );

    mlSolver->setZeroInitialGuess( false );

    mlSolver->apply( mechRhsVec, mechSolVec );

#ifdef USE_EXT_SILO
    // Create the silo writer and register the data
    auto siloWriter = AMP::IO::Writer::buildWriter( "Silo" );
    siloWriter->registerVector( mechSolVec, meshAdapter, AMP::Mesh::GeomType::Vertex, "Solution" );
    siloWriter->writeFile( exeName, 0 );
#endif

    bvpOperator->residual( mechRhsVec, mechSolVec, mechResVec );

    double finalResidualNorm = static_cast<double>( mechResVec->L2Norm() );

    std::cout << "Final Residual Norm: " << finalResidualNorm << std::endl;

    // BP: convergence assumes convergence rate of ~0.6 and 40 iterations
    if ( finalResidualNorm > ( 1.0e-8 * initResidualNorm ) ) {
        ut->failure( "BoomerAMGSolver fails to solve a linear elasticity problem" );
    } else {
        ut->passes( "BoomerAMGSolver successfully solves a linear elasticity problem" );
    }

    input_db.reset();
}

int main( int argc, char *argv[] )
{
    AMP::AMPManager::startup( argc, argv );
    AMP::UnitTest ut;

    linearElasticTest( &ut );

    ut.report();

    int num_failed = ut.NumFailGlobal();
    AMP::AMPManager::shutdown();
    return num_failed;
}
