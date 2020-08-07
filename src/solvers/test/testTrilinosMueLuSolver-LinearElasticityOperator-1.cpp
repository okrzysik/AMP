#include "AMP/discretization/simpleDOF_Manager.h"
#include "AMP/operators/LinearBVPOperator.h"
#include "AMP/operators/OperatorBuilder.h"
#include "AMP/operators/boundary/DirichletVectorCorrection.h"
#include "AMP/solvers/trilinos/muelu/TrilinosMueLuSolver.h"
#include "AMP/utils/AMPManager.h"
#include "AMP/utils/AMP_MPI.h"
#include "AMP/utils/Database.h"
#include "AMP/utils/PIO.h"
#include "AMP/utils/UnitTest.h"
#include "AMP/utils/Utilities.h"
#include "AMP/utils/Writer.h"
#include "AMP/vectors/VectorBuilder.h"

#include <iostream>
#include <string>


void linearElasticTest( AMP::UnitTest *ut )
{
    std::string exeName( "testTrilinosMueLuSolver-LinearElasticityOperator-1" );
    std::string input_file = "input_" + exeName;
    std::string log_file   = "output_" + exeName;

    AMP::PIO::logOnlyNodeZero( log_file );

#ifdef USE_EXT_SILO
    // Create the silo writer and register the data
    AMP::Utilities::Writer::shared_ptr siloWriter = AMP::Utilities::Writer::buildWriter( "Silo" );
#endif

    auto input_db = AMP::Database::parseInputFile( input_file );
    input_db->print( AMP::plog );

    AMP_INSIST( input_db->keyExists( "Mesh" ), "Key ''Mesh'' is missing!" );
    std::shared_ptr<AMP::Database> mesh_db = input_db->getDatabase( "Mesh" );
    std::shared_ptr<AMP::Mesh::MeshParameters> meshParams(
        new AMP::Mesh::MeshParameters( mesh_db ) );
    meshParams->setComm( AMP::AMP_MPI( AMP_COMM_WORLD ) );
    AMP::Mesh::Mesh::shared_ptr meshAdapter = AMP::Mesh::Mesh::buildMesh( meshParams );

    std::shared_ptr<AMP::Operator::ElementPhysicsModel> elementPhysicsModel;
    std::shared_ptr<AMP::Operator::LinearBVPOperator> bvpOperator =
        std::dynamic_pointer_cast<AMP::Operator::LinearBVPOperator>(
            AMP::Operator::OperatorBuilder::createOperator(
                meshAdapter, "MechanicsBVPOperator", input_db, elementPhysicsModel ) );

    AMP::LinearAlgebra::Variable::shared_ptr var = bvpOperator->getOutputVariable();

    std::shared_ptr<AMP::Operator::ElementPhysicsModel> dummyModel;
    std::shared_ptr<AMP::Operator::DirichletVectorCorrection> dirichletVecOp =
        std::dynamic_pointer_cast<AMP::Operator::DirichletVectorCorrection>(
            AMP::Operator::OperatorBuilder::createOperator(
                meshAdapter, "Load_Boundary", input_db, dummyModel ) );
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

    mechSolVec->setToScalar( 0.5, mechSolVec );
    mechRhsVec->setToScalar( 0.0, mechRhsVec );
    mechResVec->setToScalar( 0.0, mechResVec );

    dirichletVecOp->apply( nullVec, mechRhsVec );

    double rhsNorm = mechRhsVec->L2Norm(mechRhsVec);

    std::cout << "RHS Norm: " << rhsNorm << std::endl;

    double initSolNorm = mechSolVec->L2Norm(mechSolVec);

    std::cout << "Initial Solution Norm: " << initSolNorm << std::endl;

    bvpOperator->residual( mechRhsVec, mechSolVec, mechResVec );

    double initResidualNorm = mechResVec->L2Norm(mechResVec);

    std::cout << "Initial Residual Norm: " << initResidualNorm << std::endl;

    std::shared_ptr<AMP::Database> mlSolver_db = input_db->getDatabase( "LinearSolver" );

    std::shared_ptr<AMP::Solver::SolverStrategyParameters> mlSolverParams(
        new AMP::Solver::SolverStrategyParameters( mlSolver_db ) );

    mlSolverParams->d_pOperator = bvpOperator;

    // create the ML solver interface
    std::shared_ptr<AMP::Solver::TrilinosMueLuSolver> mlSolver(
        new AMP::Solver::TrilinosMueLuSolver( mlSolverParams ) );

    mlSolver->setZeroInitialGuess( false );

    mlSolver->solve( mechRhsVec, mechSolVec );

#ifdef USE_EXT_SILO
    siloWriter->registerVector( mechSolVec, meshAdapter, AMP::Mesh::GeomType::Vertex, "Solution" );
    siloWriter->writeFile( exeName, 0 );
#endif

    bvpOperator->residual( mechRhsVec, mechSolVec, mechResVec );

    double finalResidualNorm = mechResVec->L2Norm(mechResVec);

    std::cout << "Final Residual Norm: " << finalResidualNorm << std::endl;

    if ( finalResidualNorm > ( 1e-10 * initResidualNorm ) ) {
        ut->failure(
            "TrilinosMueLuSolver could not solve a linear elasticity problem successfully" );
    } else {
        ut->passes( "TrilinosMueLuSolver successfully solves a linear elasticity problem" );
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
