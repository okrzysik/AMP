#include "AMP/IO/PIO.h"
#include "AMP/IO/Writer.h"
#include "AMP/discretization/simpleDOF_Manager.h"
#include "AMP/mesh/MeshFactory.h"
#include "AMP/mesh/MeshParameters.h"
#include "AMP/operators/ColumnOperator.h"
#include "AMP/operators/LinearBVPOperator.h"
#include "AMP/operators/OperatorBuilder.h"
#include "AMP/operators/mechanics/ConstructLinearMechanicsRHSVector.h"
#include "AMP/solvers/SolverStrategyParameters.h"
#include "AMP/solvers/petsc/PetscKrylovSolver.h"
#include "AMP/solvers/trilinos/ml/TrilinosMLSolver.h"
#include "AMP/utils/AMPManager.h"
#include "AMP/utils/Database.h"
#include "AMP/utils/UnitTest.h"
#include "AMP/vectors/VectorBuilder.h"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <memory>
#include <string>
#include <vector>


static void myTest( AMP::UnitTest *ut, const std::string &exeName )
{
    std::string input_file = "input_" + exeName;
    std::string log_file   = "output_" + exeName;

    AMP::logOnlyNodeZero( log_file );

    auto input_db = AMP::Database::parseInputFile( input_file );
    input_db->print( AMP::plog );

    AMP_INSIST( input_db->keyExists( "Mesh" ), "Key ''Mesh'' is missing!" );
    auto mesh_db    = input_db->getDatabase( "Mesh" );
    auto meshParams = std::make_shared<AMP::Mesh::MeshParameters>( mesh_db );
    meshParams->setComm( AMP::AMP_MPI( AMP_COMM_WORLD ) );
    auto meshAdapter = AMP::Mesh::MeshFactory::create( meshParams );

    std::shared_ptr<AMP::Operator::ElementPhysicsModel> elementPhysicsModel;
    auto bvpOperator = std::dynamic_pointer_cast<AMP::Operator::LinearBVPOperator>(
        AMP::Operator::OperatorBuilder::createOperator(
            meshAdapter, "MechanicsBVPOperator", input_db, elementPhysicsModel ) );

    AMP::pout << "Constructed BVP operator" << std::endl;

    auto dispVar = bvpOperator->getOutputVariable();
    auto tempVar = std::make_shared<AMP::LinearAlgebra::Variable>( "temperature" );

    auto tempDofMap = AMP::Discretization::simpleDOFManager::create(
        meshAdapter, AMP::Mesh::GeomType::Vertex, 1, 1, true );

    auto dispDofMap = AMP::Discretization::simpleDOFManager::create(
        meshAdapter, AMP::Mesh::GeomType::Vertex, 1, 3, true );

    AMP::LinearAlgebra::Vector::shared_ptr nullVec;
    auto mechSolVec = AMP::LinearAlgebra::createVector( dispDofMap, dispVar, true );
    auto mechRhsVec = mechSolVec->cloneVector();
    auto mechResVec = mechSolVec->cloneVector();

    auto currTempVec = AMP::LinearAlgebra::createVector( tempDofMap, tempVar, true );
    auto prevTempVec = currTempVec->cloneVector();

    mechSolVec->setToScalar( 0.0 );
    mechResVec->setToScalar( 0.0 );

    currTempVec->setToScalar( 500.0 );
    prevTempVec->setToScalar( 300.0 );

    auto temperatureRhsDatabase = input_db->getDatabase( "TemperatureRHS" );

    computeTemperatureRhsVector( meshAdapter,
                                 temperatureRhsDatabase,
                                 tempVar,
                                 dispVar,
                                 currTempVec,
                                 prevTempVec,
                                 mechRhsVec );
    bvpOperator->modifyRHSvector( mechRhsVec );

    auto linearSolver_db = input_db->getDatabase( "LinearSolver" );
    auto pcSolver_db     = linearSolver_db->getDatabase( "Preconditioner" );
    auto pcSolverParams  = std::make_shared<AMP::Solver::TrilinosMLSolverParameters>( pcSolver_db );
    pcSolverParams->d_pOperator = bvpOperator;
    auto pcSolver               = std::make_shared<AMP::Solver::TrilinosMLSolver>( pcSolverParams );

    auto linearSolverParams =
        std::make_shared<AMP::Solver::SolverStrategyParameters>( linearSolver_db );
    linearSolverParams->d_pOperator     = bvpOperator;
    linearSolverParams->d_comm          = AMP_COMM_WORLD;
    linearSolverParams->d_pNestedSolver = pcSolver;
    auto linearSolver = std::make_shared<AMP::Solver::PetscKrylovSolver>( linearSolverParams );

    linearSolver->apply( mechRhsVec, mechSolVec );

    // Create the silo writer and register the data
    auto siloWriter = AMP::IO::Writer::buildWriter( "Silo" );
    siloWriter->registerVector( mechSolVec, meshAdapter, AMP::Mesh::GeomType::Vertex, "Solution" );
    siloWriter->writeFile( exeName, 1 );
    meshAdapter->displaceMesh( mechSolVec );
    siloWriter->writeFile( exeName, 2 );

    ut->passes( exeName );
}

int testLinearThermalExpansion( int argc, char *argv[] )
{
    AMP::AMPManager::startup( argc, argv );
    AMP::UnitTest ut;

    std::string exeName = "testLinearThermalExpansion";

    myTest( &ut, exeName );

    ut.report();

    int num_failed = ut.NumFailGlobal();
    AMP::AMPManager::shutdown();
    return num_failed;
}
