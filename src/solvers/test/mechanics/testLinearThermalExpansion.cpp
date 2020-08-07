#include "AMP/discretization/simpleDOF_Manager.h"
#include "AMP/operators/ColumnOperator.h"
#include "AMP/operators/LinearBVPOperator.h"
#include "AMP/operators/OperatorBuilder.h"
#include "AMP/operators/mechanics/ConstructLinearMechanicsRHSVector.h"
#include "AMP/solvers/petsc/PetscKrylovSolver.h"
#include "AMP/solvers/petsc/PetscKrylovSolverParameters.h"
#include "AMP/solvers/trilinos/ml/TrilinosMLSolver.h"
#include "AMP/utils/AMPManager.h"
#include "AMP/utils/Database.h"
#include "AMP/utils/PIO.h"
#include "AMP/utils/UnitTest.h"
#include "AMP/utils/Utilities.h"
#include "AMP/utils/Writer.h"
#include "AMP/vectors/VectorBuilder.h"
#include <memory>

#include <algorithm>
#include <cmath>
#include <iostream>
#include <string>
#include <vector>


static void myTest( AMP::UnitTest *ut, const std::string &exeName )
{
    std::string input_file = "input_" + exeName;
    std::string log_file   = "output_" + exeName;

    AMP::PIO::logOnlyNodeZero( log_file );


    auto input_db = AMP::Database::parseInputFile( input_file );
    input_db->print( AMP::plog );

#ifdef USE_EXT_SILO
    // Create the silo writer and register the data
    AMP::Utilities::Writer::shared_ptr siloWriter = AMP::Utilities::Writer::buildWriter( "Silo" );
#endif

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

    AMP::pout << "Constructed BVP operator" << std::endl;

    AMP::LinearAlgebra::Variable::shared_ptr dispVar = bvpOperator->getOutputVariable();
    AMP::LinearAlgebra::Variable::shared_ptr tempVar( new AMP::LinearAlgebra::Variable( "temp" ) );

    AMP::Discretization::DOFManager::shared_ptr tempDofMap =
        AMP::Discretization::simpleDOFManager::create(
            meshAdapter, AMP::Mesh::GeomType::Vertex, 1, 1, true );

    AMP::Discretization::DOFManager::shared_ptr dispDofMap =
        AMP::Discretization::simpleDOFManager::create(
            meshAdapter, AMP::Mesh::GeomType::Vertex, 1, 3, true );

    AMP::LinearAlgebra::Vector::shared_ptr nullVec;
    AMP::LinearAlgebra::Vector::shared_ptr mechSolVec =
        AMP::LinearAlgebra::createVector( dispDofMap, dispVar, true );
    AMP::LinearAlgebra::Vector::shared_ptr mechRhsVec = mechSolVec->cloneVector();
    AMP::LinearAlgebra::Vector::shared_ptr mechResVec = mechSolVec->cloneVector();

    AMP::LinearAlgebra::Vector::shared_ptr currTempVec =
        AMP::LinearAlgebra::createVector( tempDofMap, tempVar, true );
    AMP::LinearAlgebra::Vector::shared_ptr prevTempVec = currTempVec->cloneVector();

    mechSolVec->setToScalar( 0.0, mechSolVec );
    mechResVec->setToScalar( 0.0, mechResVec );

    currTempVec->setToScalar( 500.0, currTempVec );
    prevTempVec->setToScalar( 300.0, prevTempVec );

    std::shared_ptr<AMP::Database> temperatureRhsDatabase =
        input_db->getDatabase( "TemperatureRHS" );

    computeTemperatureRhsVector( meshAdapter,
                                 temperatureRhsDatabase,
                                 tempVar,
                                 dispVar,
                                 currTempVec,
                                 prevTempVec,
                                 mechRhsVec );
    bvpOperator->modifyRHSvector( mechRhsVec );

    std::shared_ptr<AMP::Database> linearSolver_db = input_db->getDatabase( "LinearSolver" );
    std::shared_ptr<AMP::Database> pcSolver_db = linearSolver_db->getDatabase( "Preconditioner" );
    std::shared_ptr<AMP::Solver::TrilinosMLSolverParameters> pcSolverParams(
        new AMP::Solver::TrilinosMLSolverParameters( pcSolver_db ) );
    pcSolverParams->d_pOperator = bvpOperator;
    std::shared_ptr<AMP::Solver::TrilinosMLSolver> pcSolver(
        new AMP::Solver::TrilinosMLSolver( pcSolverParams ) );

    std::shared_ptr<AMP::Solver::PetscKrylovSolverParameters> linearSolverParams(
        new AMP::Solver::PetscKrylovSolverParameters( linearSolver_db ) );
    linearSolverParams->d_pOperator       = bvpOperator;
    linearSolverParams->d_comm            = AMP_COMM_WORLD;
    linearSolverParams->d_pPreconditioner = pcSolver;
    std::shared_ptr<AMP::Solver::PetscKrylovSolver> linearSolver(
        new AMP::Solver::PetscKrylovSolver( linearSolverParams ) );

    linearSolver->solve( mechRhsVec, mechSolVec );

#ifdef USE_EXT_SILO
    siloWriter->registerVector( mechSolVec, meshAdapter, AMP::Mesh::GeomType::Vertex, "Solution" );
#endif

#ifdef USE_EXT_SILO
    siloWriter->writeFile( exeName, 1 );
    meshAdapter->displaceMesh( mechSolVec );
    siloWriter->writeFile( exeName, 2 );
#endif

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
