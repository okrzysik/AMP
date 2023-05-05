#include "AMP/IO/PIO.h"
#include "AMP/IO/Writer.h"
#include "AMP/discretization/DOF_Manager.h"
#include "AMP/discretization/simpleDOF_Manager.h"
#include "AMP/matrices/Matrix.h"
#include "AMP/matrices/MatrixBuilder.h"
#include "AMP/matrices/trilinos/EpetraMatrix.h"
#include "AMP/mesh/Mesh.h"
#include "AMP/mesh/MeshFactory.h"
#include "AMP/mesh/MeshParameters.h"
#include "AMP/operators/BlockOperator.h"
#include "AMP/operators/LinearBVPOperator.h"
#include "AMP/operators/NeutronicsRhs.h"
#include "AMP/operators/OperatorBuilder.h"
#include "AMP/operators/libmesh/VolumeIntegralOperator.h"
#include "AMP/operators/petsc/PetscMatrixShellOperator.h"
#include "AMP/operators/trilinos/EpetraMatrixOperator.h"
#include "AMP/operators/trilinos/EpetraMatrixOperatorParameters.h"
#include "AMP/operators/trilinos/TrilinosMatrixShellOperator.h"
#include "AMP/solvers/SolverStrategyParameters.h"
#include "AMP/solvers/petsc/PetscKrylovSolver.h"
#include "AMP/solvers/petsc/PetscSNESSolver.h"
#include "AMP/solvers/trilinos/ml/TrilinosMLSolver.h"
#include "AMP/utils/AMPManager.h"
#include "AMP/utils/Database.h"
#include "AMP/utils/UnitTest.h"
#include "AMP/vectors/MultiVector.h"
#include "AMP/vectors/Variable.h"
#include "AMP/vectors/Vector.h"
#include "AMP/vectors/VectorBuilder.h"

#include <memory>
#include <string>


#define ITFAILS ut.failure( __LINE__ );
#define UNIT_TEST( a ) \
    if ( !( a ) )      \
        ut.failure( __LINE__ );

double h_Size = 1.46286;

static void myTest( AMP::UnitTest *ut, const std::string &exeName )

{

    std::string input_file = "input_" + exeName;

    // Read the input file
    auto input_db           = AMP::Database::parseInputFile( input_file );
    AMP::AMP_MPI globalComm = AMP::AMP_MPI( AMP_COMM_WORLD );
    input_db->print( AMP::plog );

    // Get the Mesh database and create the mesh parameters
    auto database = input_db->getDatabase( "Mesh" );
    auto params   = std::make_shared<AMP::Mesh::MeshParameters>( database );
    params->setComm( globalComm );

    // Create the meshes from the input database
    auto manager        = AMP::Mesh::MeshFactory::create( params );
    auto meshAdapterH27 = manager->Subset( "cubeH27" );
    auto meshAdapterH08 = manager->Subset( "cubeH08" );

    // Create the Conservation of Momentum Operator
    std::shared_ptr<AMP::Operator::ElementPhysicsModel> FlowTransportModel;
    AMP_INSIST( input_db->keyExists( "ConsMomentumLinearFEOperator" ), "key missing!" );
    auto ConsMomentumOperator = std::dynamic_pointer_cast<AMP::Operator::LinearBVPOperator>(
        AMP::Operator::OperatorBuilder::createOperator(
            meshAdapterH27, "ConsMomentumLinearBVPOperator", input_db, FlowTransportModel ) );

    // Create the Conservation of Mass Operator
    AMP_INSIST( input_db->keyExists( "ConsMassLinearFEOperator" ), "key missing!" );
    auto ConsMassOperator = std::dynamic_pointer_cast<AMP::Operator::LinearBVPOperator>(
        AMP::Operator::OperatorBuilder::createOperator(
            meshAdapterH08, "ConsMassLinearBVPOperator", input_db, FlowTransportModel ) );
    AMP::pout << "Finished creating Mass Operator" << std::endl;

    // Create the variables
    auto velocityVar = ConsMomentumOperator->getOutputVariable();
    auto pressureVar = ConsMassOperator->getOutputVariable();

    // Create the DOF managers
    auto DOF_scalar = AMP::Discretization::simpleDOFManager::create(
        manager, AMP::Mesh::GeomType::Vertex, 1, 1, true );
    auto DOF_vector = AMP::Discretization::simpleDOFManager::create(
        manager, AMP::Mesh::GeomType::Vertex, 1, 3, true );

    // Create the vectors
    auto globalSolVec =
        AMP::LinearAlgebra::MultiVector::create( "globalSolVec", manager->getComm() );
    globalSolVec->addVector( AMP::LinearAlgebra::createVector( DOF_scalar, pressureVar ) );
    globalSolVec->addVector( AMP::LinearAlgebra::createVector( DOF_vector, velocityVar ) );
    auto globalRhsVec = globalSolVec->clone( "globalRhsVec" );

    // Get the matrices
    auto FMat  = ConsMomentumOperator->getMatrix();
    auto BMat  = ConsMassOperator->getMatrix();
    auto BtMat = BMat->transpose();

    // Create a zero matrix over meshAdapterH08
    auto DOF_H08_scalar = AMP::Discretization::simpleDOFManager::create(
        meshAdapterH08, AMP::Mesh::GeomType::Vertex, 1, 1, true );
    auto zeroMat = AMP::LinearAlgebra::createMatrix(
        AMP::LinearAlgebra::createVector( DOF_H08_scalar, pressureVar ),
        AMP::LinearAlgebra::createVector( DOF_H08_scalar, pressureVar ) );
    zeroMat->zero();

    std::shared_ptr<AMP::Database> dummy_db;
    auto dummyParams1 = std::make_shared<AMP::Operator::EpetraMatrixOperatorParameters>( dummy_db );
    dummyParams1->d_Matrix = &( std::dynamic_pointer_cast<AMP::LinearAlgebra::EpetraMatrix>( BtMat )
                                    ->getEpetra_CrsMatrix() );
    auto bTOperator        = std::make_shared<AMP::Operator::EpetraMatrixOperator>( dummyParams1 );
    bTOperator->setVariables( ConsMassOperator->getOutputVariable(),
                              ConsMassOperator->getInputVariable() );

    auto dummyParams2 = std::make_shared<AMP::Operator::EpetraMatrixOperatorParameters>( dummy_db );
    dummyParams2->d_Matrix =
        &( std::dynamic_pointer_cast<AMP::LinearAlgebra::EpetraMatrix>( zeroMat )
               ->getEpetra_CrsMatrix() );
    auto zeroOperator = std::make_shared<AMP::Operator::EpetraMatrixOperator>( dummyParams2 );
    zeroOperator->setVariables( ConsMassOperator->getOutputVariable(),
                                ConsMassOperator->getOutputVariable() );

    AMP_INSIST( input_db->keyExists( "LinearSolver" ), "Key ''LinearSolver'' is missing!" );
    auto mlSolver_db = input_db->getDatabase( "LinearSolver" );

    auto convdiffSolverParams =
        std::make_shared<AMP::Solver::SolverStrategyParameters>( mlSolver_db );
    convdiffSolverParams->d_pOperator = ConsMomentumOperator;
    auto convdiffSolver = std::make_shared<AMP::Solver::TrilinosMLSolver>( convdiffSolverParams );
    convdiffSolver->setZeroInitialGuess( false );

    auto diagonalVec    = FMat->extractDiagonal();
    auto diagonalInvVec = diagonalVec->clone();
    diagonalInvVec->reciprocal( *diagonalVec );

    auto DMat = FMat->clone();
    DMat->zero();
    DMat->setDiagonal( diagonalVec );

    auto DInvMat = FMat->clone();
    DInvMat->zero();
    DInvMat->setDiagonal( diagonalInvVec );

    auto DInvBtMat = AMP::LinearAlgebra::Matrix::matMultiply( DInvMat, BtMat );

    auto schurMat = AMP::LinearAlgebra::Matrix::matMultiply( BtMat, DInvBtMat );

    auto dummyParams3 = std::make_shared<AMP::Operator::EpetraMatrixOperatorParameters>( dummy_db );
    dummyParams3->d_Matrix =
        &( std::dynamic_pointer_cast<AMP::LinearAlgebra::EpetraMatrix>( schurMat )
               ->getEpetra_CrsMatrix() );
    auto schurMatOperator = std::make_shared<AMP::Operator::EpetraMatrixOperator>( dummyParams3 );
    schurMatOperator->setVariables( ConsMassOperator->getOutputVariable(),
                                    ConsMassOperator->getOutputVariable() );

    auto schurMatSolverParams =
        std::make_shared<AMP::Solver::SolverStrategyParameters>( mlSolver_db );
    schurMatSolverParams->d_pOperator = schurMatOperator;
    auto schurMatSolver = std::make_shared<AMP::Solver::TrilinosMLSolver>( schurMatSolverParams );
    schurMatSolver->setZeroInitialGuess( false );

    auto velocityRhsVec = globalRhsVec->subsetVectorForVariable( velocityVar );
    auto pressureRhsVec = globalRhsVec->subsetVectorForVariable( pressureVar );

    auto velocitySolVec = globalSolVec->subsetVectorForVariable( velocityVar );
    auto pressureSolVec = globalSolVec->subsetVectorForVariable( pressureVar );

    auto pressureUpdateVec = pressureSolVec->clone();
    auto velocityUpdateVec = velocitySolVec->clone();

    auto velocityPrimeVec = velocitySolVec->clone();

    // SIMPLE(Semi Implicit Method for Pressure Linked Equations) ALGORITHM

    // STEP 1 :
    BtMat->mult( pressureSolVec, velocityRhsVec );
    convdiffSolver->apply( velocityRhsVec, velocityPrimeVec );

    // STEP 2 :
    BMat->mult( velocityPrimeVec, pressureRhsVec );
    schurMatSolver->apply( pressureRhsVec, pressureUpdateVec );

    // STEP 3 :
    DInvBtMat->mult( pressureUpdateVec, velocityUpdateVec );
    velocitySolVec->subtract( *velocityPrimeVec, *velocityUpdateVec );

    // STEP 4 :
    pressureSolVec->add( *velocityPrimeVec, *pressureUpdateVec );

    ut->passes( "Ran to completion" );
}


int testLinearFlow( int argc, char *argv[] )
{
    AMP::AMPManager::startup( argc, argv );
    AMP::UnitTest ut;
    std::vector<std::string> exeNames;
    exeNames.emplace_back( "testLinearFlow-1" );

    for ( auto name : exeNames )
        myTest( &ut, name );

    ut.report();
    int num_failed = ut.NumFailGlobal();

    AMP::AMPManager::shutdown();
    return num_failed;
}
