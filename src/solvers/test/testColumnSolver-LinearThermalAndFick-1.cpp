#include "AMP/ampmesh/Mesh.h"
#include "AMP/discretization/DOF_Manager.h"
#include "AMP/discretization/simpleDOF_Manager.h"
#include "AMP/materials/Material.h"
#include "AMP/operators/BVPOperatorParameters.h"
#include "AMP/operators/ColumnOperator.h"
#include "AMP/operators/LinearBVPOperator.h"
#include "AMP/operators/OperatorBuilder.h"
#include "AMP/operators/boundary/DirichletVectorCorrection.h"
#include "AMP/operators/diffusion/DiffusionLinearFEOperator.h"
#include "AMP/solvers/ColumnSolver.h"
#include "AMP/solvers/trilinos/ml/TrilinosMLSolver.h"
#include "AMP/utils/AMPManager.h"
#include "AMP/utils/AMP_MPI.h"
#include "AMP/utils/Database.h"
#include "AMP/utils/PIO.h"
#include "AMP/utils/UnitTest.h"
#include "AMP/utils/Utilities.h"
#include "AMP/utils/Writer.h"
#include "AMP/vectors/MultiVariable.h"
#include "AMP/vectors/VectorBuilder.h"
#include <memory>

#include <iostream>
#include <string>


void myTest( AMP::UnitTest *ut, const std::string &exeName )
{
    std::string input_file = "input_" + exeName;
    std::string log_file   = "output_" + exeName;

    AMP::PIO::logOnlyNodeZero( log_file );


    auto input_db = AMP::Database::parseInputFile( input_file );
    input_db->print( AMP::plog );

    //--------------------------------------------------
    //   Create the Mesh.
    //--------------------------------------------------
    AMP_INSIST( input_db->keyExists( "Mesh" ), "Key ''Mesh'' is missing!" );
    std::shared_ptr<AMP::Database> mesh_db = input_db->getDatabase( "Mesh" );
    std::shared_ptr<AMP::Mesh::MeshParameters> mgrParams(
        new AMP::Mesh::MeshParameters( mesh_db ) );
    mgrParams->setComm( AMP::AMP_MPI( AMP_COMM_WORLD ) );
    std::shared_ptr<AMP::Mesh::Mesh> meshAdapter = AMP::Mesh::Mesh::buildMesh( mgrParams );
    //--------------------------------------------------

    //--------------------------------------------------
    // Create a DOF manager for a nodal vector
    //--------------------------------------------------
    int DOFsPerNode     = 1;
    int nodalGhostWidth = 1;
    bool split          = true;
    AMP::Discretization::DOFManager::shared_ptr nodalDofMap =
        AMP::Discretization::simpleDOFManager::create(
            meshAdapter, AMP::Mesh::GeomType::Vertex, nodalGhostWidth, DOFsPerNode, split );
    //--------------------------------------------------

    std::shared_ptr<AMP::Operator::ElementPhysicsModel> FickMaterialModel;
    std::shared_ptr<AMP::Operator::ElementPhysicsModel> thermalTransportModel;

    //----------------------------------------------------------------------------------------------------------------------------------------------//
    // now construct the linear BVP operator for Fick
    AMP_INSIST( input_db->keyExists( "testLinearFickOperator" ), "key missing!" );
    std::shared_ptr<AMP::Operator::LinearBVPOperator> linearFickOperator =
        std::dynamic_pointer_cast<AMP::Operator::LinearBVPOperator>(
            AMP::Operator::OperatorBuilder::createOperator(
                meshAdapter, "testLinearFickOperator", input_db, FickMaterialModel ) );

    //----------------------------------------------------------------------------------------------------------------------------------------------//
    // now construct the linear BVP operator for thermal
    AMP_INSIST( input_db->keyExists( "testLinearThermalOperator" ), "key missing!" );
    std::shared_ptr<AMP::Operator::LinearBVPOperator> linearThermalOperator =
        std::dynamic_pointer_cast<AMP::Operator::LinearBVPOperator>(
            AMP::Operator::OperatorBuilder::createOperator(
                meshAdapter, "testLinearThermalOperator", input_db, thermalTransportModel ) );

    //----------------------------------------------------------------------------------------------------------------------------------------------//
    // create a column operator object for linear Thermal-Fick
    std::shared_ptr<AMP::Operator::OperatorParameters> params;
    std::shared_ptr<AMP::Operator::ColumnOperator> linearThermalFickOperator(
        new AMP::Operator::ColumnOperator( params ) );
    linearThermalFickOperator->append( linearThermalOperator );
    linearThermalFickOperator->append( linearFickOperator );

    //----------------------------------------------------------------------------------------------------------------------------------------------//
    // initialize the input multi-variable
    std::shared_ptr<AMP::Operator::DiffusionLinearFEOperator> thermalVolumeOperator =
        std::dynamic_pointer_cast<AMP::Operator::DiffusionLinearFEOperator>(
            linearThermalOperator->getVolumeOperator() );
    std::shared_ptr<AMP::Operator::DiffusionLinearFEOperator> FickVolumeOperator =
        std::dynamic_pointer_cast<AMP::Operator::DiffusionLinearFEOperator>(
            linearFickOperator->getVolumeOperator() );

    std::shared_ptr<AMP::LinearAlgebra::MultiVariable> inputVariable(
        new AMP::LinearAlgebra::MultiVariable( "inputVariable" ) );
    inputVariable->add( thermalVolumeOperator->getInputVariable() );
    inputVariable->add( FickVolumeOperator->getInputVariable() );

    // initialize the output multi-variable
    AMP::LinearAlgebra::Variable::shared_ptr outputVariable =
        linearThermalFickOperator->getOutputVariable();

    // create solution, rhs, and residual vectors
    AMP::LinearAlgebra::Vector::shared_ptr solVec =
        AMP::LinearAlgebra::createVector( nodalDofMap, inputVariable );
    AMP::LinearAlgebra::Vector::shared_ptr rhsVec =
        AMP::LinearAlgebra::createVector( nodalDofMap, outputVariable );
    AMP::LinearAlgebra::Vector::shared_ptr resVec =
        AMP::LinearAlgebra::createVector( nodalDofMap, outputVariable );

    // create the following shared pointers for ease of use
    AMP::LinearAlgebra::Vector::shared_ptr nullVec;

    //----------------------------------------------------------------------------------------------------------------------------------------------//
    // Initial-Guess for thermal
    std::shared_ptr<AMP::Operator::ElementPhysicsModel> dummyThermalModel;
    std::shared_ptr<AMP::Operator::DirichletVectorCorrection> dirichletThermalInVecOp =
        std::dynamic_pointer_cast<AMP::Operator::DirichletVectorCorrection>(
            AMP::Operator::OperatorBuilder::createOperator(
                meshAdapter, "ThermalInitialGuess", input_db, dummyThermalModel ) );
    dirichletThermalInVecOp->setVariable( thermalVolumeOperator->getInputVariable() );

    //----------------------------------------------------------------------------------------------------------------------------------------------//
    // Initial-Guess for Fick
    std::shared_ptr<AMP::Operator::ElementPhysicsModel> dummyFickModel;
    std::shared_ptr<AMP::Operator::DirichletVectorCorrection> dirichletFickInVecOp =
        std::dynamic_pointer_cast<AMP::Operator::DirichletVectorCorrection>(
            AMP::Operator::OperatorBuilder::createOperator(
                meshAdapter, "FickInitialGuess", input_db, dummyFickModel ) );
    dirichletFickInVecOp->setVariable( FickVolumeOperator->getInputVariable() );

    //----------------------------------------------------------------------------------------------------------------------------------------------//
    // Random initial guess
    solVec->setRandomValues();

    // Initial guess for thermal must satisfy the thermal Dirichlet boundary conditions
    dirichletThermalInVecOp->apply( nullVec, solVec );

    // Initial guess for Fick must satisfy the Fick Dirichlet boundary conditions
    dirichletFickInVecOp->apply( nullVec, solVec );

    rhsVec->setToScalar( 0.0 );

    //----------------------------------------------------------------------------------------------------------------------------------------------//
    // ---- initialize the solvers
    std::shared_ptr<AMP::Database> columnSolver_db = input_db->getDatabase( "ColumnSolver" );
    std::shared_ptr<AMP::Solver::SolverStrategyParameters> columnSolverParams(
        new AMP::Solver::SolverStrategyParameters( columnSolver_db ) );
    std::shared_ptr<AMP::Solver::ColumnSolver> columnSolver(
        new AMP::Solver::ColumnSolver( columnSolverParams ) );

    std::shared_ptr<AMP::Database> thermalSolver_db =
        columnSolver_db->getDatabase( "ThermalSolver" );
    std::shared_ptr<AMP::Solver::SolverStrategyParameters> thermalSolverParams(
        new AMP::Solver::SolverStrategyParameters( thermalSolver_db ) );
    thermalSolverParams->d_pOperator = linearThermalOperator;
    std::shared_ptr<AMP::Solver::TrilinosMLSolver> linearThermalSolver(
        new AMP::Solver::TrilinosMLSolver( thermalSolverParams ) );

    std::shared_ptr<AMP::Database> FickSolver_db = columnSolver_db->getDatabase( "FickSolver" );
    std::shared_ptr<AMP::Solver::SolverStrategyParameters> FickSolverParams(
        new AMP::Solver::SolverStrategyParameters( FickSolver_db ) );
    FickSolverParams->d_pOperator = linearFickOperator;
    std::shared_ptr<AMP::Solver::TrilinosMLSolver> linearFickSolver(
        new AMP::Solver::TrilinosMLSolver( FickSolverParams ) );

    columnSolver->append( linearThermalSolver );
    columnSolver->append( linearFickSolver );

    double initialResidualNorm = resVec->L2Norm();

    AMP::pout << "Initial Residual Norm: " << initialResidualNorm << std::endl;

    columnSolver->solve( rhsVec, solVec );

    linearThermalFickOperator->residual( rhsVec, solVec, resVec );

    double finalResidualNorm = resVec->L2Norm();

    std::cout << "Final Residual Norm: " << finalResidualNorm << std::endl;

    if ( finalResidualNorm > 1.0e-08 ) {
        ut->failure( "ColumnSolver unsuccessfully solves two linear operators" );
    } else {
        ut->passes( "ColumnSolver successfully solves two linear operators" );
    }
    ut->passes( exeName );
}

int main( int argc, char *argv[] )
{
    AMP::AMPManager::startup( argc, argv );
    AMP::UnitTest ut;

    std::vector<std::string> exeNames;
    exeNames.emplace_back( "testColumnSolver-LinearThermalAndFick-1" );

    for ( auto &exeName : exeNames )
        myTest( &ut, exeName );

    ut.report();

    int num_failed = ut.NumFailGlobal();
    AMP::AMPManager::shutdown();
    return num_failed;
}
