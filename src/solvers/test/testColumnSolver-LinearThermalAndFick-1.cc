#include "utils/AMPManager.h"
#include "utils/UnitTest.h"
#include "utils/Utilities.h"
#include <iostream>
#include <string>

#include "utils/shared_ptr.h"

#include "materials/Material.h"
#include "utils/AMPManager.h"
#include "utils/AMP_MPI.h"
#include "utils/Database.h"
#include "utils/InputDatabase.h"
#include "utils/InputManager.h"
#include "utils/PIO.h"


#include "utils/Writer.h"

#include "ampmesh/Mesh.h"
#include "discretization/DOF_Manager.h"
#include "discretization/simpleDOF_Manager.h"
#include "vectors/MultiVariable.h"
#include "vectors/VectorBuilder.h"

#include "operators/boundary/DirichletVectorCorrection.h"

#include "operators/BVPOperatorParameters.h"
#include "operators/ColumnOperator.h"
#include "operators/LinearBVPOperator.h"
#include "operators/OperatorBuilder.h"
#include "operators/diffusion/DiffusionLinearFEOperator.h"

#include "solvers/ColumnSolver.h"
#include "solvers/trilinos/ml/TrilinosMLSolver.h"


void myTest( AMP::UnitTest *ut, const std::string& exeName )
{
    std::string input_file = "input_" + exeName;
    std::string log_file   = "output_" + exeName;

    AMP::PIO::logOnlyNodeZero( log_file );

    AMP::shared_ptr<AMP::InputDatabase> input_db( new AMP::InputDatabase( "input_db" ) );
    AMP::InputManager::getManager()->parseInputFile( input_file, input_db );
    input_db->printClassData( AMP::plog );

    //--------------------------------------------------
    //   Create the Mesh.
    //--------------------------------------------------
    AMP_INSIST( input_db->keyExists( "Mesh" ), "Key ''Mesh'' is missing!" );
    AMP::shared_ptr<AMP::Database> mesh_db = input_db->getDatabase( "Mesh" );
    AMP::shared_ptr<AMP::Mesh::MeshParameters> mgrParams(
        new AMP::Mesh::MeshParameters( mesh_db ) );
    mgrParams->setComm( AMP::AMP_MPI( AMP_COMM_WORLD ) );
    AMP::shared_ptr<AMP::Mesh::Mesh> meshAdapter = AMP::Mesh::Mesh::buildMesh( mgrParams );
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

    AMP::shared_ptr<AMP::Operator::ElementPhysicsModel> FickMaterialModel;
    AMP::shared_ptr<AMP::Operator::ElementPhysicsModel> thermalTransportModel;

    //----------------------------------------------------------------------------------------------------------------------------------------------//
    // now construct the linear BVP operator for Fick
    AMP_INSIST( input_db->keyExists( "testLinearFickOperator" ), "key missing!" );
    AMP::shared_ptr<AMP::Operator::LinearBVPOperator> linearFickOperator =
        AMP::dynamic_pointer_cast<AMP::Operator::LinearBVPOperator>(
            AMP::Operator::OperatorBuilder::createOperator(
                meshAdapter, "testLinearFickOperator", input_db, FickMaterialModel ) );

    //----------------------------------------------------------------------------------------------------------------------------------------------//
    // now construct the linear BVP operator for thermal
    AMP_INSIST( input_db->keyExists( "testLinearThermalOperator" ), "key missing!" );
    AMP::shared_ptr<AMP::Operator::LinearBVPOperator> linearThermalOperator =
        AMP::dynamic_pointer_cast<AMP::Operator::LinearBVPOperator>(
            AMP::Operator::OperatorBuilder::createOperator(
                meshAdapter, "testLinearThermalOperator", input_db, thermalTransportModel ) );

    //----------------------------------------------------------------------------------------------------------------------------------------------//
    // create a column operator object for linear Thermal-Fick
    AMP::shared_ptr<AMP::Operator::OperatorParameters> params;
    AMP::shared_ptr<AMP::Operator::ColumnOperator> linearThermalFickOperator(
        new AMP::Operator::ColumnOperator( params ) );
    linearThermalFickOperator->append( linearThermalOperator );
    linearThermalFickOperator->append( linearFickOperator );

    //----------------------------------------------------------------------------------------------------------------------------------------------//
    // initialize the input multi-variable
    AMP::shared_ptr<AMP::Operator::DiffusionLinearFEOperator> thermalVolumeOperator =
        AMP::dynamic_pointer_cast<AMP::Operator::DiffusionLinearFEOperator>(
            linearThermalOperator->getVolumeOperator() );
    AMP::shared_ptr<AMP::Operator::DiffusionLinearFEOperator> FickVolumeOperator =
        AMP::dynamic_pointer_cast<AMP::Operator::DiffusionLinearFEOperator>(
            linearFickOperator->getVolumeOperator() );

    AMP::shared_ptr<AMP::LinearAlgebra::MultiVariable> inputVariable(
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
    AMP::shared_ptr<AMP::Operator::ElementPhysicsModel> dummyThermalModel;
    AMP::shared_ptr<AMP::Operator::DirichletVectorCorrection> dirichletThermalInVecOp =
        AMP::dynamic_pointer_cast<AMP::Operator::DirichletVectorCorrection>(
            AMP::Operator::OperatorBuilder::createOperator(
                meshAdapter, "ThermalInitialGuess", input_db, dummyThermalModel ) );
    dirichletThermalInVecOp->setVariable( thermalVolumeOperator->getInputVariable() );

    //----------------------------------------------------------------------------------------------------------------------------------------------//
    // Initial-Guess for Fick
    AMP::shared_ptr<AMP::Operator::ElementPhysicsModel> dummyFickModel;
    AMP::shared_ptr<AMP::Operator::DirichletVectorCorrection> dirichletFickInVecOp =
        AMP::dynamic_pointer_cast<AMP::Operator::DirichletVectorCorrection>(
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
    AMP::shared_ptr<AMP::Database> columnSolver_db = input_db->getDatabase( "ColumnSolver" );
    AMP::shared_ptr<AMP::Solver::SolverStrategyParameters> columnSolverParams(
        new AMP::Solver::SolverStrategyParameters( columnSolver_db ) );
    AMP::shared_ptr<AMP::Solver::ColumnSolver> columnSolver(
        new AMP::Solver::ColumnSolver( columnSolverParams ) );

    AMP::shared_ptr<AMP::Database> thermalSolver_db =
        columnSolver_db->getDatabase( "ThermalSolver" );
    AMP::shared_ptr<AMP::Solver::SolverStrategyParameters> thermalSolverParams(
        new AMP::Solver::SolverStrategyParameters( thermalSolver_db ) );
    thermalSolverParams->d_pOperator = linearThermalOperator;
    AMP::shared_ptr<AMP::Solver::TrilinosMLSolver> linearThermalSolver(
        new AMP::Solver::TrilinosMLSolver( thermalSolverParams ) );

    AMP::shared_ptr<AMP::Database> FickSolver_db = columnSolver_db->getDatabase( "FickSolver" );
    AMP::shared_ptr<AMP::Solver::SolverStrategyParameters> FickSolverParams(
        new AMP::Solver::SolverStrategyParameters( FickSolver_db ) );
    FickSolverParams->d_pOperator = linearFickOperator;
    AMP::shared_ptr<AMP::Solver::TrilinosMLSolver> linearFickSolver(
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
