#include "AMP/discretization/simpleDOF_Manager.h"
#include "AMP/operators/BVPOperatorParameters.h"
#include "AMP/operators/LinearBVPOperator.h"
#include "AMP/operators/NonlinearBVPOperator.h"
#include "AMP/operators/OperatorBuilder.h"
#include "AMP/operators/boundary/DirichletMatrixCorrection.h"
#include "AMP/operators/boundary/DirichletVectorCorrection.h"
#include "AMP/operators/mechanics/MechanicsLinearElement.h"
#include "AMP/operators/mechanics/MechanicsLinearFEOperator.h"
#include "AMP/operators/mechanics/MechanicsNonlinearElement.h"
#include "AMP/operators/mechanics/MechanicsNonlinearFEOperator.h"
#include "AMP/operators/mechanics/ThermalStrainMaterialModel.h"
#include "AMP/solvers/petsc/PetscKrylovSolver.h"
#include "AMP/solvers/petsc/PetscKrylovSolverParameters.h"
#include "AMP/solvers/petsc/PetscSNESSolver.h"
#include "AMP/solvers/petsc/PetscSNESSolverParameters.h"
#include "AMP/solvers/trilinos/ml/TrilinosMLSolver.h"
#include "AMP/utils/AMPManager.h"
#include "AMP/utils/Database.h"
#include "AMP/utils/InputDatabase.h"
#include "AMP/utils/InputManager.h"
#include "AMP/utils/PIO.h"
#include "AMP/utils/UnitTest.h"
#include "AMP/utils/Utilities.h"
#include "AMP/utils/Writer.h"
#include "AMP/vectors/VectorBuilder.h"

#include <iostream>
#include <string>

static void myTest( AMP::UnitTest *ut, const std::string &exeName )
{
    std::string input_file  = "input_" + exeName;
    std::string output_file = "output_" + exeName + ".txt";
    std::string log_file    = "log_" + exeName;

    AMP::PIO::logOnlyNodeZero( log_file );
    AMP::AMP_MPI globalComm( AMP_COMM_WORLD );

    // Read the input file
    AMP::shared_ptr<AMP::InputDatabase> input_db( new AMP::InputDatabase( "input_db" ) );
    AMP::InputManager::getManager()->parseInputFile( input_file, input_db );
    input_db->printClassData( AMP::plog );

    // Read the mesh
    AMP_INSIST( input_db->keyExists( "Mesh" ), "Key ''Mesh'' is missing!" );
    AMP::shared_ptr<AMP::Database> mesh_db = input_db->getDatabase( "Mesh" );
    AMP::shared_ptr<AMP::Mesh::MeshParameters> meshParams(
        new AMP::Mesh::MeshParameters( mesh_db ) );
    meshParams->setComm( AMP::AMP_MPI( AMP_COMM_WORLD ) );
    AMP::Mesh::Mesh::shared_ptr meshAdapter = AMP::Mesh::Mesh::buildMesh( meshParams );

    // Create a nonlinear BVP operator for mechanics
    AMP_INSIST( input_db->keyExists( "NonlinearMechanicsOperator" ), "key missing!" );
    AMP::shared_ptr<AMP::Operator::NonlinearBVPOperator> nonlinearMechanicsBVPoperator =
        AMP::dynamic_pointer_cast<AMP::Operator::NonlinearBVPOperator>(
            AMP::Operator::OperatorBuilder::createOperator(
                meshAdapter, "NonlinearMechanicsOperator", input_db ) );
    AMP::shared_ptr<AMP::Operator::MechanicsNonlinearFEOperator> nonlinearMechanicsVolumeOperator =
        AMP::dynamic_pointer_cast<AMP::Operator::MechanicsNonlinearFEOperator>(
            nonlinearMechanicsBVPoperator->getVolumeOperator() );
    AMP::shared_ptr<AMP::Operator::ElementPhysicsModel> mechanicsMaterialModel =
        nonlinearMechanicsVolumeOperator->getMaterialModel();

    // Create a Linear BVP operator for mechanics
    AMP_INSIST( input_db->keyExists( "LinearMechanicsOperator" ), "key missing!" );
    AMP::shared_ptr<AMP::Operator::LinearBVPOperator> linearMechanicsBVPoperator =
        AMP::dynamic_pointer_cast<AMP::Operator::LinearBVPOperator>(
            AMP::Operator::OperatorBuilder::createOperator(
                meshAdapter, "LinearMechanicsOperator", input_db, mechanicsMaterialModel ) );

    // Create the variables
    AMP::shared_ptr<AMP::Operator::MechanicsNonlinearFEOperator> mechanicsNonlinearVolumeOperator =
        AMP::dynamic_pointer_cast<AMP::Operator::MechanicsNonlinearFEOperator>(
            nonlinearMechanicsBVPoperator->getVolumeOperator() );
    AMP::LinearAlgebra::Variable::shared_ptr dispVar =
        mechanicsNonlinearVolumeOperator->getOutputVariable();
    AMP::LinearAlgebra::Variable::shared_ptr tempVar( new AMP::LinearAlgebra::Variable( "temp" ) );

    AMP::Discretization::DOFManager::shared_ptr tempDofMap =
        AMP::Discretization::simpleDOFManager::create(
            meshAdapter, AMP::Mesh::GeomType::Vertex, 1, 1, true );

    AMP::Discretization::DOFManager::shared_ptr dispDofMap =
        AMP::Discretization::simpleDOFManager::create(
            meshAdapter, AMP::Mesh::GeomType::Vertex, 1, 3, true );

    // Create the vectors
    AMP::LinearAlgebra::Vector::shared_ptr nullVec;
    AMP::LinearAlgebra::Vector::shared_ptr solVec =
        AMP::LinearAlgebra::createVector( dispDofMap, dispVar, true );
    AMP::LinearAlgebra::Vector::shared_ptr rhsVec = solVec->cloneVector();
    AMP::LinearAlgebra::Vector::shared_ptr resVec = solVec->cloneVector();
    AMP::LinearAlgebra::Vector::shared_ptr initTempVec =
        AMP::LinearAlgebra::createVector( tempDofMap, tempVar, true );
    AMP::LinearAlgebra::Vector::shared_ptr finalTempVec = initTempVec->cloneVector();

    // Set Initial Temperature
    AMP_INSIST( input_db->keyExists( "INIT_TEMP_CONST" ), "key missing!" );
    double initTempVal = input_db->getDouble( "INIT_TEMP_CONST" );
    initTempVec->setToScalar( initTempVal );
    mechanicsNonlinearVolumeOperator->setReferenceTemperature( initTempVec );

    // Set Final Temperature
    AMP_INSIST( input_db->keyExists( "FINAL_TEMP_CONST" ), "key missing!" );
    double finalTempVal = input_db->getDouble( "FINAL_TEMP_CONST" );
    finalTempVec->setToScalar( finalTempVal );
    mechanicsNonlinearVolumeOperator->setVector( AMP::Operator::Mechanics::TEMPERATURE,
                                                 finalTempVec );

    // Initial guess
    solVec->zero();
    nonlinearMechanicsBVPoperator->modifyInitialSolutionVector( solVec );

    // RHS
    rhsVec->zero();
    nonlinearMechanicsBVPoperator->modifyRHSvector( rhsVec );

    AMP::shared_ptr<AMP::Database> nonlinearSolver_db = input_db->getDatabase( "NonlinearSolver" );
    AMP::shared_ptr<AMP::Solver::PetscSNESSolver> nonlinearSolver;

    AMP_INSIST( nonlinearSolver_db->keyExists( "LinearSolver" ), "key missing!" );
    AMP::shared_ptr<AMP::Database> linearSolver_db =
        nonlinearSolver_db->getDatabase( "LinearSolver" );
    AMP::shared_ptr<AMP::Solver::PetscKrylovSolver> linearSolver;

    AMP_INSIST( nonlinearSolver_db->keyExists( "usesJacobian" ), "key missing!" );
    bool usesJacobian = nonlinearSolver_db->getBool( "usesJacobian" );

    if ( usesJacobian ) {
        AMP::shared_ptr<AMP::Solver::PetscKrylovSolverParameters> linearSolverParams(
            new AMP::Solver::PetscKrylovSolverParameters( linearSolver_db ) );
        linearSolverParams->d_pOperator = linearMechanicsBVPoperator;
        linearSolverParams->d_comm      = globalComm;
        linearSolver.reset( new AMP::Solver::PetscKrylovSolver( linearSolverParams ) );
    }

    AMP::shared_ptr<AMP::Solver::PetscSNESSolverParameters> nonlinearSolverParams(
        new AMP::Solver::PetscSNESSolverParameters( nonlinearSolver_db ) );
    nonlinearSolverParams->d_comm          = globalComm;
    nonlinearSolverParams->d_pOperator     = nonlinearMechanicsBVPoperator;
    nonlinearSolverParams->d_pInitialGuess = solVec;
    if ( usesJacobian ) {
        nonlinearSolverParams->d_pKrylovSolver = linearSolver;
    }
    nonlinearSolver.reset( new AMP::Solver::PetscSNESSolver( nonlinearSolverParams ) );
    nonlinearSolver->setZeroInitialGuess( false );

    if ( !usesJacobian ) {
        linearSolver = nonlinearSolver->getKrylovSolver();
    }

    // We need to reset the linear operator before the solve since TrilinosML does
    // the factorization of the matrix during construction and so the matrix must
    // be correct before constructing the TrilinosML object.
    nonlinearMechanicsBVPoperator->apply( solVec, resVec );
    linearMechanicsBVPoperator->reset(
        nonlinearMechanicsBVPoperator->getParameters( "Jacobian", solVec ) );

    AMP_INSIST( linearSolver_db->keyExists( "Preconditioner" ), "key missing!" );
    AMP::shared_ptr<AMP::Database> preconditioner_db =
        linearSolver_db->getDatabase( "Preconditioner" );
    AMP::shared_ptr<AMP::Solver::SolverStrategyParameters> preconditionerParams(
        new AMP::Solver::SolverStrategyParameters( preconditioner_db ) );
    preconditionerParams->d_pOperator = linearMechanicsBVPoperator;
    AMP::shared_ptr<AMP::Solver::TrilinosMLSolver> preconditioner(
        new AMP::Solver::TrilinosMLSolver( preconditionerParams ) );

    linearSolver->setPreconditioner( preconditioner );

    nonlinearSolver->solve( rhsVec, solVec );

    mechanicsNonlinearVolumeOperator->printStressAndStrain( solVec, output_file );

    ut->passes( "testElasticThermalPatch" );
}


int testElasticThermalPatch( int argc, char *argv[] )
{
    AMP::AMPManager::startup( argc, argv );
    AMP::UnitTest ut;

    std::vector<std::string> mesh_files;
    mesh_files.emplace_back( "testElasticThermalPatch-1" );
    mesh_files.emplace_back( "testElasticThermalPatch-2" );
    mesh_files.emplace_back( "testElasticThermalPatch-3" );
    mesh_files.emplace_back( "testElasticThermalPatch-4" );

    for ( auto &mesh_file : mesh_files )
        myTest( &ut, mesh_file );

    ut.report();

    int num_failed = ut.NumFailGlobal();
    AMP::AMPManager::shutdown();
    return num_failed;
}
