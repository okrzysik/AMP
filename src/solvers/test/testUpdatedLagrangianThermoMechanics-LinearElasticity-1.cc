#include "discretization/simpleDOF_Manager.h"
#include "vectors/VectorBuilder.h"

#include "operators/mechanics/MechanicsLinearFEOperator.h"
#include "operators/mechanics/MechanicsNonlinearFEOperator.h"
#include "operators/mechanics/ThermalStrainMaterialModel.h"

#include "operators/boundary/DirichletMatrixCorrection.h"
#include "operators/boundary/DirichletVectorCorrection.h"

#include "operators/BVPOperatorParameters.h"
#include "operators/LinearBVPOperator.h"
#include "operators/NonlinearBVPOperator.h"
#include "operators/OperatorBuilder.h"

#include "solvers/petsc/PetscKrylovSolver.h"
#include "solvers/petsc/PetscKrylovSolverParameters.h"
#include "solvers/petsc/PetscSNESSolver.h"
#include "solvers/petsc/PetscSNESSolverParameters.h"
#include "solvers/trilinos/ml/TrilinosMLSolver.h"

#include "utils/ReadTestMesh.h"
#include "utils/AMPManager.h"
#include "utils/InputManager.h"
#include "utils/UnitTest.h"
#include "utils/Utilities.h"
#include "utils/Writer.h"
#include "utils/PIO.h"

#include <iostream>
#include <string>


void myTest( AMP::UnitTest *ut, std::string exeName )
{
    std::string input_file = "input_" + exeName;
    std::string log_file   = "log_" + exeName;

    AMP::PIO::logOnlyNodeZero( log_file );
    AMP::AMP_MPI globalComm( AMP_COMM_WORLD );

    // Read the input file
    AMP::shared_ptr<AMP::InputDatabase> input_db( new AMP::InputDatabase( "input_db" ) );
    AMP::InputManager::getManager()->parseInputFile( input_file, input_db );
    input_db->printClassData( AMP::plog );

    AMP_INSIST( input_db->keyExists( "Mesh" ), "Key ''Mesh'' is missing!" );
    AMP::shared_ptr<AMP::Database> mesh_db = input_db->getDatabase( "Mesh" );
    AMP::shared_ptr<AMP::Mesh::MeshParameters> meshParams(
        new AMP::Mesh::MeshParameters( mesh_db ) );
    meshParams->setComm( AMP::AMP_MPI( AMP_COMM_WORLD ) );
    AMP::Mesh::Mesh::shared_ptr meshAdapter = AMP::Mesh::Mesh::buildMesh( meshParams );

    AMP_INSIST( input_db->keyExists( "NumberOfLoadingSteps" ),
                "Key ''NumberOfLoadingSteps'' is missing!" );
    int NumberOfLoadingSteps = input_db->getInteger( "NumberOfLoadingSteps" );

    // Create a nonlinear BVP operator for mechanics
    AMP_INSIST( input_db->keyExists( "NonlinearMechanicsOperator" ), "key missing!" );
    AMP::shared_ptr<AMP::Operator::ElementPhysicsModel> mechanicsMaterialModel;
    AMP::shared_ptr<AMP::Operator::NonlinearBVPOperator> nonlinearMechanicsBVPoperator =
        AMP::dynamic_pointer_cast<AMP::Operator::NonlinearBVPOperator>(
            AMP::Operator::OperatorBuilder::createOperator(
                meshAdapter, "NonlinearMechanicsOperator", input_db, mechanicsMaterialModel ) );

    // Create the variables
    AMP::shared_ptr<AMP::Operator::MechanicsNonlinearFEOperator> mechanicsNonlinearVolumeOperator =
        AMP::dynamic_pointer_cast<AMP::Operator::MechanicsNonlinearFEOperator>(
            nonlinearMechanicsBVPoperator->getVolumeOperator() );
    AMP::LinearAlgebra::Variable::shared_ptr dispVar =
        mechanicsNonlinearVolumeOperator->getOutputVariable();
    AMP::shared_ptr<AMP::LinearAlgebra::MultiVariable> inputVars =
        AMP::dynamic_pointer_cast<AMP::LinearAlgebra::MultiVariable>(
            mechanicsNonlinearVolumeOperator->getInputVariable() );
    AMP::LinearAlgebra::Variable::shared_ptr tempVar =
        inputVars->getVariable( AMP::Operator::Mechanics::TEMPERATURE );

    // For RHS (Point Forces)
    AMP::shared_ptr<AMP::Operator::ElementPhysicsModel> dummyModel;
    AMP::shared_ptr<AMP::Operator::DirichletVectorCorrection> dirichletLoadVecOp =
        AMP::dynamic_pointer_cast<AMP::Operator::DirichletVectorCorrection>(
            AMP::Operator::OperatorBuilder::createOperator(
                meshAdapter, "Load_Boundary", input_db, dummyModel ) );
    dirichletLoadVecOp->setVariable( dispVar );

    AMP::Discretization::DOFManager::shared_ptr dispDofMap =
        AMP::Discretization::simpleDOFManager::create( meshAdapter, AMP::Mesh::GeomType::Vertex, 1, 3, true );

    AMP::Discretization::DOFManager::shared_ptr tempDofMap =
        AMP::Discretization::simpleDOFManager::create( meshAdapter, AMP::Mesh::GeomType::Vertex, 1, 1, true );

    // Create the vectors
    AMP::LinearAlgebra::Vector::shared_ptr nullVec;
    AMP::LinearAlgebra::Vector::shared_ptr solVec =
        AMP::LinearAlgebra::createVector( dispDofMap, dispVar, true );
    AMP::LinearAlgebra::Vector::shared_ptr rhsVec       = solVec->cloneVector();
    AMP::LinearAlgebra::Vector::shared_ptr resVec       = solVec->cloneVector();
    AMP::LinearAlgebra::Vector::shared_ptr scaledRhsVec = solVec->cloneVector();
    AMP::LinearAlgebra::Vector::shared_ptr refTempVec =
        AMP::LinearAlgebra::createVector( tempDofMap, tempVar, true );
    AMP::LinearAlgebra::Vector::shared_ptr curTempVec = refTempVec->cloneVector();

    // Initial guess
    solVec->zero();
    nonlinearMechanicsBVPoperator->modifyInitialSolutionVector( solVec );

    // RHS
    rhsVec->zero();
    dirichletLoadVecOp->apply( nullVec, rhsVec );
    nonlinearMechanicsBVPoperator->modifyRHSvector( rhsVec );

    // Set the temperatures
    refTempVec->setToScalar( 300.0 );
    curTempVec->setToScalar( 700.0 );

    // Set the reference temperature
    ( AMP::dynamic_pointer_cast<AMP::Operator::MechanicsNonlinearFEOperator>(
          nonlinearMechanicsBVPoperator->getVolumeOperator() ) )
        ->setReferenceTemperature( refTempVec );
    ( AMP::dynamic_pointer_cast<AMP::Operator::MechanicsNonlinearFEOperator>(
          nonlinearMechanicsBVPoperator->getVolumeOperator() ) )
        ->setVector( AMP::Operator::Mechanics::TEMPERATURE, curTempVec );

    // Create a Linear BVP operator for mechanics
    AMP_INSIST( input_db->keyExists( "LinearMechanicsOperator" ), "key missing!" );
    AMP::shared_ptr<AMP::Operator::LinearBVPOperator> linearMechanicsBVPoperator =
        AMP::dynamic_pointer_cast<AMP::Operator::LinearBVPOperator>(
            AMP::Operator::OperatorBuilder::createOperator(
                meshAdapter, "LinearMechanicsOperator", input_db, mechanicsMaterialModel ) );

    // We need to reset the linear operator before the solve since TrilinosML does
    // the factorization of the matrix during construction and so the matrix must
    // be correct before constructing the TrilinosML object.
    nonlinearMechanicsBVPoperator->apply( solVec, resVec );
    linearMechanicsBVPoperator->reset(
        nonlinearMechanicsBVPoperator->getParameters( "Jacobian", solVec ) );

    double epsilon =
        1.0e-13 * ( ( ( linearMechanicsBVPoperator->getMatrix() )->extractDiagonal() )->L1Norm() );

    AMP::shared_ptr<AMP::Database> nonlinearSolver_db = input_db->getDatabase( "NonlinearSolver" );
    AMP::shared_ptr<AMP::Database> linearSolver_db =
        nonlinearSolver_db->getDatabase( "LinearSolver" );

    // ---- first initialize the preconditioner
    AMP::shared_ptr<AMP::Database> pcSolver_db = linearSolver_db->getDatabase( "Preconditioner" );
    AMP::shared_ptr<AMP::Solver::TrilinosMLSolverParameters> pcSolverParams(
        new AMP::Solver::TrilinosMLSolverParameters( pcSolver_db ) );
    pcSolverParams->d_pOperator = linearMechanicsBVPoperator;
    AMP::shared_ptr<AMP::Solver::TrilinosMLSolver> pcSolver(
        new AMP::Solver::TrilinosMLSolver( pcSolverParams ) );

    // HACK to prevent a double delete on Petsc Vec
    AMP::shared_ptr<AMP::Solver::PetscSNESSolver> nonlinearSolver;

    // initialize the linear solver
    AMP::shared_ptr<AMP::Solver::PetscKrylovSolverParameters> linearSolverParams(
        new AMP::Solver::PetscKrylovSolverParameters( linearSolver_db ) );
    linearSolverParams->d_pOperator       = linearMechanicsBVPoperator;
    linearSolverParams->d_comm            = globalComm;
    linearSolverParams->d_pPreconditioner = pcSolver;
    AMP::shared_ptr<AMP::Solver::PetscKrylovSolver> linearSolver(
        new AMP::Solver::PetscKrylovSolver( linearSolverParams ) );

    // initialize the nonlinear solver
    AMP::shared_ptr<AMP::Solver::PetscSNESSolverParameters> nonlinearSolverParams(
        new AMP::Solver::PetscSNESSolverParameters( nonlinearSolver_db ) );
    // change the next line to get the correct communicator out
    nonlinearSolverParams->d_comm          = globalComm;
    nonlinearSolverParams->d_pOperator     = nonlinearMechanicsBVPoperator;
    nonlinearSolverParams->d_pKrylovSolver = linearSolver;
    nonlinearSolverParams->d_pInitialGuess = solVec;
    nonlinearSolver.reset( new AMP::Solver::PetscSNESSolver( nonlinearSolverParams ) );

    nonlinearSolver->setZeroInitialGuess( false );

    int loadingSubSteps = 1;
    double scaleValue;

    for ( int step = 0; step < NumberOfLoadingSteps; step++ ) {
        AMP::pout << "########################################" << std::endl;
        AMP::pout << "The current loading step is " << ( step + 1 ) << std::endl;

        if ( step > 3 ) {
            loadingSubSteps = 1;
        }
        for ( int subStep = 0; subStep < loadingSubSteps; subStep++ ) {
            if ( step <= 3 ) {
                scaleValue = ( (double) step + 1.0 ) / NumberOfLoadingSteps;
            } else {
                scaleValue =
                    ( ( (double) step ) / NumberOfLoadingSteps ) +
                    ( ( (double) subStep + 1.0 ) / ( NumberOfLoadingSteps * loadingSubSteps ) );
            }

            AMP::pout << "########################################" << std::endl;
            AMP::pout << "The current loading sub step is " << ( subStep + 1 )
                      << " and scaleValue = " << scaleValue << std::endl;

            scaledRhsVec->scale( scaleValue, rhsVec );
            AMP::pout << "L2 Norm of RHS at loading step " << ( step + 1 ) << " is "
                      << scaledRhsVec->L2Norm() << std::endl;

            nonlinearMechanicsBVPoperator->residual( scaledRhsVec, solVec, resVec );
            double initialResidualNorm = resVec->L2Norm();
            AMP::pout << "Initial Residual Norm for loading step " << ( step + 1 ) << " is "
                      << initialResidualNorm << std::endl;

            nonlinearSolver->solve( scaledRhsVec, solVec );

            nonlinearMechanicsBVPoperator->residual( scaledRhsVec, solVec, resVec );
            double finalResidualNorm = resVec->L2Norm();
            AMP::pout << "Final Residual Norm for loading step " << ( step + 1 ) << " is "
                      << finalResidualNorm << std::endl;

            if ( finalResidualNorm > ( 1.0e-10 * initialResidualNorm ) ) {
                ut->failure( "Nonlinear solve for current loading step" );
            } else {
                ut->passes( "Nonlinear solve for current loading step" );
            }

            double finalSolNorm = solVec->L2Norm();

            AMP::pout << "Final Solution Norm: " << finalSolNorm << std::endl;

            AMP::LinearAlgebra::Vector::shared_ptr mechUvec =
                solVec->select( AMP::LinearAlgebra::VS_Stride( 0, 3 ), "U" );
            AMP::LinearAlgebra::Vector::shared_ptr mechVvec =
                solVec->select( AMP::LinearAlgebra::VS_Stride( 1, 3 ), "V" );
            AMP::LinearAlgebra::Vector::shared_ptr mechWvec =
                solVec->select( AMP::LinearAlgebra::VS_Stride( 2, 3 ), "W" );

            double finalMaxU = mechUvec->maxNorm();
            double finalMaxV = mechVvec->maxNorm();
            double finalMaxW = mechWvec->maxNorm();

            AMP::pout << "Maximum U displacement: " << finalMaxU << std::endl;
            AMP::pout << "Maximum V displacement: " << finalMaxV << std::endl;
            AMP::pout << "Maximum W displacement: " << finalMaxW << std::endl;

            AMP::shared_ptr<AMP::InputDatabase> tmp_db( new AMP::InputDatabase( "Dummy" ) );
            AMP::shared_ptr<AMP::Operator::MechanicsNonlinearFEOperatorParameters> tmpParams(
                new AMP::Operator::MechanicsNonlinearFEOperatorParameters( tmp_db ) );
            ( nonlinearMechanicsBVPoperator->getVolumeOperator() )->reset( tmpParams );
            nonlinearSolver->setZeroInitialGuess( false );
        }
    }

    AMP::pout << "epsilon = " << epsilon << std::endl;

    ut->passes( exeName );
}

int main( int argc, char *argv[] )
{
    AMP::AMPManager::startup( argc, argv );
    AMP::UnitTest ut;

    std::vector<std::string> exeNames;
    exeNames.push_back( "testUpdatedLagrangianThermoMechanics-LinearElasticity-1" );

    for ( auto &exeName : exeNames )
        myTest( &ut, exeName );

    ut.report();

    int num_failed = ut.NumFailGlobal();
    AMP::AMPManager::shutdown();
    return num_failed;
}
