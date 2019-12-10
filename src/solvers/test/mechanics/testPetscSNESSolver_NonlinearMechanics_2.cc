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
#include "AMP/operators/mechanics/ThermalVonMisesMatModel.h"
#include "AMP/solvers/petsc/PetscKrylovSolver.h"
#include "AMP/solvers/petsc/PetscKrylovSolverParameters.h"
#include "AMP/solvers/petsc/PetscSNESSolver.h"
#include "AMP/solvers/petsc/PetscSNESSolverParameters.h"
#include "AMP/solvers/trilinos/ml/TrilinosMLSolver.h"
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


static void myTest( AMP::UnitTest *ut )
{
    std::string exeName( "testPetscSNESSolver-NonlinearMechanics-2" );
    std::string input_file = "input_" + exeName;
    std::string log_file   = "output_" + exeName;

    AMP::PIO::logOnlyNodeZero( log_file );
    AMP::AMP_MPI globalComm( AMP_COMM_WORLD );

#ifdef USE_EXT_SILO
    // Create the silo writer and register the data
    auto siloWriter = AMP::Utilities::Writer::buildWriter( "Silo" );
#endif

    auto input_db = AMP::Database::parseInputFile( input_file );
    input_db->print( AMP::plog );

    AMP_INSIST( input_db->keyExists( "Mesh" ), "Key ''Mesh'' is missing!" );
    auto mesh_db    = input_db->getDatabase( "Mesh" );
    auto meshParams = AMP::make_shared<AMP::Mesh::MeshParameters>( mesh_db );
    meshParams->setComm( AMP::AMP_MPI( AMP_COMM_WORLD ) );
    auto meshAdapter = AMP::Mesh::Mesh::buildMesh( meshParams );

    AMP_INSIST( input_db->keyExists( "NumberOfLoadingSteps" ),
                "Key ''NumberOfLoadingSteps'' is missing!" );
    int NumberOfLoadingSteps = input_db->getScalar<int>( "NumberOfLoadingSteps" );

    auto nonlinBvpOperator = AMP::dynamic_pointer_cast<AMP::Operator::NonlinearBVPOperator>(
        AMP::Operator::OperatorBuilder::createOperator(
            meshAdapter, "nonlinearMechanicsBVPOperator", input_db ) );
    auto nonlinearMechanicsVolumeOperator =
        AMP::dynamic_pointer_cast<AMP::Operator::MechanicsNonlinearFEOperator>(
            nonlinBvpOperator->getVolumeOperator() );
    auto elementPhysicsModel = nonlinearMechanicsVolumeOperator->getMaterialModel();

    auto temperatureVariable  = AMP::make_shared<AMP::LinearAlgebra::Variable>( "temp" );
    auto displacementVariable = nonlinBvpOperator->getOutputVariable();

    auto tempDofMap = AMP::Discretization::simpleDOFManager::create(
        meshAdapter, AMP::Mesh::GeomType::Vertex, 1, 1, true );

    auto initTempVec  = AMP::LinearAlgebra::createVector( tempDofMap, temperatureVariable, true );
    auto finalTempVec = initTempVec->cloneVector();

    initTempVec->setRandomValues();
    initTempVec->abs( initTempVec );
    double initTempConst = input_db->getWithDefault<double>( "INIT_TEMP_CONST", 10.0 );
    initTempVec->scale( initTempConst );
    initTempVec->makeConsistent( AMP::LinearAlgebra::Vector::ScatterType::CONSISTENT_SET );

    bool setFinalTempEqualsInitialTemp =
        input_db->getWithDefault( "SET_FINAL_TEMP_EQUALS_INIT_TEMP", false );

    if ( setFinalTempEqualsInitialTemp ) {
        finalTempVec->copyVector( initTempVec );
    } else {
        finalTempVec->setRandomValues();
        double finalTempConst = input_db->getWithDefault<double>( "FINAL_TEMP_CONST", 12.0 );
        finalTempVec->scale( finalTempConst );
    }
    finalTempVec->makeConsistent( AMP::LinearAlgebra::Vector::ScatterType::CONSISTENT_SET );

    ( AMP::dynamic_pointer_cast<AMP::Operator::MechanicsNonlinearFEOperator>(
          nonlinBvpOperator->getVolumeOperator() ) )
        ->setReferenceTemperature( initTempVec );
    ( AMP::dynamic_pointer_cast<AMP::Operator::MechanicsNonlinearFEOperator>(
          nonlinBvpOperator->getVolumeOperator() ) )
        ->setVector( AMP::Operator::Mechanics::TEMPERATURE, finalTempVec );

    auto linBvpOperator = AMP::dynamic_pointer_cast<AMP::Operator::LinearBVPOperator>(
        AMP::Operator::OperatorBuilder::createOperator(
            meshAdapter, "linearMechanicsBVPOperator", input_db, elementPhysicsModel ) );

    // For RHS (Point Forces)
    AMP::shared_ptr<AMP::Operator::ElementPhysicsModel> dummyModel;
    auto dirichletLoadVecOp = AMP::dynamic_pointer_cast<AMP::Operator::DirichletVectorCorrection>(
        AMP::Operator::OperatorBuilder::createOperator(
            meshAdapter, "Load_Boundary", input_db, dummyModel ) );
    dirichletLoadVecOp->setVariable( displacementVariable );

    // For Initial-Guess
    auto dirichletDispInVecOp = AMP::dynamic_pointer_cast<AMP::Operator::DirichletVectorCorrection>(
        AMP::Operator::OperatorBuilder::createOperator(
            meshAdapter, "Displacement_Boundary", input_db, dummyModel ) );
    dirichletDispInVecOp->setVariable( displacementVariable );

    AMP::Discretization::DOFManager::shared_ptr dispDofMap =
        AMP::Discretization::simpleDOFManager::create(
            meshAdapter, AMP::Mesh::GeomType::Vertex, 1, 3, true );

    AMP::LinearAlgebra::Vector::shared_ptr nullVec;
    auto mechNlSolVec = AMP::LinearAlgebra::createVector( dispDofMap, displacementVariable, true );
    auto mechNlRhsVec = mechNlSolVec->cloneVector();
    auto mechNlResVec = mechNlSolVec->cloneVector();
    auto mechNlScaledRhsVec = mechNlSolVec->cloneVector();


    // Initial guess for NL solver must satisfy the displacement boundary conditions
    mechNlSolVec->setToScalar( 0.0 );
    dirichletDispInVecOp->apply( nullVec, mechNlSolVec );
    mechNlSolVec->makeConsistent( AMP::LinearAlgebra::Vector::ScatterType::CONSISTENT_SET );

    nonlinBvpOperator->apply( mechNlSolVec, mechNlResVec );
    linBvpOperator->reset( nonlinBvpOperator->getParameters( "Jacobian", mechNlSolVec ) );

    mechNlRhsVec->setToScalar( 0.0 );
    dirichletLoadVecOp->apply( nullVec, mechNlRhsVec );

    double initSolNorm = mechNlSolVec->L2Norm();

    std::cout << "Initial Solution Norm: " << initSolNorm << std::endl;

    auto nonlinearSolver_db = input_db->getDatabase( "NonlinearSolver" );
    auto linearSolver_db    = nonlinearSolver_db->getDatabase( "LinearSolver" );

    // ---- first initialize the preconditioner
    auto pcSolver_db    = linearSolver_db->getDatabase( "Preconditioner" );
    auto pcSolverParams = AMP::make_shared<AMP::Solver::SolverStrategyParameters>( pcSolver_db );
    pcSolverParams->d_pOperator = linBvpOperator;
    auto pcSolver               = AMP::make_shared<AMP::Solver::TrilinosMLSolver>( pcSolverParams );

    // HACK to prevent a double delete on Petsc Vec
    AMP::shared_ptr<AMP::Solver::PetscSNESSolver> nonlinearSolver;

    // initialize the linear solver
    auto linearSolverParams =
        AMP::make_shared<AMP::Solver::PetscKrylovSolverParameters>( linearSolver_db );
    linearSolverParams->d_pOperator       = linBvpOperator;
    linearSolverParams->d_comm            = globalComm;
    linearSolverParams->d_pPreconditioner = pcSolver;
    auto linearSolver = AMP::make_shared<AMP::Solver::PetscKrylovSolver>( linearSolverParams );
    auto nonlinearSolverParams =
        AMP::make_shared<AMP::Solver::PetscSNESSolverParameters>( nonlinearSolver_db );

    // change the next line to get the correct communicator out
    nonlinearSolverParams->d_comm          = globalComm;
    nonlinearSolverParams->d_pOperator     = nonlinBvpOperator;
    nonlinearSolverParams->d_pKrylovSolver = linearSolver;
    nonlinearSolverParams->d_pInitialGuess = mechNlSolVec;

    nonlinearSolver.reset( new AMP::Solver::PetscSNESSolver( nonlinearSolverParams ) );

    nonlinearSolver->setZeroInitialGuess( false );

#ifdef USE_EXT_SILO
    siloWriter->registerVector(
        mechNlSolVec, meshAdapter, AMP::Mesh::GeomType::Vertex, "Solution" );
#endif

    for ( int step = 0; step < NumberOfLoadingSteps; step++ ) {
        AMP::pout << "########################################" << std::endl;
        AMP::pout << "The current loading step is " << ( step + 1 ) << std::endl;

        double scaleValue = ( (double) step + 1.0 ) / NumberOfLoadingSteps;
        mechNlScaledRhsVec->scale( scaleValue, mechNlRhsVec );
        mechNlScaledRhsVec->makeConsistent(
            AMP::LinearAlgebra::Vector::ScatterType::CONSISTENT_SET );
        AMP::pout << "L2 Norm at loading step " << ( step + 1 ) << " is "
                  << mechNlScaledRhsVec->L2Norm() << std::endl;

        nonlinBvpOperator->residual( mechNlScaledRhsVec, mechNlSolVec, mechNlResVec );
        double initialResidualNorm = mechNlResVec->L2Norm();
        AMP::pout << "Initial Residual Norm for loading step " << ( step + 1 ) << " is "
                  << initialResidualNorm << std::endl;

        if ( initialResidualNorm < 1.0e-2 ) {
            ut->passes( "Nonlinear solve for current loading step" );
        } else {
            AMP::pout << "Starting Nonlinear Solve..." << std::endl;
            nonlinearSolver->solve( mechNlScaledRhsVec, mechNlSolVec );

            nonlinBvpOperator->residual( mechNlScaledRhsVec, mechNlSolVec, mechNlResVec );
            double finalResidualNorm = mechNlResVec->L2Norm();
            AMP::pout << "Final Residual Norm for loading step " << ( step + 1 ) << " is "
                      << finalResidualNorm << std::endl;

            if ( finalResidualNorm > ( 1.0e-8 * initialResidualNorm ) ) {
                ut->failure( "Nonlinear solve for current loading step" );
            } else {
                ut->passes( "Nonlinear solve for current loading step" );
            }
        }

        auto tmp_db = AMP::make_shared<AMP::Database>( "Dummy" );
        auto tmpParams =
            AMP::make_shared<AMP::Operator::MechanicsNonlinearFEOperatorParameters>( tmp_db );
        ( nonlinBvpOperator->getVolumeOperator() )->reset( tmpParams );
        nonlinearSolver->setZeroInitialGuess( false );

#ifdef USE_EXT_SILO
        siloWriter->writeFile( "FrozenTemp_NonlinearMechExample", step );
#endif
    }

    double finalSolNorm = mechNlSolVec->L2Norm();
    AMP::pout << "Final Solution Norm: " << finalSolNorm << std::endl;

    ut->passes( exeName );
}


int testPetscSNESSolver_NonlinearMechanics_2( int argc, char *argv[] )
{
    AMP::AMPManager::startup( argc, argv );
    AMP::UnitTest ut;

    myTest( &ut );

    ut.report();

    int num_failed = ut.NumFailGlobal();
    AMP::AMPManager::shutdown();
    return num_failed;
}
