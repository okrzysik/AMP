#include "AMP/ampmesh/Mesh.h"
#include "AMP/discretization/DOF_Manager.h"
#include "AMP/discretization/simpleDOF_Manager.h"
#include "AMP/materials/Material.h"
#include "AMP/operators/BVPOperatorParameters.h"
#include "AMP/operators/LinearBVPOperator.h"
#include "AMP/operators/NonlinearBVPOperator.h"
#include "AMP/operators/OperatorBuilder.h"
#include "AMP/operators/boundary/DirichletVectorCorrection.h"
#include "AMP/operators/boundary/libmesh/PressureBoundaryOperator.h"
#include "AMP/operators/mechanics/MechanicsLinearFEOperator.h"
#include "AMP/operators/mechanics/MechanicsNonlinearFEOperator.h"
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
#include <memory>

#include <iostream>
#include <string>


static void myTest( AMP::UnitTest *ut, const std::string &exeName )
{
    std::string input_file = "input_" + exeName;
    std::string log_file   = "output_" + exeName;

    AMP::PIO::logOnlyNodeZero( log_file );
    AMP::AMP_MPI globalComm( AMP_COMM_WORLD );

#ifdef USE_EXT_SILO
    // Create the silo writer and register the data
    AMP::Utilities::Writer::shared_ptr siloWriter = AMP::Utilities::Writer::buildWriter( "Silo" );
#endif


    auto input_db = AMP::Database::parseInputFile( input_file );
    input_db->print( AMP::plog );

    //--------------------------------------------------
    //   Create the Mesh.
    //--------------------------------------------------
    AMP_INSIST( input_db->keyExists( "Mesh" ), "Key ''Mesh'' is missing!" );
    std::shared_ptr<AMP::Database> mesh_db = input_db->getDatabase( "Mesh" );
    std::shared_ptr<AMP::Mesh::MeshParameters> meshParams(
        new AMP::Mesh::MeshParameters( mesh_db ) );
    meshParams->setComm( AMP::AMP_MPI( AMP_COMM_WORLD ) );
    AMP::Mesh::Mesh::shared_ptr meshAdapter = AMP::Mesh::Mesh::buildMesh( meshParams );

    AMP_INSIST( input_db->keyExists( "NumberOfLoadingSteps" ),
                "Key ''NumberOfLoadingSteps'' is missing!" );
    int NumberOfLoadingSteps = input_db->getScalar<int>( "NumberOfLoadingSteps" );

    std::shared_ptr<AMP::Operator::NonlinearBVPOperator> nonlinBvpOperator =
        std::dynamic_pointer_cast<AMP::Operator::NonlinearBVPOperator>(
            AMP::Operator::OperatorBuilder::createOperator(
                meshAdapter, "nonlinearMechanicsBVPOperator", input_db ) );
    std::shared_ptr<AMP::Operator::MechanicsNonlinearFEOperator> nonlinearMechanicsVolumeOperator =
        std::dynamic_pointer_cast<AMP::Operator::MechanicsNonlinearFEOperator>(
            nonlinBvpOperator->getVolumeOperator() );
    std::shared_ptr<AMP::Operator::ElementPhysicsModel> elementPhysicsModel =
        nonlinearMechanicsVolumeOperator->getMaterialModel();

    std::shared_ptr<AMP::Operator::LinearBVPOperator> linBvpOperator =
        std::dynamic_pointer_cast<AMP::Operator::LinearBVPOperator>(
            AMP::Operator::OperatorBuilder::createOperator(
                meshAdapter, "linearMechanicsBVPOperator", input_db, elementPhysicsModel ) );

    AMP::LinearAlgebra::Variable::shared_ptr displacementVariable =
        std::dynamic_pointer_cast<AMP::Operator::MechanicsNonlinearFEOperator>(
            nonlinBvpOperator->getVolumeOperator() )
            ->getOutputVariable();

    // For RHS (Point Forces)
    std::shared_ptr<AMP::Operator::ElementPhysicsModel> dummyModel;
    std::shared_ptr<AMP::Operator::DirichletVectorCorrection> dirichletLoadVecOp =
        std::dynamic_pointer_cast<AMP::Operator::DirichletVectorCorrection>(
            AMP::Operator::OperatorBuilder::createOperator(
                meshAdapter, "Load_Boundary", input_db, dummyModel ) );
    dirichletLoadVecOp->setVariable( displacementVariable );

    // Pressure RHS
    std::shared_ptr<AMP::Operator::PressureBoundaryOperator> pressureLoadVecOp =
        std::dynamic_pointer_cast<AMP::Operator::PressureBoundaryOperator>(
            AMP::Operator::OperatorBuilder::createOperator(
                meshAdapter, "Pressure_Boundary", input_db, dummyModel ) );

    // For Initial-Guess
    std::shared_ptr<AMP::Operator::DirichletVectorCorrection> dirichletDispInVecOp =
        std::dynamic_pointer_cast<AMP::Operator::DirichletVectorCorrection>(
            AMP::Operator::OperatorBuilder::createOperator(
                meshAdapter, "Displacement_Boundary", input_db, dummyModel ) );
    dirichletDispInVecOp->setVariable( displacementVariable );

    AMP::Discretization::DOFManager::shared_ptr nodalDofMap =
        AMP::Discretization::simpleDOFManager::create(
            meshAdapter, AMP::Mesh::GeomType::Vertex, 1, 3, true );

    AMP::LinearAlgebra::Vector::shared_ptr nullVec;

    AMP::LinearAlgebra::Vector::shared_ptr mechNlSolVec =
        AMP::LinearAlgebra::createVector( nodalDofMap, displacementVariable, true );
    AMP::LinearAlgebra::Vector::shared_ptr mechNlRhsVec       = mechNlSolVec->cloneVector();
    AMP::LinearAlgebra::Vector::shared_ptr mechNlResVec       = mechNlSolVec->cloneVector();
    AMP::LinearAlgebra::Vector::shared_ptr mechNlScaledRhsVec = mechNlSolVec->cloneVector();
    AMP::LinearAlgebra::Vector::shared_ptr mechPressureVec    = mechNlSolVec->cloneVector();

    // Initial guess for NL solver must satisfy the displacement boundary conditions
    mechNlSolVec->setToScalar( 0.0, mechNlSolVec );
    mechPressureVec->setToScalar( 0.0, mechPressureVec );

    dirichletDispInVecOp->apply( nullVec, mechNlSolVec );
    mechNlSolVec->makeConsistent( AMP::LinearAlgebra::Vector::ScatterType::CONSISTENT_SET );

    nonlinBvpOperator->apply( mechNlSolVec, mechNlResVec );
    linBvpOperator->reset( nonlinBvpOperator->getParameters( "Jacobian", mechNlSolVec ) );

    // Point forces
    mechNlRhsVec->setToScalar( 0.0, mechNlRhsVec );
    dirichletLoadVecOp->apply( nullVec, mechNlRhsVec );

    // Applying the pressure load
    pressureLoadVecOp->addRHScorrection( mechPressureVec );
    mechNlRhsVec->add( mechNlRhsVec, mechPressureVec, mechNlRhsVec );

    std::shared_ptr<AMP::Database> nonlinearSolver_db = input_db->getDatabase( "NonlinearSolver" );
    std::shared_ptr<AMP::Database> linearSolver_db =
        nonlinearSolver_db->getDatabase( "LinearSolver" );

    // ---- first initialize the preconditioner
    std::shared_ptr<AMP::Database> pcSolver_db = linearSolver_db->getDatabase( "Preconditioner" );
    std::shared_ptr<AMP::Solver::TrilinosMLSolverParameters> pcSolverParams(
        new AMP::Solver::TrilinosMLSolverParameters( pcSolver_db ) );
    pcSolverParams->d_pOperator = linBvpOperator;
    std::shared_ptr<AMP::Solver::TrilinosMLSolver> pcSolver(
        new AMP::Solver::TrilinosMLSolver( pcSolverParams ) );

    // HACK to prevent a double delete on Petsc Vec
    std::shared_ptr<AMP::Solver::PetscSNESSolver> nonlinearSolver;

    // initialize the linear solver
    std::shared_ptr<AMP::Solver::PetscKrylovSolverParameters> linearSolverParams(
        new AMP::Solver::PetscKrylovSolverParameters( linearSolver_db ) );
    linearSolverParams->d_pOperator       = linBvpOperator;
    linearSolverParams->d_comm            = globalComm;
    linearSolverParams->d_pPreconditioner = pcSolver;
    std::shared_ptr<AMP::Solver::PetscKrylovSolver> linearSolver(
        new AMP::Solver::PetscKrylovSolver( linearSolverParams ) );

    // initialize the nonlinear solver
    std::shared_ptr<AMP::Solver::PetscSNESSolverParameters> nonlinearSolverParams(
        new AMP::Solver::PetscSNESSolverParameters( nonlinearSolver_db ) );
    // change the next line to get the correct communicator out
    nonlinearSolverParams->d_comm          = globalComm;
    nonlinearSolverParams->d_pOperator     = nonlinBvpOperator;
    nonlinearSolverParams->d_pKrylovSolver = linearSolver;
    nonlinearSolverParams->d_pInitialGuess = mechNlSolVec;
    nonlinearSolver.reset( new AMP::Solver::PetscSNESSolver( nonlinearSolverParams ) );

    nonlinearSolver->setZeroInitialGuess( false );

    for ( int step = 0; step < NumberOfLoadingSteps; step++ ) {
        AMP::pout << "########################################" << std::endl;
        AMP::pout << "The current loading step is " << ( step + 1 ) << std::endl;

        double scaleValue = ( (double) step + 1.0 ) / NumberOfLoadingSteps;
        mechNlScaledRhsVec->scale( scaleValue, mechNlRhsVec, mechNlScaledRhsVec );
        mechNlScaledRhsVec->makeConsistent(
            AMP::LinearAlgebra::Vector::ScatterType::CONSISTENT_SET );
        AMP::pout << "L2 Norm of RHS at loading step " << ( step + 1 ) << " is "
                  << mechNlScaledRhsVec->L2Norm(mechNlScaledRhsVec) << std::endl;

        nonlinBvpOperator->residual( mechNlScaledRhsVec, mechNlSolVec, mechNlResVec );
        double initialResidualNorm = mechNlResVec->L2Norm(mechNlResVec);
        AMP::pout << "Initial Residual Norm for loading step " << ( step + 1 ) << " is "
                  << initialResidualNorm << std::endl;

        AMP::pout << "Starting Nonlinear Solve..." << std::endl;
        nonlinearSolver->solve( mechNlScaledRhsVec, mechNlSolVec );

        nonlinBvpOperator->residual( mechNlScaledRhsVec, mechNlSolVec, mechNlResVec );
        double finalResidualNorm = mechNlResVec->L2Norm(mechNlResVec);
        AMP::pout << "Final Residual Norm for loading step " << ( step + 1 ) << " is "
                  << finalResidualNorm << std::endl;

        if ( finalResidualNorm > ( 1.0e-10 * initialResidualNorm ) ) {
            ut->failure( "Nonlinear solve for current loading step" );
        } else {
            ut->passes( "Nonlinear solve for current loading step" );
        }

        double finalSolNorm = mechNlSolVec->L2Norm(mechNlSolVec);

        AMP::pout << "Final Solution Norm: " << finalSolNorm << std::endl;

        AMP::LinearAlgebra::Vector::shared_ptr mechUvec =
            mechNlSolVec->select( AMP::LinearAlgebra::VS_Stride( 0, 3 ), "U" );
        AMP::LinearAlgebra::Vector::shared_ptr mechVvec =
            mechNlSolVec->select( AMP::LinearAlgebra::VS_Stride( 1, 3 ), "V" );
        AMP::LinearAlgebra::Vector::shared_ptr mechWvec =
            mechNlSolVec->select( AMP::LinearAlgebra::VS_Stride( 2, 3 ), "W" );

        double finalMaxU = mechUvec->maxNorm(mechUvec);
        double finalMaxV = mechVvec->maxNorm(mechVvec);
        double finalMaxW = mechWvec->maxNorm(mechWvec);

        AMP::pout << "Maximum U displacement: " << finalMaxU << std::endl;
        AMP::pout << "Maximum V displacement: " << finalMaxV << std::endl;
        AMP::pout << "Maximum W displacement: " << finalMaxW << std::endl;

        std::shared_ptr<AMP::Database> tmp_db( new AMP::Database( "Dummy" ) );
        std::shared_ptr<AMP::Operator::MechanicsNonlinearFEOperatorParameters> tmpParams(
            new AMP::Operator::MechanicsNonlinearFEOperatorParameters( tmp_db ) );
        ( nonlinBvpOperator->getVolumeOperator() )->reset( tmpParams );
        nonlinearSolver->setZeroInitialGuess( false );

#ifdef USE_EXT_SILO
        siloWriter->registerVector(
            mechNlSolVec, meshAdapter, AMP::Mesh::GeomType::Vertex, "Solution_Vector" );
        meshAdapter->displaceMesh( mechNlSolVec );
        char outFileName2[256];
        sprintf( outFileName2, "PressurePrescribed-DeformedBrick-LinearElasticity_%d", step );
        siloWriter->writeFile( outFileName2, 1 );
#endif
    }

    double finalSolNorm = mechNlSolVec->L2Norm(mechNlSolVec);
    AMP::pout << "Final Solution Norm: " << finalSolNorm << std::endl;

#ifdef USE_EXT_SILO
    siloWriter->writeFile( exeName, 1 );
#endif

    ut->passes( exeName );
}

int testBrickPressureBoundary_1( int argc, char *argv[] )
{
    AMP::AMPManager::startup( argc, argv );
    AMP::UnitTest ut;

    std::vector<std::string> exeNames;
    exeNames.emplace_back( "testBrickPressureBoundary-1" );

    for ( auto &exeName : exeNames )
        myTest( &ut, exeName );

    ut.report();

    int num_failed = ut.NumFailGlobal();
    AMP::AMPManager::shutdown();
    return num_failed;
}
