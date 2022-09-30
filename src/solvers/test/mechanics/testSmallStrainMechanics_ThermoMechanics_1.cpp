#include "AMP/IO/PIO.h"
#include "AMP/IO/Writer.h"
#include "AMP/discretization/simpleDOF_Manager.h"
#include "AMP/mesh/MeshParameters.h"
#include "AMP/mesh/libmesh/ReadTestMesh.h"
#include "AMP/operators/BVPOperatorParameters.h"
#include "AMP/operators/LinearBVPOperator.h"
#include "AMP/operators/NonlinearBVPOperator.h"
#include "AMP/operators/OperatorBuilder.h"
#include "AMP/operators/boundary/DirichletMatrixCorrection.h"
#include "AMP/operators/boundary/DirichletVectorCorrection.h"
#include "AMP/operators/mechanics/MechanicsLinearFEOperator.h"
#include "AMP/operators/mechanics/MechanicsNonlinearFEOperator.h"
#include "AMP/operators/mechanics/ThermalStrainMaterialModel.h"
#include "AMP/solvers/SolverStrategyParameters.h"
#include "AMP/solvers/petsc/PetscKrylovSolver.h"
#include "AMP/solvers/SolverStrategyParameters.h"
#include "AMP/solvers/petsc/PetscSNESSolver.h"
#include "AMP/solvers/trilinos/ml/TrilinosMLSolver.h"
#include "AMP/utils/AMPManager.h"
#include "AMP/utils/AMP_MPI.h"
#include "AMP/utils/Database.h"
#include "AMP/utils/UnitTest.h"
#include "AMP/utils/Utilities.h"
#include "AMP/vectors/Vector.h"
#include "AMP/vectors/VectorBuilder.h"
#include "AMP/vectors/VectorSelector.h"

#include <iostream>
#include <memory>
#include <string>


static void myTest( AMP::UnitTest *ut, const std::string &exeName )
{
    std::string input_file = "input_" + exeName;
    std::string log_file   = "log_" + exeName;

    AMP::logOnlyNodeZero( log_file );
    AMP::AMP_MPI globalComm( AMP_COMM_WORLD );

    // Read the input file
    auto input_db = AMP::Database::parseInputFile( input_file );
    input_db->print( AMP::plog );

    AMP_INSIST( input_db->keyExists( "Mesh" ), "Key ''Mesh'' is missing!" );
    auto mesh_db    = input_db->getDatabase( "Mesh" );
    auto meshParams = std::make_shared<AMP::Mesh::MeshParameters>( mesh_db );
    meshParams->setComm( AMP::AMP_MPI( AMP_COMM_WORLD ) );
    auto meshAdapter = AMP::Mesh::MeshFactory::create( meshParams );

    AMP_INSIST( input_db->keyExists( "NumberOfLoadingSteps" ),
                "Key ''NumberOfLoadingSteps'' is missing!" );
    int NumberOfLoadingSteps = input_db->getScalar<int>( "NumberOfLoadingSteps" );
    AMP::pout << "NumberOfLoadingSteps = " << NumberOfLoadingSteps << std::endl;

    // Create a nonlinear BVP operator for mechanics
    AMP_INSIST( input_db->keyExists( "NonlinearMechanicsOperator" ), "key missing!" );
    std::shared_ptr<AMP::Operator::ElementPhysicsModel> mechanicsMaterialModel;
    auto nonlinearMechanicsBVPoperator =
        std::dynamic_pointer_cast<AMP::Operator::NonlinearBVPOperator>(
            AMP::Operator::OperatorBuilder::createOperator(
                meshAdapter, "NonlinearMechanicsOperator", input_db, mechanicsMaterialModel ) );

    // Create the variables
    auto mechanicsNonlinearVolumeOperator =
        std::dynamic_pointer_cast<AMP::Operator::MechanicsNonlinearFEOperator>(
            nonlinearMechanicsBVPoperator->getVolumeOperator() );
    auto dispVar = nonlinearMechanicsBVPoperator->getOutputVariable();
    auto tempVar = std::make_shared<AMP::LinearAlgebra::Variable>( "Temperature" );
    auto burnVar = std::make_shared<AMP::LinearAlgebra::Variable>( "Burnup" );
    auto lhgrVar = std::make_shared<AMP::LinearAlgebra::Variable>( "LHGR" );

    // For RHS (Point Forces)
    std::shared_ptr<AMP::Operator::ElementPhysicsModel> dummyModel;
    auto dirichletLoadVecOp = std::dynamic_pointer_cast<AMP::Operator::DirichletVectorCorrection>(
        AMP::Operator::OperatorBuilder::createOperator(
            meshAdapter, "Load_Boundary", input_db, dummyModel ) );
    dirichletLoadVecOp->setVariable( dispVar );

    auto vectorDofMap = AMP::Discretization::simpleDOFManager::create(
        meshAdapter, AMP::Mesh::GeomType::Vertex, 1, 3, true );

    auto scalarDofMap = AMP::Discretization::simpleDOFManager::create(
        meshAdapter, AMP::Mesh::GeomType::Vertex, 1, 1, true );

    // Create the vectors
    AMP::LinearAlgebra::Vector::shared_ptr nullVec;
    auto solVec       = AMP::LinearAlgebra::createVector( vectorDofMap, dispVar, true );
    auto tempVec      = AMP::LinearAlgebra::createVector( scalarDofMap, tempVar, true );
    auto burnVec      = AMP::LinearAlgebra::createVector( scalarDofMap, burnVar, true );
    auto lhgrVec      = AMP::LinearAlgebra::createVector( scalarDofMap, lhgrVar, true );
    auto rhsVec       = solVec->cloneVector();
    auto resVec       = solVec->cloneVector();
    auto scaledRhsVec = solVec->cloneVector();
    auto tempVecRef   = tempVec->cloneVector();

    // Initial guess
    solVec->zero();
    nonlinearMechanicsBVPoperator->modifyInitialSolutionVector( solVec );

    // RHS
    rhsVec->zero();
    dirichletLoadVecOp->apply( nullVec, rhsVec );
    nonlinearMechanicsBVPoperator->modifyRHSvector( rhsVec );

    tempVecRef->setToScalar( 301.0 );
    tempVec->setToScalar( 301.0 );

    burnVec->setRandomValues();
    burnVec->abs( *burnVec );
    double maxBurnVec        = static_cast<double>( burnVec->max() );
    double oneOverMaxBurnVec = 1.0 / maxBurnVec;
    burnVec->scale( oneOverMaxBurnVec );
    burnVec->scale( 0.2 );

    lhgrVec->setRandomValues();
    lhgrVec->abs( *lhgrVec );
    double maxLHGRVec        = static_cast<double>( lhgrVec->max() );
    double oneOverMaxLHGRVec = 1.0 / maxLHGRVec;
    lhgrVec->scale( oneOverMaxLHGRVec );
    lhgrVec->scale( 0.2 );
    double constLHGR = 0.4;
    lhgrVec->addScalar( *lhgrVec, constLHGR );

    mechanicsNonlinearVolumeOperator->setReferenceTemperature( tempVecRef );
    mechanicsNonlinearVolumeOperator->setVector( AMP::Operator::Mechanics::TEMPERATURE, tempVec );
    mechanicsNonlinearVolumeOperator->setVector( AMP::Operator::Mechanics::BURNUP, burnVec );
    mechanicsNonlinearVolumeOperator->setVector( AMP::Operator::Mechanics::LHGR, lhgrVec );

    // Create a Linear BVP operator for mechanics
    AMP_INSIST( input_db->keyExists( "LinearMechanicsOperator" ), "key missing!" );
    auto linearMechanicsBVPoperator = std::dynamic_pointer_cast<AMP::Operator::LinearBVPOperator>(
        AMP::Operator::OperatorBuilder::createOperator(
            meshAdapter, "LinearMechanicsOperator", input_db, mechanicsMaterialModel ) );

    AMP::pout << "Created the linearMechanicsOperator." << std::endl;

    // We need to reset the linear operator before the solve since TrilinosML does
    // the factorization of the matrix during construction and so the matrix must
    // be correct before constructing the TrilinosML object.
    AMP::pout << "About to call the first apply." << std::endl;
    nonlinearMechanicsBVPoperator->apply( solVec, resVec );
    AMP::pout << "About to call the first reset." << std::endl;
    linearMechanicsBVPoperator->reset(
        nonlinearMechanicsBVPoperator->getParameters( "Jacobian", solVec ) );

    double epsilon =
        1.0e-13 *
        static_cast<double>(
            ( ( linearMechanicsBVPoperator->getMatrix() )->extractDiagonal() )->L1Norm() );

    auto nonlinearSolver_db = input_db->getDatabase( "NonlinearSolver" );
    auto linearSolver_db    = nonlinearSolver_db->getDatabase( "LinearSolver" );

    // ---- first initialize the preconditioner
    auto pcSolver_db    = linearSolver_db->getDatabase( "Preconditioner" );
    auto pcSolverParams = std::make_shared<AMP::Solver::TrilinosMLSolverParameters>( pcSolver_db );
    pcSolverParams->d_pOperator = linearMechanicsBVPoperator;
    auto pcSolver               = std::make_shared<AMP::Solver::TrilinosMLSolver>( pcSolverParams );

    // HACK to prevent a double delete on Petsc Vec
    std::shared_ptr<AMP::Solver::PetscSNESSolver> nonlinearSolver;

    // initialize the linear solver
    auto linearSolverParams =
        std::make_shared<AMP::Solver::SolverStrategyParameters>( linearSolver_db );
    linearSolverParams->d_pOperator       = linearMechanicsBVPoperator;
    linearSolverParams->d_comm            = globalComm;
    linearSolverParams->d_pPreconditioner = pcSolver;
    auto linearSolver = std::make_shared<AMP::Solver::PetscKrylovSolver>( linearSolverParams );

    // initialize the nonlinear solver
    auto nonlinearSolverParams =
        std::make_shared<AMP::Solver::SolverStrategyParameters>( nonlinearSolver_db );
    // change the next line to get the correct communicator out
    nonlinearSolverParams->d_comm          = globalComm;
    nonlinearSolverParams->d_pOperator     = nonlinearMechanicsBVPoperator;
    nonlinearSolverParams->d_pNestedSolver = linearSolver;
    nonlinearSolverParams->d_pInitialGuess = solVec;
    nonlinearSolver.reset( new AMP::Solver::PetscSNESSolver( nonlinearSolverParams ) );

    nonlinearSolver->setZeroInitialGuess( false );

    for ( int step = 0; step < NumberOfLoadingSteps; step++ ) {
        AMP::pout << "########################################" << std::endl;
        AMP::pout << "The current loading step is " << ( step + 1 ) << std::endl;

        double finalTemperature = 301.0 + ( ( (double) ( step + 1 ) ) * 100.0 );
        tempVec->setToScalar( finalTemperature );

        double scaleValue = ( (double) step + 1.0 ) / NumberOfLoadingSteps;
        scaledRhsVec->scale( scaleValue, *rhsVec );
        AMP::pout << "L2 Norm of RHS at loading step " << ( step + 1 ) << " is "
                  << scaledRhsVec->L2Norm() << std::endl;

        nonlinearMechanicsBVPoperator->residual( scaledRhsVec, solVec, resVec );
        double initialResidualNorm = static_cast<double>( resVec->L2Norm() );
        AMP::pout << "Initial Residual Norm for loading step " << ( step + 1 ) << " is "
                  << initialResidualNorm << std::endl;

        nonlinearSolver->apply( scaledRhsVec, solVec );

        nonlinearMechanicsBVPoperator->residual( scaledRhsVec, solVec, resVec );
        double finalResidualNorm = static_cast<double>( resVec->L2Norm() );
        AMP::pout << "Final Residual Norm for loading step " << ( step + 1 ) << " is "
                  << finalResidualNorm << std::endl;

        if ( finalResidualNorm > ( 1.0e-10 * initialResidualNorm ) ) {
            ut->failure( "Nonlinear solve for current loading step" );
        } else {
            ut->passes( "Nonlinear solve for current loading step" );
        }

        double finalSolNorm = static_cast<double>( solVec->L2Norm() );

        AMP::pout << "Final Solution Norm: " << finalSolNorm << std::endl;

        auto mechUvec = solVec->select( AMP::LinearAlgebra::VS_Stride( 0, 3 ), "U" );
        auto mechVvec = solVec->select( AMP::LinearAlgebra::VS_Stride( 1, 3 ), "V" );
        auto mechWvec = solVec->select( AMP::LinearAlgebra::VS_Stride( 2, 3 ), "W" );

        AMP::pout << "Maximum U displacement: " << mechUvec->maxNorm() << std::endl;
        AMP::pout << "Maximum V displacement: " << mechVvec->maxNorm() << std::endl;
        AMP::pout << "Maximum W displacement: " << mechWvec->maxNorm() << std::endl;

        auto tmp_db = std::make_shared<AMP::Database>( "Dummy" );
        auto tmpParams =
            std::make_shared<AMP::Operator::MechanicsNonlinearFEOperatorParameters>( tmp_db );
        nonlinearMechanicsBVPoperator->getVolumeOperator()->reset( tmpParams );
        nonlinearSolver->setZeroInitialGuess( false );
    }

    AMP::pout << "epsilon = " << epsilon << std::endl;

    auto siloWriter = AMP::IO::Writer::buildWriter( "Silo" );
    siloWriter->registerVector( solVec, meshAdapter, AMP::Mesh::GeomType::Vertex, "Solution" );
    siloWriter->registerVector(
        burnVec, meshAdapter, AMP::Mesh::GeomType::Vertex, "InitialDamageThreshold" );
    siloWriter->registerVector(
        lhgrVec, meshAdapter, AMP::Mesh::GeomType::Vertex, "CriticalDamageThreshold" );
    siloWriter->writeFile( "undeformedMesh_SS_3", 1 );
    meshAdapter->displaceMesh( solVec );
    siloWriter->writeFile( "deformedMesh_SS_3", 1 );

    ut->passes( exeName );
}

int testSmallStrainMechanics_ThermoMechanics_1( int argc, char *argv[] )
{
    AMP::AMPManager::startup( argc, argv );
    AMP::UnitTest ut;

    std::vector<std::string> exeNames;
    exeNames.emplace_back( "testSmallStrainMechanics-ThermoMechanics-1" );

    for ( auto &exeName : exeNames )
        myTest( &ut, exeName );

    ut.report();

    int num_failed = ut.NumFailGlobal();
    AMP::AMPManager::shutdown();
    return num_failed;
}
