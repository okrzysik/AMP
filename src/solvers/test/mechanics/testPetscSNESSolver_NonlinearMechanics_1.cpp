#include "AMP/IO/PIO.h"
#include "AMP/IO/Writer.h"
#include "AMP/discretization/simpleDOF_Manager.h"
#include "AMP/mesh/MeshParameters.h"
#include "AMP/operators/BVPOperatorParameters.h"
#include "AMP/operators/NonlinearBVPOperator.h"
#include "AMP/operators/OperatorBuilder.h"
#include "AMP/operators/boundary/DirichletVectorCorrection.h"
#include "AMP/operators/mechanics/MechanicsNonlinearFEOperator.h"
#include "AMP/solvers/SolverFactory.h"
#include "AMP/solvers/SolverStrategy.h"
#include "AMP/solvers/SolverStrategyParameters.h"
#include "AMP/utils/AMPManager.h"
#include "AMP/utils/AMP_MPI.h"
#include "AMP/utils/Database.h"
#include "AMP/utils/UnitTest.h"
#include "AMP/utils/Utilities.h"
#include "AMP/vectors/VectorBuilder.h"
#include "AMP/vectors/VectorSelector.h"

#include <iostream>
#include <string>


std::shared_ptr<AMP::Solver::SolverStrategy>
buildSolver( const std::string &solver_name,
             std::shared_ptr<AMP::Database> input_db,
             const AMP::AMP_MPI &comm,
             std::shared_ptr<AMP::LinearAlgebra::Vector> initialGuess,
             std::shared_ptr<AMP::Operator::Operator> op )
{

    AMP_INSIST( input_db->keyExists( solver_name ), "Key " + solver_name + " is missing!" );

    auto db = input_db->getDatabase( solver_name );
    AMP_INSIST( db->keyExists( "name" ), "Key name does not exist in solver database" );

    auto parameters             = std::make_shared<AMP::Solver::SolverStrategyParameters>( db );
    parameters->d_pOperator     = op;
    parameters->d_comm          = comm;
    parameters->d_pInitialGuess = initialGuess;
    parameters->d_global_db     = input_db;

    return AMP::Solver::SolverFactory::create( parameters );
}

static void myTest( AMP::UnitTest *ut, const std::string &inputName )
{
    std::string input_file = inputName;
    std::string log_file   = "output_" + inputName;

    AMP::logOnlyNodeZero( log_file );
    AMP::AMP_MPI globalComm( AMP_COMM_WORLD );

    auto siloWriter = AMP::IO::Writer::buildWriter( "Silo" );

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

    std::shared_ptr<AMP::Operator::ElementPhysicsModel> elementPhysicsModel;
    auto nonlinBvpOperator = std::dynamic_pointer_cast<AMP::Operator::NonlinearBVPOperator>(
        AMP::Operator::OperatorBuilder::createOperator(
            meshAdapter, "nonlinearMechanicsBVPOperator", input_db, elementPhysicsModel ) );

    auto displacementVariable = nonlinBvpOperator->getOutputVariable();

    // For RHS (Point Forces)
    std::shared_ptr<AMP::Operator::ElementPhysicsModel> dummyModel;
    auto dirichletLoadVecOp = std::dynamic_pointer_cast<AMP::Operator::DirichletVectorCorrection>(
        AMP::Operator::OperatorBuilder::createOperator(
            meshAdapter, "Load_Boundary", input_db, dummyModel ) );
    dirichletLoadVecOp->setVariable( displacementVariable );

    // For Initial-Guess
    auto dirichletDispInVecOp = std::dynamic_pointer_cast<AMP::Operator::DirichletVectorCorrection>(
        AMP::Operator::OperatorBuilder::createOperator(
            meshAdapter, "Displacement_Boundary", input_db, dummyModel ) );
    dirichletDispInVecOp->setVariable( displacementVariable );

    auto dofMap = AMP::Discretization::simpleDOFManager::create(
        meshAdapter, AMP::Mesh::GeomType::Vertex, 1, 3, true );

    AMP::LinearAlgebra::Vector::shared_ptr nullVec;
    auto mechNlSolVec = AMP::LinearAlgebra::createVector( dofMap, displacementVariable, true );
    auto mechNlRhsVec = mechNlSolVec->clone();
    auto mechNlResVec = mechNlSolVec->clone();
    auto mechNlScaledRhsVec = mechNlSolVec->clone();

    // Initial guess for NL solver must satisfy the displacement boundary conditions
    mechNlSolVec->setToScalar( 0.0 );
    dirichletDispInVecOp->apply( nullVec, mechNlSolVec );
    mechNlSolVec->makeConsistent( AMP::LinearAlgebra::ScatterType::CONSISTENT_SET );

    nonlinBvpOperator->apply( mechNlSolVec, mechNlResVec );

    // Point forces
    mechNlRhsVec->setToScalar( 0.0 );
    dirichletLoadVecOp->apply( nullVec, mechNlRhsVec );

    // Create the solver
    auto nonlinearSolver =
        buildSolver( "NonlinearSolver", input_db, globalComm, mechNlSolVec, nonlinBvpOperator );

    nonlinearSolver->setZeroInitialGuess( false );

    for ( int step = 0; step < NumberOfLoadingSteps; step++ ) {
        AMP::pout << "########################################" << std::endl;
        AMP::pout << "The current loading step is " << ( step + 1 ) << std::endl;

        double scaleValue = ( (double) step + 1.0 ) / NumberOfLoadingSteps;
        mechNlScaledRhsVec->scale( scaleValue, *mechNlRhsVec );
        mechNlScaledRhsVec->makeConsistent( AMP::LinearAlgebra::ScatterType::CONSISTENT_SET );
        AMP::pout << "L2 Norm of RHS at loading step " << ( step + 1 ) << " is "
                  << mechNlScaledRhsVec->L2Norm() << std::endl;

        nonlinBvpOperator->residual( mechNlScaledRhsVec, mechNlSolVec, mechNlResVec );
        double initialResidualNorm = static_cast<double>( mechNlResVec->L2Norm() );
        AMP::pout << "Initial Residual Norm for loading step " << ( step + 1 ) << " is "
                  << initialResidualNorm << std::endl;

        AMP::pout << "Starting Nonlinear Solve..." << std::endl;
        nonlinearSolver->apply( mechNlScaledRhsVec, mechNlSolVec );

        nonlinBvpOperator->residual( mechNlScaledRhsVec, mechNlSolVec, mechNlResVec );
        double finalResidualNorm = static_cast<double>( mechNlResVec->L2Norm() );
        AMP::pout << "Final Residual Norm for loading step " << ( step + 1 ) << " is "
                  << finalResidualNorm << std::endl;

        if ( finalResidualNorm > ( 1.0e-10 * initialResidualNorm ) ) {
            ut->failure( "Nonlinear solve for current loading step" );
        } else {
            ut->passes( "Nonlinear solve for current loading step" );
        }

        double finalSolNorm = static_cast<double>( mechNlSolVec->L2Norm() );

        AMP::pout << "Final Solution Norm: " << finalSolNorm << std::endl;

        auto mechUvec = mechNlSolVec->select( AMP::LinearAlgebra::VS_Stride( 0, 3 ), "U" );
        auto mechVvec = mechNlSolVec->select( AMP::LinearAlgebra::VS_Stride( 1, 3 ), "V" );
        auto mechWvec = mechNlSolVec->select( AMP::LinearAlgebra::VS_Stride( 2, 3 ), "W" );

        double finalMaxU = static_cast<double>( mechUvec->maxNorm() );
        double finalMaxV = static_cast<double>( mechVvec->maxNorm() );
        double finalMaxW = static_cast<double>( mechWvec->maxNorm() );

        AMP::pout << "Maximum U displacement: " << finalMaxU << std::endl;
        AMP::pout << "Maximum V displacement: " << finalMaxV << std::endl;
        AMP::pout << "Maximum W displacement: " << finalMaxW << std::endl;

        auto tmp_db = std::make_shared<AMP::Database>( "Dummy" );
        auto tmpParams =
            std::make_shared<AMP::Operator::MechanicsNonlinearFEOperatorParameters>( tmp_db );
        ( nonlinBvpOperator->getVolumeOperator() )->reset( tmpParams );
        nonlinearSolver->setZeroInitialGuess( false );

        meshAdapter->displaceMesh( mechNlSolVec );

        siloWriter->registerVector(
            mechNlSolVec, meshAdapter, AMP::Mesh::GeomType::Vertex, "Solution" );
        auto outFileName = AMP::Utilities::stringf(
            "LoadPrescribed-DeformedPlateWithHole-LinearElasticity_%d", step );
        siloWriter->writeFile( outFileName, 0 );
    }

    double finalSolNorm = static_cast<double>( mechNlSolVec->L2Norm() );
    AMP::pout << "Final Solution Norm: " << finalSolNorm << std::endl;

    ut->passes( inputName );
}

int testPetscSNESSolver_NonlinearMechanics_1( int argc, char *argv[] )
{
    AMP::AMPManager::startup( argc, argv );
    AMP::UnitTest ut;

    std::vector<std::string> inputNames;
    inputNames.emplace_back( "input_testPetscSNESSolver-NonlinearMechanics-PlateWithHole-1" );
    inputNames.emplace_back( "input_testPetscSNESSolver-LU-NonlinearMechanics-1-normal" );
    inputNames.emplace_back( "input_testPetscSNESSolver-ML-NonlinearMechanics-1-normal" );
    inputNames.emplace_back( "input_testPetscSNESSolver-LU-NonlinearMechanics-1-reduced" );
    inputNames.emplace_back( "input_testPetscSNESSolver-ML-NonlinearMechanics-1-reduced" );

    for ( auto &inputName : inputNames )
        myTest( &ut, inputName );

    ut.report();

    int num_failed = ut.NumFailGlobal();
    AMP::AMPManager::shutdown();
    return num_failed;
}
