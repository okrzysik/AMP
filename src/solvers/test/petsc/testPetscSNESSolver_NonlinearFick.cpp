#include "AMP/ampmesh/Mesh.h"
#include "AMP/ampmesh/MeshParameters.h"
#include "AMP/discretization/DOF_Manager.h"
#include "AMP/discretization/simpleDOF_Manager.h"

#include "AMP/operators/BVPOperatorParameters.h"
#include "AMP/operators/ColumnOperator.h"
#include "AMP/operators/LinearBVPOperator.h"
#include "AMP/operators/NeutronicsRhs.h"
#include "AMP/operators/NonlinearBVPOperator.h"
#include "AMP/operators/OperatorBuilder.h"
#include "AMP/operators/boundary/DirichletVectorCorrection.h"
#include "AMP/operators/diffusion/DiffusionLinearFEOperator.h"
#include "AMP/operators/diffusion/DiffusionNonlinearFEOperator.h"
#include "AMP/operators/libmesh/VolumeIntegralOperator.h"
#include "AMP/operators/mechanics/MechanicsLinearFEOperator.h"
#include "AMP/operators/mechanics/MechanicsNonlinearFEOperator.h"
#include "AMP/solvers/ColumnSolver.h"
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
#include "AMP/vectors/MultiVariable.h"
#include "AMP/vectors/Variable.h"
#include "AMP/vectors/Vector.h"
#include "AMP/vectors/VectorBuilder.h"

#include <iostream>
#include <memory>
#include <string>


static void myTest( AMP::UnitTest *ut, const std::string &exeName )
{
    std::string input_file = "input_" + exeName;
    std::string log_file   = "output_" + exeName;

    AMP::logOnlyNodeZero( log_file );
    AMP::AMP_MPI globalComm( AMP_COMM_WORLD );

    // Read the input file
    auto input_db = AMP::Database::parseInputFile( input_file );
    input_db->print( AMP::plog );

    // Get the Mesh database and create the mesh parameters
    auto database = input_db->getDatabase( "Mesh" );
    auto params   = std::make_shared<AMP::Mesh::MeshParameters>( database );
    params->setComm( globalComm );

    // Create the meshes from the input database
    auto meshAdapter = AMP::Mesh::Mesh::buildMesh( params );

    auto DOF_scalar = AMP::Discretization::simpleDOFManager::create(
        meshAdapter, AMP::Mesh::GeomType::Vertex, 1, 1, true );

    // create a nonlinear BVP operator for nonlinear fick diffusion
    AMP_INSIST( input_db->keyExists( "testNonlinearFickOperator" ), "key missing!" );

    std::shared_ptr<AMP::Operator::ElementPhysicsModel> fickTransportModel;
    auto nonlinearFickOperator = std::dynamic_pointer_cast<AMP::Operator::NonlinearBVPOperator>(
        AMP::Operator::OperatorBuilder::createOperator(
            meshAdapter, "testNonlinearFickOperator", input_db, fickTransportModel ) );

    // initialize the input variable
    auto fickVolumeOperator =
        std::dynamic_pointer_cast<AMP::Operator::DiffusionNonlinearFEOperator>(
            nonlinearFickOperator->getVolumeOperator() );
    auto fickVariable = fickVolumeOperator->getOutputVariable();

    // create solution, rhs, and residual vectors
    auto solVec = AMP::LinearAlgebra::createVector( DOF_scalar, fickVariable );
    auto rhsVec = AMP::LinearAlgebra::createVector( DOF_scalar, fickVariable );
    auto resVec = AMP::LinearAlgebra::createVector( DOF_scalar, fickVariable );

// register some variables for plotting
#ifdef USE_EXT_SILO
    auto siloWriter = AMP::Utilities::Writer::buildWriter( "Silo" );
    siloWriter->registerVector( solVec, meshAdapter, AMP::Mesh::GeomType::Vertex, "Solution" );
    siloWriter->registerVector( resVec, meshAdapter, AMP::Mesh::GeomType::Vertex, "Residual" );
#endif

    // now construct the linear BVP operator for fick
    AMP_INSIST( input_db->keyExists( "testLinearFickOperator" ), "key missing!" );
    auto linearFickOperator = std::dynamic_pointer_cast<AMP::Operator::LinearBVPOperator>(
        AMP::Operator::OperatorBuilder::createOperator(
            meshAdapter, "testLinearFickOperator", input_db, fickTransportModel ) );

    // Initial guess
    solVec->setToScalar( .05 );
    std::cout << "initial guess norm = " << solVec->L2Norm() << "\n";
    nonlinearFickOperator->modifyInitialSolutionVector( solVec );
    std::cout << "initial guess norm  after apply = " << solVec->L2Norm() << "\n";

    rhsVec->setToScalar( 0.0 );
    nonlinearFickOperator->modifyRHSvector( rhsVec );

    auto nonlinearSolver_db = input_db->getDatabase( "NonlinearSolver" );
    auto linearSolver_db    = nonlinearSolver_db->getDatabase( "LinearSolver" );

    // initialize the nonlinear solver
    auto nonlinearSolverParams =
        std::make_shared<AMP::Solver::PetscSNESSolverParameters>( nonlinearSolver_db );

    // change the next line to get the correct communicator out
    nonlinearSolverParams->d_comm          = globalComm;
    nonlinearSolverParams->d_pOperator     = nonlinearFickOperator;
    nonlinearSolverParams->d_pInitialGuess = solVec;
    auto nonlinearSolver = std::make_shared<AMP::Solver::PetscSNESSolver>( nonlinearSolverParams );

    auto fickPreconditioner_db = linearSolver_db->getDatabase( "Preconditioner" );
    auto fickPreconditionerParams =
        std::make_shared<AMP::Solver::SolverStrategyParameters>( fickPreconditioner_db );
    fickPreconditionerParams->d_pOperator = linearFickOperator;
    auto linearFickPreconditioner =
        std::make_shared<AMP::Solver::TrilinosMLSolver>( fickPreconditionerParams );

    // register the preconditioner with the Jacobian free Krylov solver
    auto linearSolver = nonlinearSolver->getKrylovSolver();

    linearSolver->setPreconditioner( linearFickPreconditioner );

    nonlinearFickOperator->residual( rhsVec, solVec, resVec );

    AMP::pout << "Initial Residual Norm: " << resVec->L2Norm() << std::endl;

    nonlinearSolver->setZeroInitialGuess( false );

    nonlinearSolver->apply( rhsVec, solVec );

    nonlinearFickOperator->residual( rhsVec, solVec, resVec );

    double finalResidualNorm = static_cast<double>( resVec->L2Norm() );

    std::cout << "Final Residual Norm: " << finalResidualNorm << std::endl;

    solVec->makeConsistent( AMP::LinearAlgebra::VectorData::ScatterType::CONSISTENT_SET );
    resVec->makeConsistent( AMP::LinearAlgebra::VectorData::ScatterType::CONSISTENT_SET );

#ifdef USE_EXT_SILO
    siloWriter->writeFile( exeName, 0 );
#endif

    if ( finalResidualNorm > 1.0e-08 ) {
        ut->failure( exeName );
    } else {
        ut->passes( exeName );
    }
}

int testPetscSNESSolver_NonlinearFick( int argc, char *argv[] )
{
    AMP::AMPManager::startup( argc, argv );
    AMP::UnitTest ut;

    std::vector<std::string> exeNames;
    exeNames.emplace_back( "testPetscSNESSolver-NonlinearFick-cylinder-1a" );
    exeNames.emplace_back( "testPetscSNESSolver-NonlinearFick-cylinder-1b" );
    exeNames.emplace_back( "testPetscSNESSolver-NonlinearFick-cylinder-1c" );
    exeNames.emplace_back( "testPetscSNESSolver-NonlinearFick-cylinder-1d" );

    for ( auto &exeName : exeNames )
        myTest( &ut, exeName );

    ut.report();

    int num_failed = ut.NumFailGlobal();
    AMP::AMPManager::shutdown();
    return num_failed;
}
