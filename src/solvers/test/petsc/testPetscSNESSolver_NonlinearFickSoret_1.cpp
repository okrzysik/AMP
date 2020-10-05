#include "AMP/ampmesh/Mesh.h"
#include "AMP/discretization/DOF_Manager.h"
#include "AMP/discretization/simpleDOF_Manager.h"
#include "AMP/materials/Material.h"
#include "AMP/operators/BVPOperatorParameters.h"
#include "AMP/operators/ColumnOperator.h"
#include "AMP/operators/LinearBVPOperator.h"
#include "AMP/operators/NeutronicsRhs.h"
#include "AMP/operators/NonlinearBVPOperator.h"
#include "AMP/operators/OperatorBuilder.h"
#include "AMP/operators/boundary/DirichletVectorCorrection.h"
#include "AMP/operators/diffusion/DiffusionLinearFEOperator.h"
#include "AMP/operators/diffusion/DiffusionNonlinearFEOperator.h"
#include "AMP/operators/diffusion/FickSoretNonlinearFEOperator.h"
#include "AMP/operators/libmesh/VolumeIntegralOperator.h"
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
#include "AMP/vectors/Variable.h"
#include "AMP/vectors/Vector.h"
#include "AMP/vectors/VectorBuilder.h"
#include "AMP/vectors/VectorSelector.h"
#include <memory>

#include <cmath>
#include <iostream>
#include <string>


static void fickTest( AMP::UnitTest *ut, std::string exeName, std::vector<double> &results )
{
    std::string input_file = "input_" + exeName;
    std::string log_file   = "output_" + exeName;

    AMP::PIO::logOnlyNodeZero( log_file );
    AMP::AMP_MPI globalComm( AMP_COMM_WORLD );

    auto input_db = AMP::Database::parseInputFile( input_file );
    input_db->print( AMP::plog );

    // Get the Mesh database and create the mesh parameters
    auto database = input_db->getDatabase( "Mesh" );
    auto params   = std::make_shared<AMP::Mesh::MeshParameters>( database );
    params->setComm( globalComm );

    // Create the meshes from the input database
    auto manager     = AMP::Mesh::Mesh::buildMesh( params );
    auto meshAdapter = manager->Subset( "cylinder" );

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
    auto nodalScalarDOF = AMP::Discretization::simpleDOFManager::create(
        meshAdapter, AMP::Mesh::GeomType::Vertex, 1, 1, true );
    auto solVec = AMP::LinearAlgebra::createVector( nodalScalarDOF, fickVariable, true );
    auto rhsVec = AMP::LinearAlgebra::createVector( nodalScalarDOF, fickVariable, true );
    auto resVec = AMP::LinearAlgebra::createVector( nodalScalarDOF, fickVariable, true );

#ifdef USE_EXT_SILO
    // register some variables for plotting
    AMP::Utilities::Writer::shared_ptr siloWriter = AMP::Utilities::Writer::buildWriter( "Silo" );
    siloWriter->registerVector( solVec, meshAdapter, AMP::Mesh::GeomType::Vertex, "Solution" );
    siloWriter->registerVector( resVec, meshAdapter, AMP::Mesh::GeomType::Vertex, "Residual" );
#endif

    // now construct the linear BVP operator for fick
    AMP_INSIST( input_db->keyExists( "testLinearFickOperator" ), "key missing!" );
    std::shared_ptr<AMP::Operator::LinearBVPOperator> linearFickOperator =
        std::dynamic_pointer_cast<AMP::Operator::LinearBVPOperator>(
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
    nonlinearSolver->solve( rhsVec, solVec );
    nonlinearFickOperator->residual( rhsVec, solVec, resVec );
    std::cout << "Final Residual Norm: " << resVec->L2Norm() << std::endl;

    solVec->makeConsistent( AMP::LinearAlgebra::VectorData::ScatterType::CONSISTENT_SET );
    resVec->makeConsistent( AMP::LinearAlgebra::VectorData::ScatterType::CONSISTENT_SET );

#ifdef USE_EXT_SILO
    siloWriter->writeFile( exeName, 0 );
#endif

    // store result
    {
        AMP::Mesh::MeshIterator iterator =
            meshAdapter->getIterator( AMP::Mesh::GeomType::Vertex, 0 );
        size_t numNodes = iterator.size();
        results.resize( numNodes );
        std::vector<size_t> dofs;
        for ( size_t iNode = 0; iNode < numNodes; iNode++ ) {
            nodalScalarDOF->getDOFs( iterator->globalID(), dofs );
            size_t gid     = dofs[0];
            results[iNode] = solVec->getValueByGlobalID( gid );
            ++iterator;
        }
    }

    ut->passes( exeName );
}


static void fickSoretTest( AMP::UnitTest *ut, std::string exeName, std::vector<double> &results )
{
    std::string input_file = "input_" + exeName;
    std::string log_file   = "output_" + exeName;

    AMP::PIO::logOnlyNodeZero( log_file );
    AMP::AMP_MPI globalComm( AMP_COMM_WORLD );

    auto input_db = AMP::Database::parseInputFile( input_file );
    input_db->print( AMP::plog );

    // Get the Mesh database and create the mesh parameters
    auto database = input_db->getDatabase( "Mesh" );
    auto params   = std::make_shared<AMP::Mesh::MeshParameters>( database );
    params->setComm( globalComm );

    // Create the meshes from the input database
    auto manager     = AMP::Mesh::Mesh::buildMesh( params );
    auto meshAdapter = manager->Subset( "cylinder" );

    // create a nonlinear BVP operator for nonlinear Fick-Soret diffusion
    AMP_INSIST( input_db->keyExists( "testNonlinearFickSoretBVPOperator" ), "key missing!" );

    // Create nonlinear FickSoret BVP operator and access volume nonlinear FickSoret operator
    std::shared_ptr<AMP::Operator::ElementPhysicsModel> elementPhysicsModel;
    auto nlinBVPOperator = AMP::Operator::OperatorBuilder::createOperator(
        meshAdapter, "testNonlinearFickSoretBVPOperator", input_db, elementPhysicsModel );
    auto nlinBVPOp =
        std::dynamic_pointer_cast<AMP::Operator::NonlinearBVPOperator>( nlinBVPOperator );
    auto nlinOp = std::dynamic_pointer_cast<AMP::Operator::FickSoretNonlinearFEOperator>(
        nlinBVPOp->getVolumeOperator() );
    auto fickOp = std::dynamic_pointer_cast<AMP::Operator::DiffusionNonlinearFEOperator>(
        nlinOp->getFickOperator() );
    auto soretOp = std::dynamic_pointer_cast<AMP::Operator::DiffusionNonlinearFEOperator>(
        nlinOp->getSoretOperator() );

    // use the linear BVP operator to create a Fick linear operator with bc's
    AMP_INSIST( input_db->keyExists( "testLinearFickBVPOperator" ), "key missing!" );

    auto linBVPOperator = AMP::Operator::OperatorBuilder::createOperator(
        meshAdapter, "testLinearFickBVPOperator", input_db, elementPhysicsModel );
    auto linBVPOp = std::dynamic_pointer_cast<AMP::Operator::LinearBVPOperator>( linBVPOperator );

    auto tVar = std::make_shared<AMP::LinearAlgebra::Variable>( "temp" );
    auto cVar = std::make_shared<AMP::LinearAlgebra::Variable>( *fickOp->getOutputVariable() );
    auto fsOutVar =
        std::make_shared<AMP::LinearAlgebra::Variable>( *nlinBVPOp->getOutputVariable() );

    // create solution, rhs, and residual vectors
    auto nodalScalarDOF = AMP::Discretization::simpleDOFManager::create(
        meshAdapter, AMP::Mesh::GeomType::Vertex, 1, 1, true );
    auto solVec = AMP::LinearAlgebra::createVector( nodalScalarDOF, cVar, true );
    auto rhsVec = AMP::LinearAlgebra::createVector( nodalScalarDOF, fsOutVar, true );
    auto resVec = AMP::LinearAlgebra::createVector( nodalScalarDOF, fsOutVar, true );

    // create parameters for reset test and reset fick and soret operators
    auto tVec = AMP::LinearAlgebra::createVector( nodalScalarDOF, tVar, true );
    tVec->setToScalar( 300. );

    fickOp->setVector( 0, tVec );
    soretOp->setVector( 0, tVec );

#ifdef USE_EXT_SILO
    // register some variables for plotting
    auto siloWriter = AMP::Utilities::Writer::buildWriter( "Silo" );
    siloWriter->registerVector( solVec, meshAdapter, AMP::Mesh::GeomType::Vertex, "Solution" );
    siloWriter->registerVector( resVec, meshAdapter, AMP::Mesh::GeomType::Vertex, "Residual" );
#endif

    // Initial guess

    solVec->setToScalar( .05 );
    std::cout << "initial guess norm = " << solVec->L2Norm() << "\n";
    nlinBVPOp->modifyInitialSolutionVector( solVec );
    std::cout << "initial guess norm  after apply = " << solVec->L2Norm() << "\n";

    rhsVec->setToScalar( 0.0 );
    nlinBVPOp->modifyRHSvector( rhsVec );

    auto nonlinearSolver_db = input_db->getDatabase( "NonlinearSolver" );
    auto linearSolver_db    = nonlinearSolver_db->getDatabase( "LinearSolver" );

    // initialize the nonlinear solver
    auto nonlinearSolverParams =
        std::make_shared<AMP::Solver::PetscSNESSolverParameters>( nonlinearSolver_db );

    // change the next line to get the correct communicator out
    nonlinearSolverParams->d_comm          = globalComm;
    nonlinearSolverParams->d_pOperator     = nlinBVPOp;
    nonlinearSolverParams->d_pInitialGuess = solVec;

    auto nonlinearSolver = std::make_shared<AMP::Solver::PetscSNESSolver>( nonlinearSolverParams );

    auto fickPreconditioner_db = linearSolver_db->getDatabase( "Preconditioner" );
    auto fickPreconditionerParams =
        std::make_shared<AMP::Solver::SolverStrategyParameters>( fickPreconditioner_db );
    fickPreconditionerParams->d_pOperator = linBVPOp;
    auto linearFickPreconditioner =
        std::make_shared<AMP::Solver::TrilinosMLSolver>( fickPreconditionerParams );

    // register the preconditioner with the Jacobian free Krylov solver
    auto linearSolver = nonlinearSolver->getKrylovSolver();

    linearSolver->setPreconditioner( linearFickPreconditioner );

    nlinBVPOp->residual( rhsVec, solVec, resVec );
    AMP::pout << "Initial Residual Norm: " << resVec->L2Norm() << std::endl;

    nonlinearSolver->setZeroInitialGuess( false );
    nonlinearSolver->solve( rhsVec, solVec );
    nlinBVPOp->residual( rhsVec, solVec, resVec );
    std::cout << "Final Residual Norm: " << resVec->L2Norm() << std::endl;

    solVec->makeConsistent( AMP::LinearAlgebra::VectorData::ScatterType::CONSISTENT_SET );
    resVec->makeConsistent( AMP::LinearAlgebra::VectorData::ScatterType::CONSISTENT_SET );

#ifdef USE_EXT_SILO
    siloWriter->writeFile( exeName, 0 );
#endif

    // store result
    {
        auto iterator   = meshAdapter->getIterator( AMP::Mesh::GeomType::Vertex, 0 );
        size_t numNodes = iterator.size();
        results.resize( numNodes );
        std::vector<size_t> dofs;
        for ( size_t iNode = 0; iNode < numNodes; iNode++ ) {
            nodalScalarDOF->getDOFs( iterator->globalID(), dofs );
            size_t gid     = dofs[0];
            results[iNode] = solVec->getValueByGlobalID( gid );
            ++iterator;
        }
    }

    ut->passes( exeName );
}


int testPetscSNESSolver_NonlinearFickSoret_1( int argc, char *argv[] )
{
    AMP::AMPManager::startup( argc, argv );
    AMP::UnitTest ut;

    std::vector<double> fickOnly, fickSoretOff, fickSoretZero, fickOnlyReal, fickSoretOffReal;

    fickTest( &ut, "testPetscSNESSolver-NonlinearFick-cylinder-TUI-1", fickOnly );
    fickSoretTest( &ut, "testPetscSNESSolver-NonlinearFickSoret-cylinder-TUI-1", fickSoretOff );
    fickSoretTest( &ut, "testPetscSNESSolver-NonlinearFickSoret-cylinder-TUI-2", fickSoretZero );
    fickTest( &ut, "testPetscSNESSolver-NonlinearFick-cylinder-TUI-2", fickOnlyReal );
    fickSoretTest( &ut, "testPetscSNESSolver-NonlinearFickSoret-cylinder-TUI-3", fickSoretOffReal );
    AMP_INSIST( fickOnly.size() == fickSoretOff.size() and
                    fickSoretOff.size() == fickSoretZero.size() and
                    fickOnlyReal.size() == fickSoretOffReal.size(),
                "sizes of results do not match" );

    double l2err1 = 0., l2err2 = 0.;
    for ( size_t i = 0; i < fickOnly.size(); i++ ) {
        double err = fickOnly[i] - fickSoretOff[i];
        l2err1 += err * err;
        err = fickSoretOff[i] - fickSoretZero[i];
        l2err2 += err * err;
    }
    l2err1 = sqrt( l2err1 );
    l2err2 = sqrt( l2err2 );

    std::cout << "fick/soretOff err = " << l2err1 << "  soretOff/soretZero err = " << l2err2
              << std::endl;

    double l2err3 = 0.;
    for ( size_t i = 0; i < fickOnlyReal.size(); i++ ) {
        double err = fickOnlyReal[i] - fickSoretOffReal[i];
        l2err3 += err * err;
    }
    l2err3 = sqrt( l2err3 );

    std::cout << "fick/soretOff real err = " << l2err3 << std::endl;

    if ( ( l2err1 < 1.e-6 ) && ( l2err2 < 1.e-6 ) && ( l2err3 < 1.e-6 ) ) {
        ut.passes( "fick, fick-soret/off, and fick-soret/zero all agree" );
    }

    ut.report();

    int num_failed = ut.NumFailGlobal();
    AMP::AMPManager::shutdown();
    return num_failed;
}
