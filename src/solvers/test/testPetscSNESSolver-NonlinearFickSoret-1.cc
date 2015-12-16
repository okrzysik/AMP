#include "utils/AMPManager.h"
#include "utils/UnitTest.h"
#include "utils/Utilities.h"

#include <cmath>
#include <iostream>
#include <string>

#include "utils/shared_ptr.h"

#include "operators/NeutronicsRhs.h"
#include "operators/libmesh/VolumeIntegralOperator.h"

#include "materials/Material.h"
#include "utils/AMPManager.h"
#include "utils/AMP_MPI.h"
#include "utils/Database.h"
#include "utils/InputDatabase.h"
#include "utils/InputManager.h"
#include "utils/PIO.h"

#include "ampmesh/Mesh.h"
#include "utils/Writer.h"

#include "discretization/DOF_Manager.h"
#include "discretization/simpleDOF_Manager.h"
#include "vectors/Variable.h"
#include "vectors/Variable.h"
#include "vectors/Vector.h"
#include "vectors/VectorBuilder.h"
#include "vectors/VectorSelector.h"

#include "operators/boundary/DirichletVectorCorrection.h"

#include "operators/BVPOperatorParameters.h"
#include "operators/ColumnOperator.h"
#include "operators/LinearBVPOperator.h"
#include "operators/NonlinearBVPOperator.h"
#include "operators/OperatorBuilder.h"
#include "operators/diffusion/DiffusionLinearFEOperator.h"
#include "operators/diffusion/DiffusionNonlinearFEOperator.h"
#include "operators/diffusion/FickSoretNonlinearFEOperator.h"

#include "solvers/ColumnSolver.h"
#include "solvers/petsc/PetscKrylovSolver.h"
#include "solvers/petsc/PetscKrylovSolverParameters.h"
#include "solvers/petsc/PetscSNESSolver.h"
#include "solvers/petsc/PetscSNESSolverParameters.h"

#include "solvers/trilinos/TrilinosMLSolver.h"


void fickTest( AMP::UnitTest *ut, std::string exeName, std::vector<double> &results )
{
    std::string input_file = "input_" + exeName;
    std::string log_file   = "output_" + exeName;

    AMP::PIO::logOnlyNodeZero( log_file );
    AMP::AMP_MPI globalComm( AMP_COMM_WORLD );

    AMP::shared_ptr<AMP::InputDatabase> input_db( new AMP::InputDatabase( "input_db" ) );
    AMP::InputManager::getManager()->parseInputFile( input_file, input_db );
    input_db->printClassData( AMP::plog );

    // Get the Mesh database and create the mesh parameters
    AMP::shared_ptr<AMP::Database> database = input_db->getDatabase( "Mesh" );
    AMP::shared_ptr<AMP::Mesh::MeshParameters> params( new AMP::Mesh::MeshParameters( database ) );
    params->setComm( globalComm );

    // Create the meshes from the input database
    AMP::Mesh::Mesh::shared_ptr manager     = AMP::Mesh::Mesh::buildMesh( params );
    AMP::Mesh::Mesh::shared_ptr meshAdapter = manager->Subset( "cylinder" );

    //----------------------------------------------------------------------------------------------------------------------------------------------//
    // create a nonlinear BVP operator for nonlinear fick diffusion
    AMP_INSIST( input_db->keyExists( "testNonlinearFickOperator" ), "key missing!" );

    AMP::shared_ptr<AMP::Operator::ElementPhysicsModel> fickTransportModel;
    AMP::shared_ptr<AMP::Operator::NonlinearBVPOperator> nonlinearFickOperator =
        AMP::dynamic_pointer_cast<AMP::Operator::NonlinearBVPOperator>(
            AMP::Operator::OperatorBuilder::createOperator(
                meshAdapter, "testNonlinearFickOperator", input_db, fickTransportModel ) );

    //----------------------------------------------------------------------------------------------------------------------------------------------//
    // initialize the input variable
    AMP::shared_ptr<AMP::Operator::DiffusionNonlinearFEOperator> fickVolumeOperator =
        AMP::dynamic_pointer_cast<AMP::Operator::DiffusionNonlinearFEOperator>(
            nonlinearFickOperator->getVolumeOperator() );

    AMP::shared_ptr<AMP::LinearAlgebra::Variable> fickVariable =
        fickVolumeOperator->getOutputVariable();

    // create solution, rhs, and residual vectors
    AMP::Discretization::DOFManager::shared_ptr nodalScalarDOF =
        AMP::Discretization::simpleDOFManager::create( meshAdapter, AMP::Mesh::Vertex, 1, 1, true );
    AMP::LinearAlgebra::Vector::shared_ptr solVec =
        AMP::LinearAlgebra::createVector( nodalScalarDOF, fickVariable, true );
    AMP::LinearAlgebra::Vector::shared_ptr rhsVec =
        AMP::LinearAlgebra::createVector( nodalScalarDOF, fickVariable, true );
    AMP::LinearAlgebra::Vector::shared_ptr resVec =
        AMP::LinearAlgebra::createVector( nodalScalarDOF, fickVariable, true );

#ifdef USE_EXT_SILO
    //----------------------------------------------------------------------------------------------------------------------------------------------//
    // register some variables for plotting
    AMP::Utilities::Writer::shared_ptr siloWriter = AMP::Utilities::Writer::buildWriter( "Silo" );
    siloWriter->registerVector( solVec, meshAdapter, AMP::Mesh::Vertex, "Solution" );
    siloWriter->registerVector( resVec, meshAdapter, AMP::Mesh::Vertex, "Residual" );
#endif

    //----------------------------------------------------------------------------------------------------------------------------------------------//
    // now construct the linear BVP operator for fick
    AMP_INSIST( input_db->keyExists( "testLinearFickOperator" ), "key missing!" );
    AMP::shared_ptr<AMP::Operator::LinearBVPOperator> linearFickOperator =
        AMP::dynamic_pointer_cast<AMP::Operator::LinearBVPOperator>(
            AMP::Operator::OperatorBuilder::createOperator(
                meshAdapter, "testLinearFickOperator", input_db, fickTransportModel ) );

    //----------------------------------------------------------------------------------------------------------------------------------------------//
    // Initial guess

    solVec->setToScalar( .05 );
    double initialGuessNorm = solVec->L2Norm();
    std::cout << "initial guess norm = " << initialGuessNorm << "\n";

    nonlinearFickOperator->modifyInitialSolutionVector( solVec );

    initialGuessNorm = solVec->L2Norm();
    std::cout << "initial guess norm  after apply = " << initialGuessNorm << "\n";

    rhsVec->setToScalar( 0.0 );
    nonlinearFickOperator->modifyRHSvector( rhsVec );

    //----------------------------------------------------------------------------------------------------------------------------------------------/

    AMP::shared_ptr<AMP::Database> nonlinearSolver_db = input_db->getDatabase( "NonlinearSolver" );
    AMP::shared_ptr<AMP::Database> linearSolver_db =
        nonlinearSolver_db->getDatabase( "LinearSolver" );

    //----------------------------------------------------------------------------------------------------------------------------------------------//
    // initialize the nonlinear solver
    AMP::shared_ptr<AMP::Solver::PetscSNESSolverParameters> nonlinearSolverParams(
        new AMP::Solver::PetscSNESSolverParameters( nonlinearSolver_db ) );

    // change the next line to get the correct communicator out
    nonlinearSolverParams->d_comm          = globalComm;
    nonlinearSolverParams->d_pOperator     = nonlinearFickOperator;
    nonlinearSolverParams->d_pInitialGuess = solVec;

    AMP::shared_ptr<AMP::Solver::PetscSNESSolver> nonlinearSolver(
        new AMP::Solver::PetscSNESSolver( nonlinearSolverParams ) );

    //----------------------------------------------------------------------------------------------------------------------------------------------//
    AMP::shared_ptr<AMP::Database> fickPreconditioner_db =
        linearSolver_db->getDatabase( "Preconditioner" );
    AMP::shared_ptr<AMP::Solver::SolverStrategyParameters> fickPreconditionerParams(
        new AMP::Solver::SolverStrategyParameters( fickPreconditioner_db ) );
    fickPreconditionerParams->d_pOperator = linearFickOperator;
    AMP::shared_ptr<AMP::Solver::TrilinosMLSolver> linearFickPreconditioner(
        new AMP::Solver::TrilinosMLSolver( fickPreconditionerParams ) );

    //----------------------------------------------------------------------------------------------------------------------------------------------//
    // register the preconditioner with the Jacobian free Krylov solver
    AMP::shared_ptr<AMP::Solver::PetscKrylovSolver> linearSolver =
        nonlinearSolver->getKrylovSolver();

    linearSolver->setPreconditioner( linearFickPreconditioner );

    nonlinearFickOperator->residual( rhsVec, solVec, resVec );
    double initialResidualNorm = resVec->L2Norm();

    AMP::pout << "Initial Residual Norm: " << initialResidualNorm << std::endl;

    nonlinearSolver->setZeroInitialGuess( false );

    nonlinearSolver->solve( rhsVec, solVec );

    nonlinearFickOperator->residual( rhsVec, solVec, resVec );

    double finalResidualNorm = resVec->L2Norm();

    std::cout << "Final Residual Norm: " << finalResidualNorm << std::endl;

    solVec->makeConsistent( AMP::LinearAlgebra::Vector::CONSISTENT_SET );
    resVec->makeConsistent( AMP::LinearAlgebra::Vector::CONSISTENT_SET );

#ifdef USE_EXT_SILO
    siloWriter->writeFile( exeName, 0 );
#endif

    // store result
    {
        AMP::Mesh::MeshIterator iterator = meshAdapter->getIterator( AMP::Mesh::Vertex, 0 );
        size_t numNodes                  = iterator.size();
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


void fickSoretTest( AMP::UnitTest *ut, std::string exeName, std::vector<double> &results )
{
    std::string input_file = "input_" + exeName;
    std::string log_file   = "output_" + exeName;

    AMP::PIO::logOnlyNodeZero( log_file );
    AMP::AMP_MPI globalComm( AMP_COMM_WORLD );

    AMP::shared_ptr<AMP::InputDatabase> input_db( new AMP::InputDatabase( "input_db" ) );
    AMP::InputManager::getManager()->parseInputFile( input_file, input_db );
    input_db->printClassData( AMP::plog );

    // Get the Mesh database and create the mesh parameters
    AMP::shared_ptr<AMP::Database> database = input_db->getDatabase( "Mesh" );
    AMP::shared_ptr<AMP::Mesh::MeshParameters> params( new AMP::Mesh::MeshParameters( database ) );
    params->setComm( globalComm );

    // Create the meshes from the input database
    AMP::Mesh::Mesh::shared_ptr manager     = AMP::Mesh::Mesh::buildMesh( params );
    AMP::Mesh::Mesh::shared_ptr meshAdapter = manager->Subset( "cylinder" );

    //----------------------------------------------------------------------------------------------------------------------------------------------//
    // create a nonlinear BVP operator for nonlinear Fick-Soret diffusion
    AMP_INSIST( input_db->keyExists( "testNonlinearFickSoretBVPOperator" ), "key missing!" );

    // Create nonlinear FickSoret BVP operator and access volume nonlinear FickSoret operator
    AMP::shared_ptr<AMP::Operator::ElementPhysicsModel> elementPhysicsModel;
    AMP::shared_ptr<AMP::Operator::Operator> nlinBVPOperator =
        AMP::Operator::OperatorBuilder::createOperator(
            meshAdapter, "testNonlinearFickSoretBVPOperator", input_db, elementPhysicsModel );
    AMP::shared_ptr<AMP::Operator::NonlinearBVPOperator> nlinBVPOp =
        AMP::dynamic_pointer_cast<AMP::Operator::NonlinearBVPOperator>( nlinBVPOperator );
    AMP::shared_ptr<AMP::Operator::FickSoretNonlinearFEOperator> nlinOp =
        AMP::dynamic_pointer_cast<AMP::Operator::FickSoretNonlinearFEOperator>(
            nlinBVPOp->getVolumeOperator() );
    AMP::shared_ptr<AMP::Operator::DiffusionNonlinearFEOperator> fickOp =
        AMP::dynamic_pointer_cast<AMP::Operator::DiffusionNonlinearFEOperator>(
            nlinOp->getFickOperator() );
    AMP::shared_ptr<AMP::Operator::DiffusionNonlinearFEOperator> soretOp =
        AMP::dynamic_pointer_cast<AMP::Operator::DiffusionNonlinearFEOperator>(
            nlinOp->getSoretOperator() );

    //----------------------------------------------------------------------------------------------------------------------------------------------//
    // use the linear BVP operator to create a Fick linear operator with bc's
    AMP_INSIST( input_db->keyExists( "testLinearFickBVPOperator" ), "key missing!" );

    AMP::shared_ptr<AMP::Operator::Operator> linBVPOperator =
        AMP::Operator::OperatorBuilder::createOperator(
            meshAdapter, "testLinearFickBVPOperator", input_db, elementPhysicsModel );
    AMP::shared_ptr<AMP::Operator::LinearBVPOperator> linBVPOp =
        AMP::dynamic_pointer_cast<AMP::Operator::LinearBVPOperator>( linBVPOperator );
    // AMP::shared_ptr<AMP::Operator::DiffusionLinearFEOperator> linOp =
    //		 AMP::dynamic_pointer_cast<AMP::Operator::DiffusionLinearFEOperator>(linBVPOp->getVolumeOperator());

    //----------------------------------------------------------------------------------------------------------------------------------------------//
    // Set up input and output variables
    // AMP::LinearAlgebra::Variable::shared_ptr
    // tVar(fickOp->getInputVariable(AMP::Operator::Diffusion::TEMPERATURE));
    // AMP::LinearAlgebra::Variable::shared_ptr
    // cVar(fickOp->getInputVariable(AMP::Operator::Diffusion::CONCENTRATION));
    AMP::LinearAlgebra::Variable::shared_ptr tVar( new AMP::LinearAlgebra::Variable( "temp" ) );
    AMP::LinearAlgebra::Variable::shared_ptr cVar( fickOp->getOutputVariable() );
    AMP::shared_ptr<AMP::LinearAlgebra::Variable> fsOutVar( nlinBVPOp->getOutputVariable() );

    //----------------------------------------------------------------------------------------------------------------------------------------------//
    // create solution, rhs, and residual vectors
    AMP::Discretization::DOFManager::shared_ptr nodalScalarDOF =
        AMP::Discretization::simpleDOFManager::create( meshAdapter, AMP::Mesh::Vertex, 1, 1, true );
    AMP::LinearAlgebra::Vector::shared_ptr solVec =
        AMP::LinearAlgebra::createVector( nodalScalarDOF, cVar, true );
    AMP::LinearAlgebra::Vector::shared_ptr rhsVec =
        AMP::LinearAlgebra::createVector( nodalScalarDOF, fsOutVar, true );
    AMP::LinearAlgebra::Vector::shared_ptr resVec =
        AMP::LinearAlgebra::createVector( nodalScalarDOF, fsOutVar, true );

    //----------------------------------------------------------------------------------------------------------------------------------------------//
    // create parameters for reset test and reset fick and soret operators
    AMP::LinearAlgebra::Vector::shared_ptr tVec =
        AMP::LinearAlgebra::createVector( nodalScalarDOF, tVar, true );
    tVec->setToScalar( 300. );

    fickOp->setVector( 0, tVec );
    soretOp->setVector( 0, tVec );

#ifdef USE_EXT_SILO
    //----------------------------------------------------------------------------------------------------------------------------------------------//
    // register some variables for plotting
    AMP::Utilities::Writer::shared_ptr siloWriter = AMP::Utilities::Writer::buildWriter( "Silo" );
    siloWriter->registerVector( solVec, meshAdapter, AMP::Mesh::Vertex, "Solution" );
    siloWriter->registerVector( resVec, meshAdapter, AMP::Mesh::Vertex, "Residual" );
#endif

    //----------------------------------------------------------------------------------------------------------------------------------------------//
    // Initial guess

    solVec->setToScalar( .05 );
    double initialGuessNorm = solVec->L2Norm();
    std::cout << "initial guess norm = " << initialGuessNorm << "\n";

    nlinBVPOp->modifyInitialSolutionVector( solVec );

    initialGuessNorm = solVec->L2Norm();
    std::cout << "initial guess norm  after apply = " << initialGuessNorm << "\n";

    rhsVec->setToScalar( 0.0 );
    nlinBVPOp->modifyRHSvector( rhsVec );

    //----------------------------------------------------------------------------------------------------------------------------------------------/

    AMP::shared_ptr<AMP::Database> nonlinearSolver_db = input_db->getDatabase( "NonlinearSolver" );
    AMP::shared_ptr<AMP::Database> linearSolver_db =
        nonlinearSolver_db->getDatabase( "LinearSolver" );

    //----------------------------------------------------------------------------------------------------------------------------------------------//
    // initialize the nonlinear solver
    AMP::shared_ptr<AMP::Solver::PetscSNESSolverParameters> nonlinearSolverParams(
        new AMP::Solver::PetscSNESSolverParameters( nonlinearSolver_db ) );

    // change the next line to get the correct communicator out
    nonlinearSolverParams->d_comm          = globalComm;
    nonlinearSolverParams->d_pOperator     = nlinBVPOp;
    nonlinearSolverParams->d_pInitialGuess = solVec;

    AMP::shared_ptr<AMP::Solver::PetscSNESSolver> nonlinearSolver(
        new AMP::Solver::PetscSNESSolver( nonlinearSolverParams ) );

    //----------------------------------------------------------------------------------------------------------------------------------------------//
    AMP::shared_ptr<AMP::Database> fickPreconditioner_db =
        linearSolver_db->getDatabase( "Preconditioner" );
    AMP::shared_ptr<AMP::Solver::SolverStrategyParameters> fickPreconditionerParams(
        new AMP::Solver::SolverStrategyParameters( fickPreconditioner_db ) );
    fickPreconditionerParams->d_pOperator = linBVPOp;
    AMP::shared_ptr<AMP::Solver::TrilinosMLSolver> linearFickPreconditioner(
        new AMP::Solver::TrilinosMLSolver( fickPreconditionerParams ) );

    //----------------------------------------------------------------------------------------------------------------------------------------------//
    // register the preconditioner with the Jacobian free Krylov solver
    AMP::shared_ptr<AMP::Solver::PetscKrylovSolver> linearSolver =
        nonlinearSolver->getKrylovSolver();

    linearSolver->setPreconditioner( linearFickPreconditioner );

    nlinBVPOp->residual( rhsVec, solVec, resVec );
    double initialResidualNorm = resVec->L2Norm();

    AMP::pout << "Initial Residual Norm: " << initialResidualNorm << std::endl;

    nonlinearSolver->setZeroInitialGuess( false );

    nonlinearSolver->solve( rhsVec, solVec );

    nlinBVPOp->residual( rhsVec, solVec, resVec );

    double finalResidualNorm = resVec->L2Norm();

    std::cout << "Final Residual Norm: " << finalResidualNorm << std::endl;

    solVec->makeConsistent( AMP::LinearAlgebra::Vector::CONSISTENT_SET );
    resVec->makeConsistent( AMP::LinearAlgebra::Vector::CONSISTENT_SET );

#ifdef USE_EXT_SILO
    siloWriter->writeFile( exeName, 0 );
#endif

    // store result
    {
        AMP::Mesh::MeshIterator iterator = meshAdapter->getIterator( AMP::Mesh::Vertex, 0 );
        size_t numNodes                  = iterator.size();
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


int main( int argc, char *argv[] )
{
    AMP::AMPManager::startup( argc, argv );
    AMP::UnitTest ut;

    try {
        std::vector<double> fickOnly, fickSoretOff, fickSoretZero, fickOnlyReal, fickSoretOffReal;

        fickTest( &ut, "testPetscSNESSolver-NonlinearFick-cylinder-TUI-1", fickOnly );
        fickSoretTest( &ut, "testPetscSNESSolver-NonlinearFickSoret-cylinder-TUI-1", fickSoretOff );
        fickSoretTest(
            &ut, "testPetscSNESSolver-NonlinearFickSoret-cylinder-TUI-2", fickSoretZero );
        fickTest( &ut, "testPetscSNESSolver-NonlinearFick-cylinder-TUI-2", fickOnlyReal );
        fickSoretTest(
            &ut, "testPetscSNESSolver-NonlinearFickSoret-cylinder-TUI-3", fickSoretOffReal );
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

        if ( l2err1 < 1.e-6 and l2err2 < 1.e-6 and l2err3 < 1.e-6 ) {
            ut.passes( "fick, fick-soret/off, and fick-soret/zero all agree" );
        }
    } catch ( std::exception &err ) {
        std::cout << "ERROR: While testing " << argv[0] << err.what() << std::endl;
        ut.failure( "ERROR: While testing" );
    } catch ( ... ) {
        std::cout << "ERROR: While testing " << argv[0] << "An unknown exception was thrown."
                  << std::endl;
        ut.failure( "ERROR: While testing" );
    }

    ut.report();

    int num_failed = ut.NumFailGlobal();
    AMP::AMPManager::shutdown();
    return num_failed;
}
