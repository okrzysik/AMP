#include "AMP/IO/PIO.h"
#include "AMP/IO/Writer.h"
#include "AMP/discretization/DOF_Manager.h"
#include "AMP/discretization/simpleDOF_Manager.h"
#include "AMP/mesh/Mesh.h"
#include "AMP/mesh/MeshFactory.h"
#include "AMP/mesh/MeshParameters.h"
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
#include "AMP/operators/mechanics/MechanicsLinearFEOperator.h"
#include "AMP/operators/mechanics/MechanicsNonlinearFEOperator.h"
#include "AMP/solvers/SolverFactory.h"
#include "AMP/solvers/SolverStrategy.h"
#include "AMP/solvers/SolverStrategyParameters.h"
#include "AMP/utils/AMPManager.h"
#include "AMP/utils/AMP_MPI.h"
#include "AMP/utils/Database.h"
#include "AMP/utils/UnitTest.h"
#include "AMP/vectors/VectorBuilder.h"
#include <memory>

#include <iostream>
#include <string>


struct null_deleter {
    void operator()( void const * ) const {}
};


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

void fickSoretTest( AMP::UnitTest *ut, std::string fileName )
{
    std::string input_file = fileName;
    std::string log_file   = "output_" + fileName;

    AMP::logOnlyNodeZero( log_file );
    AMP::AMP_MPI globalComm( AMP_COMM_WORLD );

    auto input_db = AMP::Database::parseInputFile( input_file );
    input_db->print( AMP::plog );

    // Create the Mesh
    AMP_INSIST( input_db->keyExists( "Mesh" ), "Key ''Mesh'' is missing!" );
    auto mesh_db   = input_db->getDatabase( "Mesh" );
    auto mgrParams = std::make_shared<AMP::Mesh::MeshParameters>( mesh_db );
    mgrParams->setComm( AMP::AMP_MPI( AMP_COMM_WORLD ) );
    auto meshAdapter = AMP::Mesh::MeshFactory::create( mgrParams );

    // Create a DOF manager for a nodal vector
    int DOFsPerNode     = 1;
    int nodalGhostWidth = 1;
    bool split          = true;
    auto nodalDofMap    = AMP::Discretization::simpleDOFManager::create(
        meshAdapter, AMP::Mesh::GeomType::Vertex, nodalGhostWidth, DOFsPerNode, split );

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

    auto tVar = std::make_shared<AMP::LinearAlgebra::Variable>( "temperature" );
    auto cVar( fickOp->getOutputVariable() );
    auto fsOutVar = nlinBVPOp->getOutputVariable();

    // create solution, rhs, and residual vectors
    auto solVec = AMP::LinearAlgebra::createVector( nodalDofMap, cVar );
    auto rhsVec = AMP::LinearAlgebra::createVector( nodalDofMap, fsOutVar );
    auto resVec = AMP::LinearAlgebra::createVector( nodalDofMap, fsOutVar );

    // create parameters for reset test and reset fick and soret operators
    auto tVec = AMP::LinearAlgebra::createVector( nodalDofMap, tVar );

    fickOp->setVector( "temperature", tVec );
    soretOp->setVector( "temperature", tVec );

    auto fickFrozen  = fickOp->getFrozen();
    auto soretFrozen = soretOp->getFrozen();

    double lenscale = input_db->getScalar<double>( "LengthScale" );
    soretFrozen["temperature"]->setToScalar( 300. );
    auto iterator = meshAdapter->getIterator( AMP::Mesh::GeomType::Vertex, 0 );
    for ( ; iterator != iterator.end(); ++iterator ) {
        double x, y;
        x = ( iterator->coord() )[0];
        y = ( iterator->coord() )[1];
        std::vector<size_t> gid;
        nodalDofMap->getDOFs( iterator->globalID(), gid );
        double value =
            300. + 450 * ( 1. - ( x * x / lenscale / lenscale + y * y / lenscale / lenscale ) );
        fickFrozen["temperature"]->setValueByGlobalID( gid[0], value );
        soretFrozen["temperature"]->setValueByGlobalID( gid[0], value );
    }

    // Initial guess
    auto initialValue = input_db->getScalar<double>( "InitialValue" );
    solVec->setToScalar( initialValue );
    auto initialGuessNorm = static_cast<double>( solVec->L2Norm() );
    std::cout << "initial guess norm = " << initialGuessNorm << "\n";

    nlinBVPOp->modifyInitialSolutionVector( solVec );

    initialGuessNorm = static_cast<double>( solVec->L2Norm() );
    std::cout << "initial guess norm  after apply = " << initialGuessNorm << "\n";

    rhsVec->setToScalar( 0.0 );
    nlinBVPOp->modifyRHSvector( rhsVec );

    // Create the solver
    auto nonlinearSolver =
        buildSolver( "NonlinearSolver", input_db, globalComm, solVec, nlinBVPOp );

    nlinBVPOp->residual( rhsVec, solVec, resVec );
    double initialResidualNorm = static_cast<double>( resVec->L2Norm() );

    AMP::pout << "Initial Residual Norm: " << initialResidualNorm << std::endl;

    nonlinearSolver->setZeroInitialGuess( false );

    nonlinearSolver->apply( rhsVec, solVec );

    nlinBVPOp->residual( rhsVec, solVec, resVec );

    double finalResidualNorm = static_cast<double>( resVec->L2Norm() );

    std::cout << "Final Residual Norm: " << finalResidualNorm << std::endl;

    solVec->makeConsistent( AMP::LinearAlgebra::ScatterType::CONSISTENT_SET );

    // evaluate and register material coefficients for graphical output
    auto fickCoeffVar  = std::make_shared<AMP::LinearAlgebra::Variable>( "FickCoefficient" );
    auto soretCoeffVar = std::make_shared<AMP::LinearAlgebra::Variable>( "SoretCoefficient" );
    auto fickCoeffVec  = AMP::LinearAlgebra::createVector( nodalDofMap, fickCoeffVar );
    auto soretCoeffVec = AMP::LinearAlgebra::createVector( nodalDofMap, soretCoeffVar );
    auto fickModel     = fickOp->getTransportModel();
    auto soretModel    = soretOp->getTransportModel();

    {
        iterator      = meshAdapter->getIterator( AMP::Mesh::GeomType::Vertex, 0 );
        size_t nnodes = fickCoeffVec->getLocalSize(), node;
        std::vector<size_t> gids( nnodes );
        std::vector<double> temp( nnodes ), conc( nnodes ), fickCoeff( nnodes ),
            soretCoeff( nnodes ), burn( nnodes );
        for ( node = 0; iterator != iterator.end(); iterator++ ) {
            std::vector<size_t> gid;
            nodalDofMap->getDOFs( iterator->globalID(), gid );
            gids[node] = gid[0];
            node++;
        }
        AMP_INSIST( node == nnodes, "invalid count" );
        fickFrozen["temperature"]->getValuesByGlobalID( nnodes, &gids[0], &temp[0] );
        solVec->getValuesByGlobalID( nnodes, &gids[0], &conc[0] );
        // This is kevin - i found out because the vector wasn't used when silo is not enabled.
        std::map<std::string, std::shared_ptr<std::vector<double>>> args;
        args.insert( std::make_pair(
            "temperature", std::shared_ptr<std::vector<double>>( &temp, null_deleter() ) ) );
        args.insert( std::make_pair(
            "concentration", std::shared_ptr<std::vector<double>>( &conc, null_deleter() ) ) );
        args.insert( std::make_pair(
            "burnup", std::shared_ptr<std::vector<double>>( &burn, null_deleter() ) ) );
        fickModel->getTransport( fickCoeff, args );
        soretModel->getTransport( soretCoeff, args );
        fickCoeffVec->setValuesByGlobalID( nnodes, &gids[0], &fickCoeff[0] );
        soretCoeffVec->setValuesByGlobalID( nnodes, &gids[0], &soretCoeff[0] );
    }

#ifdef AMP_USE_SILO
    // write graphical output
    auto siloWriter = AMP::IO::Writer::buildWriter( "Silo" );
    siloWriter->registerMesh( meshAdapter );
    siloWriter->registerVector( solVec, meshAdapter, AMP::Mesh::GeomType::Vertex, "Solution" );
    siloWriter->registerVector( resVec, meshAdapter, AMP::Mesh::GeomType::Vertex, "Residual" );
    siloWriter->registerVector(
        fickFrozen["temperature"], meshAdapter, AMP::Mesh::GeomType::Vertex, "Temperature" );
    siloWriter->registerVector(
        fickCoeffVec, meshAdapter, AMP::Mesh::GeomType::Vertex, "FickCoefficient" );
    siloWriter->registerVector(
        soretCoeffVec, meshAdapter, AMP::Mesh::GeomType::Vertex, "ThermalDiffusionCoefficient" );
    siloWriter->writeFile( fileName, 0 );
#endif

    if ( finalResidualNorm > 1.0e-08 ) {
        ut->failure( fileName );
    } else {
        ut->passes( "PetscSNES Solver successfully solves a nonlinear Fick-Soret equation with "
                    "Jacobian provided, "
                    "FGMRES for Krylov" );
    }
    ut->passes( fileName );
}


int main( int argc, char *argv[] )
{
    AMP::AMPManager::startup( argc, argv );
    AMP::UnitTest ut;
    AMP::Solver::registerSolverFactories();

    std::vector<std::string> inputNames;

    if ( argc > 1 ) {
        inputNames.push_back( argv[1] );
    } else {
#ifdef AMP_USE_PETSC
    #ifdef AMP_USE_HYPRE
        inputNames.emplace_back(
            "input_testPetscSNESSolver-NonlinearFickSoret-cylinder-OxMSRZC09-1" );
    #endif
#else
        AMP_ERROR( "Test requires both PETSC and HYPRE at present" );
#endif
    }

    for ( auto &inputName : inputNames )
        fickSoretTest( &ut, inputName );

    ut.report();

    int num_failed = ut.NumFailGlobal();
    AMP::AMPManager::shutdown();
    return num_failed;
}
