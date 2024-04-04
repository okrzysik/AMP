#include "AMP/AMP_TPLs.h"
#include "AMP/IO/PIO.h"
#include "AMP/discretization/DOF_Manager.h"
#include "AMP/discretization/simpleDOF_Manager.h"
#include "AMP/matrices/CSRMatrix.h"
#include "AMP/matrices/CSRMatrixParameters.h"
#include "AMP/matrices/data/hypre/HypreCSRPolicy.h"
#include "AMP/matrices/testHelpers/MatrixDataTransforms.h"
#include "AMP/mesh/Mesh.h"
#include "AMP/mesh/MeshFactory.h"
#include "AMP/mesh/MeshParameters.h"
#include "AMP/operators/ElementOperationFactory.h"
#include "AMP/operators/ElementPhysicsModelFactory.h"
#include "AMP/operators/LinearBVPOperator.h"
#include "AMP/operators/NeutronicsRhs.h"
#include "AMP/operators/OperatorBuilder.h"
#include "AMP/operators/boundary/DirichletMatrixCorrection.h"
#include "AMP/operators/boundary/DirichletVectorCorrection.h"
#include "AMP/operators/boundary/libmesh/NeumannVectorCorrection.h"
#include "AMP/operators/boundary/libmesh/RobinMatrixCorrection.h"
#include "AMP/operators/diffusion/DiffusionLinearElement.h"
#include "AMP/operators/diffusion/DiffusionLinearFEOperator.h"
#include "AMP/operators/diffusion/DiffusionTransportModel.h"
#include "AMP/operators/libmesh/VolumeIntegralOperator.h"
#include "AMP/solvers/SolverFactory.h"
#include "AMP/solvers/SolverStrategy.h"
#include "AMP/solvers/SolverStrategyParameters.h"
#include "AMP/utils/AMPManager.h"
#include "AMP/utils/Database.h"
#include "AMP/utils/UnitTest.h"
#include "AMP/vectors/Variable.h"
#include "AMP/vectors/Vector.h"
#include "AMP/vectors/VectorBuilder.h"

#ifdef USE_CUDA
    #include "AMP/utils/cuda/CudaAllocator.h"
#endif

#include <iomanip>
#include <memory>
#include <string>

#ifdef USE_CUDA
    #include <cuda_runtime_api.h>
#endif

#include "reference_solver_solutions.h"

std::shared_ptr<AMP::Solver::SolverStrategy>
buildSolver( std::shared_ptr<AMP::Database> input_db,
             const std::string &solver_name,
             const AMP::AMP_MPI &comm,
             std::shared_ptr<AMP::Operator::Operator> op )
{

    AMP_INSIST( input_db->keyExists( solver_name ), "Key " + solver_name + " is missing!" );

    auto db = input_db->getDatabase( solver_name );
    AMP_INSIST( db->keyExists( "name" ), "Key name does not exist in solver database" );

    auto parameters         = std::make_shared<AMP::Solver::SolverStrategyParameters>( db );
    parameters->d_pOperator = op;
    parameters->d_comm      = comm;

    std::shared_ptr<AMP::Solver::SolverStrategy> nestedSolver;

    // check if we need to construct a preconditioner
    auto uses_preconditioner = db->getWithDefault<bool>( "uses_preconditioner", false );
    if ( uses_preconditioner ) {
        auto pc_name = db->getWithDefault<std::string>( "pc_name", "Preconditioner" );
        nestedSolver = buildSolver( input_db, pc_name, comm, op );
        AMP_INSIST( nestedSolver, "null preconditioner" );
    }

    parameters->d_pNestedSolver = nestedSolver;

    return AMP::Solver::SolverFactory::create( parameters );
}

namespace AMP::LinearAlgebra {
std::shared_ptr<AMP::LinearAlgebra::Vector>
createVectorInSpace( std::shared_ptr<AMP::Discretization::DOFManager> DOFs,
                     std::shared_ptr<AMP::LinearAlgebra::Variable> var )
{
#ifdef USE_CUDA
    // We are ready to create a single vector
    // Create the communication list
    AMP_MPI comm = DOFs->getComm();
    AMP_ASSERT( !comm.isNull() );
    comm.barrier();
    std::shared_ptr<CommunicationList> comm_list;
    auto remote_DOFs = DOFs->getRemoteDOFs();
    bool ghosts      = comm.anyReduce( !remote_DOFs.empty() );
    if ( !ghosts ) {
        // No need for a communication list
        comm_list = std::make_shared<CommunicationList>( DOFs->numLocalDOF(), DOFs->getComm() );
    } else {
        // Construct the communication list
        auto params           = std::make_shared<CommunicationListParameters>();
        params->d_comm        = comm;
        params->d_localsize   = DOFs->numLocalDOF();
        params->d_remote_DOFs = remote_DOFs;
        comm_list             = std::make_shared<CommunicationList>( params );
    }
    comm.barrier();

    return AMP::LinearAlgebra::createSimpleVector<
        double,
        AMP::LinearAlgebra::VectorOperationsDefault<double>,
        AMP::LinearAlgebra::VectorDataDefault<double, AMP::CudaManagedAllocator<double>>>(
        var, DOFs, comm_list );
#else
    return AMP::LinearAlgebra::createVector( DOFs, var );
#endif
}
} // namespace AMP::LinearAlgebra

void linearThermalTest( AMP::UnitTest *ut, const std::string &inputFileName )
{
    // Input and output file names
    std::string input_file = inputFileName;
    std::string log_file   = "output_" + inputFileName;

    AMP::pout << "Running linearThermalTest with input " << input_file << std::endl;

    // Fill the database from the input file.
    auto input_db = AMP::Database::parseInputFile( input_file );
    input_db->print( AMP::plog );

    // Print from all cores into the output files
    AMP::logAllNodes( log_file );

    // Create the Mesh
    AMP_INSIST( input_db->keyExists( "Mesh" ), "Key ''Mesh'' is missing!" );
    auto mesh_db   = input_db->getDatabase( "Mesh" );
    auto mgrParams = std::make_shared<AMP::Mesh::MeshParameters>( mesh_db );
    mgrParams->setComm( AMP::AMP_MPI( AMP_COMM_WORLD ) );
    auto meshAdapter = AMP::Mesh::MeshFactory::create( mgrParams );

    // Create a DOF manager for a nodal vector
    int DOFsPerNode          = 1;
    int DOFsPerElement       = 8;
    int nodalGhostWidth      = 1;
    int gaussPointGhostWidth = 1;
    bool split               = true;
    auto nodalDofMap         = AMP::Discretization::simpleDOFManager::create(
        meshAdapter, AMP::Mesh::GeomType::Vertex, nodalGhostWidth, DOFsPerNode, split );
    auto gaussPointDofMap = AMP::Discretization::simpleDOFManager::create(
        meshAdapter, AMP::Mesh::GeomType::Cell, gaussPointGhostWidth, DOFsPerElement, split );

    // CREATE THE NEUTRONICS SOURCE
    AMP_INSIST( input_db->keyExists( "NeutronicsOperator" ),
                "Key ''NeutronicsOperator'' is missing!" );
    auto neutronicsOp_db = input_db->getDatabase( "NeutronicsOperator" );
    auto neutronicsParams =
        std::make_shared<AMP::Operator::NeutronicsRhsParameters>( neutronicsOp_db );
    auto neutronicsOperator = std::make_shared<AMP::Operator::NeutronicsRhs>( neutronicsParams );

    auto SpecificPowerVar = neutronicsOperator->getOutputVariable();
    auto SpecificPowerVec =
        AMP::LinearAlgebra::createVectorInSpace( gaussPointDofMap, SpecificPowerVar );

    AMP::LinearAlgebra::Vector::shared_ptr nullVec;
    neutronicsOperator->apply( nullVec, SpecificPowerVec );

    // Integrate Nuclear Source over Desnity * Volume
    AMP_INSIST( input_db->keyExists( "VolumeIntegralOperator" ), "key missing!" );
    std::shared_ptr<AMP::Operator::ElementPhysicsModel> stransportModel;
    auto sourceOperator = std::dynamic_pointer_cast<AMP::Operator::VolumeIntegralOperator>(
        AMP::Operator::OperatorBuilder::createOperator(
            meshAdapter, "VolumeIntegralOperator", input_db, stransportModel ) );

    // Create the power (heat source) vector.
    auto PowerInWattsVar = sourceOperator->getOutputVariable();
    auto PowerInWattsVec = AMP::LinearAlgebra::createVectorInSpace( nodalDofMap, PowerInWattsVar );
    PowerInWattsVec->zero();

    // convert the vector of specific power to power for a given basis.
    sourceOperator->apply( SpecificPowerVec, PowerInWattsVec );

    // CREATE THE THERMAL BVP OPERATOR
    std::shared_ptr<AMP::Operator::ElementPhysicsModel> transportModel;
    auto linearOperator = AMP::Operator::OperatorBuilder::createOperator(
        meshAdapter, "DiffusionBVPOperator", input_db, transportModel );

    auto diffusionOperator =
        std::dynamic_pointer_cast<AMP::Operator::LinearBVPOperator>( linearOperator );

    auto TemperatureInKelvinVec = AMP::LinearAlgebra::createVectorInSpace(
        nodalDofMap, diffusionOperator->getInputVariable() );
    auto RightHandSideVec = AMP::LinearAlgebra::createVectorInSpace(
        nodalDofMap, diffusionOperator->getOutputVariable() );
    auto ResidualVec = AMP::LinearAlgebra::createVectorInSpace(
        nodalDofMap, diffusionOperator->getOutputVariable() );

    RightHandSideVec->setToScalar( 0.0 );

    // Add the boundary conditions corrections
    auto boundaryOpCorrectionVec = AMP::LinearAlgebra::createVectorInSpace(
        nodalDofMap, diffusionOperator->getOutputVariable() );

    auto boundaryOp = diffusionOperator->getBoundaryOperator();
    boundaryOp->addRHScorrection( boundaryOpCorrectionVec );

    RightHandSideVec->subtract( *PowerInWattsVec, *boundaryOpCorrectionVec );
    RightHandSideVec->makeConsistent();

    // std::cout << "RHS Norm after BC Correction " << RightHandSideVec->L2Norm() << std::endl;
    // std::cout << "RHS Norm 1: " << RightHandSideVec->L2Norm() << std::endl;
    // std::cout << "RHS Norm 2: " << PowerInWattsVec->L2Norm() << std::endl;
    // std::cout << "RHS Norm 3: " << boundaryOpCorrectionVec->L2Norm() << std::endl;

    using Policy   = AMP::LinearAlgebra::HypreCSRPolicy;
    using gidx_t   = typename Policy::gidx_t;
    using lidx_t   = typename Policy::lidx_t;
    using scalar_t = typename Policy::scalar_t;

    gidx_t firstRow, endRow;
    std::vector<lidx_t> nnz;
    std::vector<gidx_t> cols;
    std::vector<scalar_t> coeffs;

    AMP::LinearAlgebra::transformDofToCSR<AMP::LinearAlgebra::HypreCSRPolicy>(
        diffusionOperator->getMatrix(), firstRow, endRow, nnz, cols, coeffs );

    lidx_t *nnz_p      = nullptr;
    gidx_t *cols_p     = nullptr;
    scalar_t *coeffs_p = nullptr;

#ifdef USE_CUDA
    cudaMallocManaged( (void **) &nnz_p, sizeof( lidx_t ) * nnz.size() );
    cudaMallocManaged( (void **) &cols_p, sizeof( gidx_t ) * cols.size() );
    cudaMallocManaged( (void **) &coeffs_p, sizeof( scalar_t ) * coeffs.size() );

    std::memcpy( nnz_p, nnz.data(), sizeof( lidx_t ) * nnz.size() );
    std::memcpy( cols_p, cols.data(), sizeof( gidx_t ) * cols.size() );
    std::memcpy( coeffs_p, coeffs.data(), sizeof( scalar_t ) * coeffs.size() );
#else
    nnz_p    = nnz.data();
    cols_p   = cols.data();
    coeffs_p = coeffs.data();
#endif

    auto csrParams = std::make_shared<AMP::LinearAlgebra::CSRMatrixParameters<Policy>>(
        firstRow, endRow, nnz_p, cols_p, coeffs_p, meshAdapter->getComm() );

    auto csrMatrix = std::make_shared<AMP::LinearAlgebra::CSRMatrix<Policy>>( csrParams );
    AMP_ASSERT( csrMatrix );

    auto csrOpParams = std::make_shared<AMP::Operator::OperatorParameters>( input_db );
    auto csrOperator = std::make_shared<AMP::Operator::LinearOperator>( csrOpParams );
    csrOperator->setMatrix( csrMatrix );
    csrOperator->setVariables( linearOperator->getInputVariable(),
                               linearOperator->getOutputVariable() );

    // make sure the database on theinput file exists for the linear solver
    AMP_INSIST( input_db->keyExists( "LinearSolver" ), "Key ''LinearSolver'' is missing!" );
    auto comm         = AMP::AMP_MPI( AMP_COMM_WORLD );
    auto linearSolver = buildSolver( input_db, "LinearSolver", comm, csrOperator );

    // Set initial guess
    TemperatureInKelvinVec->setToScalar( 1.0 );
    TemperatureInKelvinVec->makeConsistent( AMP::LinearAlgebra::ScatterType::CONSISTENT_SET );

    // Check the initial L2 norm of the solution
    double initSolNorm = static_cast<double>( TemperatureInKelvinVec->L2Norm() );
    std::cout << "Initial Solution Norm: " << initSolNorm << std::endl;
    std::cout << "RHS Norm: " << RightHandSideVec->L2Norm() << std::endl;

    // Use a random initial guess?
    linearSolver->setZeroInitialGuess( false );

    // Solve the problem.
    linearSolver->apply( RightHandSideVec, TemperatureInKelvinVec );

    // commented till necessary infrastructure in place
    checkConvergence( linearSolver.get(), inputFileName, *ut );

    input_db.reset();

#ifdef USE_CUDA
    cudaFree( nnz_p );
    cudaFree( cols_p );
    cudaFree( coeffs_p );
#endif
}


int main( int argc, char *argv[] )
{
    AMP::AMPManager::startup( argc, argv );
    AMP::UnitTest ut;
    AMP::Solver::registerSolverFactories();

    std::vector<std::string> files;

    if ( argc > 1 ) {

        files.emplace_back( argv[1] );

    } else {

        files.emplace_back( "input_testLinearSolvers-LinearThermalRobin-CG" );
        files.emplace_back( "input_testLinearSolvers-LinearThermalRobin-GMRES" );
        files.emplace_back( "input_testLinearSolvers-LinearThermalRobin-FGMRES" );
        files.emplace_back( "input_testLinearSolvers-LinearThermalRobin-BiCGSTAB" );
        files.emplace_back( "input_testLinearSolvers-LinearThermalRobin-TFQMR" );

#ifdef AMP_USE_PETSC
        files.emplace_back( "input_testLinearSolvers-LinearThermalRobin-PetscFGMRES" );
#endif

#ifdef AMP_USE_HYPRE
        files.emplace_back( "input_testLinearSolvers-LinearThermalRobin-BoomerAMG" );
        files.emplace_back( "input_testLinearSolvers-LinearThermalRobin-BoomerAMG-CG" );
        files.emplace_back( "input_testLinearSolvers-LinearThermalRobin-BoomerAMG-GMRES" );
        files.emplace_back( "input_testLinearSolvers-LinearThermalRobin-BoomerAMG-FGMRES" );
        files.emplace_back( "input_testLinearSolvers-LinearThermalRobin-BoomerAMG-BiCGSTAB" );
        //        files.emplace_back( "input_testLinearSolvers-LinearThermalRobin-BoomerAMG-TFQMR"
        //        );
        files.emplace_back( "input_testLinearSolvers-LinearThermalRobin-BoomerAMG-HypreCG" );
        files.emplace_back( "input_testLinearSolvers-LinearThermalRobin-DiagonalPC-HypreCG" );
        files.emplace_back( "input_testLinearSolvers-LinearThermalRobin-HypreCG" );
    #ifdef AMP_USE_PETSC
        files.emplace_back( "input_testLinearSolvers-LinearThermalRobin-BoomerAMG-PetscFGMRES" );
    #endif
#endif

#ifdef AMP_USE_TRILINOS_ML
        // files.emplace_back( "input_testLinearSolvers-LinearThermalRobin-ML" );
        // files.emplace_back( "input_testLinearSolvers-LinearThermalRobin-ML-CG" );
        // files.emplace_back( "input_testLinearSolvers-LinearThermalRobin-ML-GMRES" );
        // files.emplace_back( "input_testLinearSolvers-LinearThermalRobin-ML-FGMRES" );
        // files.emplace_back( "input_testLinearSolvers-LinearThermalRobin-ML-BiCGSTAB" );
        //        files.emplace_back( "input_testLinearSolvers-LinearThermalRobin-ML-TFQMR" );

    #ifdef AMP_USE_PETSC
        files.emplace_back( "input_testLinearSolvers-LinearThermalRobin-ML-PetscFGMRES" );
    #endif
#endif

#ifdef AMP_USE_TRILINOS_MUELU
        // files.emplace_back( "input_testLinearSolvers-LinearThermalRobin-MueLu" );
        // files.emplace_back( "input_testLinearSolvers-LinearThermalRobin-MueLu-GMRES" );
        // files.emplace_back( "input_testLinearSolvers-LinearThermalRobin-MueLu-BiCGSTAB" );
        // files.emplace_back( "input_testLinearSolvers-LinearThermalRobin-MueLu-TFQMR" );
#endif
    }

    for ( auto &file : files )
        linearThermalTest( &ut, file );

    ut.report();

    int num_failed = ut.NumFailGlobal();
    AMP::AMPManager::shutdown();
    return num_failed;
}
