#include "AMP/AMP_TPLs.h"
#include "AMP/IO/PIO.h"
#include "AMP/discretization/DOF_Manager.h"
#include "AMP/discretization/simpleDOF_Manager.h"
#include "AMP/matrices/CSRMatrix.h"
#include "AMP/matrices/CSRMatrixParameters.h"
#include "AMP/matrices/CSRPolicy.h"
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
#include "AMP/solvers/testHelpers/SolverTestParameters.h"
#include "AMP/utils/AMPManager.h"
#include "AMP/utils/Database.h"
#include "AMP/utils/UnitTest.h"
#include "AMP/utils/Utilities.h"
#include "AMP/vectors/Variable.h"
#include "AMP/vectors/Vector.h"
#include "AMP/vectors/VectorBuilder.h"

#ifdef AMP_USE_HYPRE
    #include "AMP/matrices/data/hypre/HypreCSRPolicy.h"
#endif

#include <iomanip>
#include <memory>
#include <string>

#include "reference_solver_solutions.h"

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
    auto SpecificPowerVec = AMP::LinearAlgebra::createVector(
        gaussPointDofMap, SpecificPowerVar, true, neutronicsOperator->getMemoryLocation() );

    AMP::LinearAlgebra::Vector::shared_ptr nullVec;
    neutronicsOperator->apply( nullVec, SpecificPowerVec );

    // Integrate Nuclear Source over Desnity * Volume
    AMP_INSIST( input_db->keyExists( "VolumeIntegralOperator" ), "key missing!" );
    auto sourceOperator = std::dynamic_pointer_cast<AMP::Operator::VolumeIntegralOperator>(
        AMP::Operator::OperatorBuilder::createOperator(
            meshAdapter, "VolumeIntegralOperator", input_db ) );

    // Create the power (heat source) vector.
    auto PowerInWattsVar = sourceOperator->getOutputVariable();
    auto PowerInWattsVec = AMP::LinearAlgebra::createVector(
        nodalDofMap, PowerInWattsVar, true, sourceOperator->getMemoryLocation() );
    PowerInWattsVec->zero();

    // convert the vector of specific power to power for a given basis.
    sourceOperator->apply( SpecificPowerVec, PowerInWattsVec );

    // CREATE THE THERMAL BVP OPERATOR
    auto linearOperator = AMP::Operator::OperatorBuilder::createOperator(
        meshAdapter, "DiffusionBVPOperator", input_db );

    auto diffusionOperator =
        std::dynamic_pointer_cast<AMP::Operator::LinearBVPOperator>( linearOperator );

    auto TemperatureInKelvinVec =
        AMP::LinearAlgebra::createVector( nodalDofMap,
                                          diffusionOperator->getInputVariable(),
                                          true,
                                          diffusionOperator->getMemoryLocation() );
    auto RightHandSideVec =
        AMP::LinearAlgebra::createVector( nodalDofMap,
                                          diffusionOperator->getOutputVariable(),
                                          true,
                                          diffusionOperator->getMemoryLocation() );
    auto ResidualVec = AMP::LinearAlgebra::createVector( nodalDofMap,
                                                         diffusionOperator->getOutputVariable(),
                                                         true,
                                                         diffusionOperator->getMemoryLocation() );

    RightHandSideVec->setToScalar( 0.0 );

    // Add the boundary conditions corrections
    auto boundaryOpCorrectionVec =
        AMP::LinearAlgebra::createVector( nodalDofMap,
                                          diffusionOperator->getOutputVariable(),
                                          true,
                                          diffusionOperator->getMemoryLocation() );

    auto boundaryOp = diffusionOperator->getBoundaryOperator();
    boundaryOp->addRHScorrection( boundaryOpCorrectionVec );

    RightHandSideVec->subtract( *PowerInWattsVec, *boundaryOpCorrectionVec );
    RightHandSideVec->makeConsistent();

#if defined( AMP_USE_HYPRE )
    using Policy = AMP::LinearAlgebra::HypreCSRPolicy;
#else
    using Policy = AMP::LinearAlgebra::CSRPolicy<size_t, int, double>;
#endif
    using gidx_t   = typename Policy::gidx_t;
    using lidx_t   = typename Policy::lidx_t;
    using scalar_t = typename Policy::scalar_t;

    gidx_t firstRow, endRow;
    lidx_t nnz_pad;
    std::vector<lidx_t> nnz_d, nnz_od;
    std::vector<lidx_t> rowstart_d, rowstart_od;
    std::vector<gidx_t> cols_d, cols_od;
    std::vector<lidx_t> cols_loc_d, cols_loc_od;
    std::vector<scalar_t> coeffs_d, coeffs_od;

    AMP::LinearAlgebra::transformDofToCSR<Policy>( diffusionOperator->getMatrix(),
                                                   firstRow,
                                                   endRow,
                                                   nnz_d,
                                                   rowstart_d,
                                                   cols_d,
                                                   cols_loc_d,
                                                   coeffs_d,
                                                   nnz_od,
                                                   rowstart_od,
                                                   cols_od,
                                                   cols_loc_od,
                                                   coeffs_od,
                                                   nnz_pad );

    AMP::LinearAlgebra::CSRMatrixParameters<Policy>::CSRSerialMatrixParameters pars_d{
        nnz_d.data(), rowstart_d.data(), cols_d.data(), cols_loc_d.data(), coeffs_d.data()
    };

    AMP::LinearAlgebra::CSRMatrixParameters<Policy>::CSRSerialMatrixParameters pars_od{
        nnz_od.data(), rowstart_od.data(), cols_od.data(), cols_loc_od.data(), coeffs_od.data()
    };

    auto csrParams = std::make_shared<AMP::LinearAlgebra::CSRMatrixParameters<Policy>>(
        firstRow, endRow, pars_d, pars_od, nnz_pad, meshAdapter->getComm() );

#ifdef USE_DEVICE
    using Alloc = AMP::ManagedAllocator<int>;
#else
    using Alloc  = AMP::HostAllocator<int>;
#endif

    auto csrMatrix = std::make_shared<AMP::LinearAlgebra::CSRMatrix<Policy, Alloc>>( csrParams );
    AMP_ASSERT( csrMatrix );

    auto csrOpParams = std::make_shared<AMP::Operator::OperatorParameters>( input_db );
    auto csrOperator = std::make_shared<AMP::Operator::LinearOperator>( csrOpParams );
    csrOperator->setMatrix( csrMatrix );
    csrOperator->setVariables( linearOperator->getInputVariable(),
                               linearOperator->getOutputVariable() );

    // make sure the database on theinput file exists for the linear solver
    AMP_INSIST( input_db->keyExists( "LinearSolver" ), "Key ''LinearSolver'' is missing!" );
    auto comm = AMP::AMP_MPI( AMP_COMM_WORLD );
    auto linearSolver =
        AMP::Solver::Test::buildSolver( "LinearSolver", input_db, comm, nullptr, csrOperator );

    // Set initial guess
    TemperatureInKelvinVec->setToScalar( 1.0 );
    TemperatureInKelvinVec->makeConsistent( AMP::LinearAlgebra::ScatterType::CONSISTENT_SET );

    // Check the initial L2 norm of the solution
    double initSolNorm = static_cast<double>( TemperatureInKelvinVec->L2Norm() );
    std::cout << "Initial Solution Norm: " << initSolNorm << std::endl;
    std::cout << "RHS Norm: " << RightHandSideVec->L2Norm() << std::endl;

    // Use a random initial guess?
    linearSolver->setZeroInitialGuess( false );

    auto b     = RightHandSideVec;
    auto bView = AMP::LinearAlgebra::createVectorAdaptor(
        b->getName(), b->getDOFManager(), b->getRawDataBlock<double>() );
    auto x     = TemperatureInKelvinVec;
    auto xView = AMP::LinearAlgebra::createVectorAdaptor(
        x->getName(), x->getDOFManager(), x->getRawDataBlock<double>() );

    // Solve the problem.
    linearSolver->apply( bView, xView );

    // commented till necessary infrastructure in place
    checkConvergence( linearSolver.get(), inputFileName, *ut );

    input_db.reset();
}


int main( int argc, char *argv[] )
{
    AMP::AMPManager::startup( argc, argv );
    AMP::UnitTest ut;
    AMP::Solver::registerSolverFactories();

    std::vector<std::string> files;
    PROFILE_ENABLE();

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
        files.emplace_back( "input_testLinearSolvers-LinearThermalRobin-BoomerAMG-TFQMR" );
        files.emplace_back( "input_testLinearSolvers-LinearThermalRobin-BoomerAMG-HypreCG" );
        files.emplace_back( "input_testLinearSolvers-LinearThermalRobin-DiagonalPC-HypreCG" );
        files.emplace_back( "input_testLinearSolvers-LinearThermalRobin-HypreCG" );
    #ifdef AMP_USE_PETSC
            // files.emplace_back(
            // "input_testLinearSolvers-LinearThermalRobin-BoomerAMG-PetscFGMRES" );
    #endif
#endif

#ifdef AMP_USE_TRILINOS_ML
        // files.emplace_back( "input_testLinearSolvers-LinearThermalRobin-ML" );
        // files.emplace_back( "input_testLinearSolvers-LinearThermalRobin-ML-CG" );
        // files.emplace_back( "input_testLinearSolvers-LinearThermalRobin-ML-GMRES" );
        // files.emplace_back( "input_testLinearSolvers-LinearThermalRobin-ML-FGMRES" );
        // files.emplace_back( "input_testLinearSolvers-LinearThermalRobin-ML-BiCGSTAB" );
    // files.emplace_back( "input_testLinearSolvers-LinearThermalRobin-ML-TFQMR" );
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

    // build unique profile name to avoid collisions
    std::ostringstream ss;
    ss << "testLinearSolvers-LinearThermalRobinCSR_r" << std::setw( 3 ) << std::setfill( '0' )
       << AMP::AMPManager::getCommWorld().getSize();
    PROFILE_SAVE( ss.str() );

    int num_failed = ut.NumFailGlobal();
    AMP::AMPManager::shutdown();
    return num_failed;
}
