#include "AMP/AMP_TPLs.h"
#include "AMP/IO/PIO.h"
#include "AMP/discretization/DOF_Manager.h"
#include "AMP/discretization/simpleDOF_Manager.h"
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
#include "AMP/vectors/Variable.h"
#include "AMP/vectors/Vector.h"
#include "AMP/vectors/VectorBuilder.h"

#include <memory>
#include <string>

#include "reference_solver_solutions_hypre.h"


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
    auto comm      = AMP::AMP_MPI( AMP_COMM_WORLD );
    mgrParams->setComm( comm );
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

    auto SpecificPowerVec =
        AMP::LinearAlgebra::createVector( gaussPointDofMap,
                                          neutronicsOperator->getOutputVariable(),
                                          true,
                                          neutronicsOperator->getMemoryLocation() );

    AMP::LinearAlgebra::Vector::shared_ptr nullVec;
    neutronicsOperator->apply( nullVec, SpecificPowerVec );

    // Integrate Nuclear Source over Density * Volume
    AMP_INSIST( input_db->keyExists( "VolumeIntegralOperator" ), "key missing!" );
    auto sourceOperator = std::dynamic_pointer_cast<AMP::Operator::VolumeIntegralOperator>(
        AMP::Operator::OperatorBuilder::createOperator(
            meshAdapter, "VolumeIntegralOperator", input_db ) );

    // Create the power (heat source) vector.
    auto PowerInWattsVec = AMP::LinearAlgebra::createVector( nodalDofMap,
                                                             sourceOperator->getOutputVariable(),
                                                             true,
                                                             sourceOperator->getMemoryLocation() );
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

    // Add the boundary conditions corrections
    auto boundaryOpCorrectionVec =
        AMP::LinearAlgebra::createVector( nodalDofMap,
                                          diffusionOperator->getOutputVariable(),
                                          true,
                                          diffusionOperator->getMemoryLocation() );

    auto boundaryOp = diffusionOperator->getBoundaryOperator();
    boundaryOp->addRHScorrection( boundaryOpCorrectionVec );

    RightHandSideVec->subtract( *PowerInWattsVec, *boundaryOpCorrectionVec );

    auto linearSolver =
        AMP::Solver::Test::buildSolver( "LinearSolver", input_db, comm, nullptr, linearOperator );

    // Set initial guess
    TemperatureInKelvinVec->setToScalar( 1.0 );

    // Check the initial L2 norm of the solution
    double initSolNorm = static_cast<double>( TemperatureInKelvinVec->L2Norm() );
    AMP::pout << "Initial Solution Norm: " << initSolNorm << std::endl;
    AMP::pout << "RHS Norm: " << RightHandSideVec->L2Norm() << std::endl;
    AMP::pout << "System size: " << RightHandSideVec->getGlobalSize() << std::endl;

    // Use a random initial guess?
    linearSolver->setZeroInitialGuess( false );

    // Solve the problem.
    linearSolver->apply( RightHandSideVec, TemperatureInKelvinVec );

    checkConvergence( linearSolver.get(), inputFileName, *ut );
}


int main( int argc, char *argv[] )
{
    AMP::AMPManager::startup( argc, argv );
    AMP::UnitTest ut;

    std::vector<std::string> files;

    if ( argc > 1 ) {

        files.emplace_back( argv[1] );

    } else {

#ifdef AMP_USE_HYPRE
        files.emplace_back( "input_testLinearSolvers-LinearThermalRobin-BoomerAMG" );
        files.emplace_back( "input_testLinearSolvers-LinearThermalRobin-BoomerAMG-CG" );
        files.emplace_back( "input_testLinearSolvers-LinearThermalRobin-CylMesh-BoomerAMG" );
        files.emplace_back( "input_testLinearSolvers-LinearThermalRobin-CylMesh-BoomerAMG-CG" );
        files.emplace_back( "input_testLinearSolvers-LinearThermalRobin-BoomerAMG-GMRES" );
        files.emplace_back( "input_testLinearSolvers-LinearThermalRobin-BoomerAMG-FGMRES" );
        files.emplace_back( "input_testLinearSolvers-LinearThermalRobin-BoomerAMG-BiCGSTAB" );
        files.emplace_back( "input_testLinearSolvers-LinearThermalRobin-BoomerAMG-TFQMR" );
        files.emplace_back( "input_testLinearSolvers-LinearThermalRobin-BoomerAMG-HypreCG" );
        files.emplace_back( "input_testLinearSolvers-LinearThermalRobin-DiagonalPC-HypreCG" );
        files.emplace_back( "input_testLinearSolvers-LinearThermalRobin-HypreCG" );
    #ifdef AMP_USE_PETSC
        files.emplace_back( "input_testLinearSolvers-LinearThermalRobin-BoomerAMG-PetscFGMRES" );
    #endif
#endif
    }

    for ( auto &file : files )
        linearThermalTest( &ut, file );

    ut.report();

    int num_failed = ut.NumFailGlobal();
    AMP::AMPManager::shutdown();
    return num_failed;
}
