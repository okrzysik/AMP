#include "AMP/materials/Material.h"
#include "AMP/utils/AMPManager.h"
#include "AMP/utils/Database.h"
#include "AMP/utils/InputDatabase.h"
#include "AMP/utils/InputManager.h"
#include "AMP/utils/PIO.h"
#include "AMP/utils/UnitTest.h"
#include "AMP/utils/Utilities.h"
#include "AMP/utils/shared_ptr.h"
#include "AMP/vectors/Variable.h"
#include <string>

#include "AMP/ampmesh/Mesh.h"
#include "AMP/discretization/DOF_Manager.h"
#include "AMP/discretization/simpleDOF_Manager.h"
#include "AMP/vectors/VectorBuilder.h"

#include "AMP/utils/Writer.h"
#include "AMP/vectors/Vector.h"

#include "AMP/operators/libmesh/MassLinearElement.h"
#include "AMP/operators/libmesh/MassLinearFEOperator.h"

#include "AMP/operators/diffusion/DiffusionLinearElement.h"
#include "AMP/operators/diffusion/DiffusionLinearFEOperator.h"
#include "AMP/operators/diffusion/DiffusionNonlinearFEOperator.h"
#include "AMP/operators/diffusion/DiffusionTransportModel.h"

#include "AMP/operators/NeutronicsRhs.h"
#include "AMP/operators/libmesh/VolumeIntegralOperator.h"
#include "AMP/operators/subchannel/FlowFrapconJacobian.h"
#include "AMP/operators/subchannel/FlowFrapconOperator.h"

#include "AMP/operators/ElementOperationFactory.h"
#include "AMP/operators/ElementPhysicsModelFactory.h"

#include "AMP/operators/boundary/ColumnBoundaryOperator.h"
#include "AMP/operators/boundary/DirichletMatrixCorrection.h"
#include "AMP/operators/boundary/DirichletVectorCorrection.h"
#include "AMP/operators/boundary/libmesh/NeumannVectorCorrection.h"
#include "AMP/operators/boundary/libmesh/NeumannVectorCorrectionParameters.h"
#include "AMP/operators/boundary/libmesh/RobinMatrixCorrection.h"
#include "AMP/operators/boundary/libmesh/RobinVectorCorrection.h"

#include "AMP/operators/CoupledOperator.h"
#include "AMP/operators/CoupledOperatorParameters.h"
#include "AMP/operators/NonlinearBVPOperator.h"

#include "AMP/operators/LinearBVPOperator.h"
#include "AMP/operators/OperatorBuilder.h"
#include "AMP/operators/map/Map1Dto3D.h"
#include "AMP/operators/map/Map3Dto1D.h"
#include "AMP/operators/map/MapOperatorParameters.h"

#include "AMP/solvers/ColumnSolver.h"
#include "AMP/solvers/libmesh/Flow1DSolver.h"
#include "AMP/solvers/petsc/PetscKrylovSolver.h"
#include "AMP/solvers/petsc/PetscKrylovSolverParameters.h"
#include "AMP/solvers/petsc/PetscSNESSolver.h"
#include "AMP/solvers/petsc/PetscSNESSolverParameters.h"
#include "AMP/solvers/trilinos/ml/TrilinosMLSolver.h"


void flowTest( AMP::UnitTest *ut, const std::string &exeName )
{
    std::string input_file = "input_" + exeName;
    std::string log_file   = "output_" + exeName;
    AMP::shared_ptr<AMP::InputDatabase> input_db( new AMP::InputDatabase( "input_db" ) );
    AMP::InputManager::getManager()->parseInputFile( input_file, input_db );
    input_db->printClassData( AMP::plog );

    AMP::PIO::logAllNodes( log_file );
    AMP::AMP_MPI globalComm = AMP::AMP_MPI( AMP_COMM_WORLD );

    //--------------------------------------------------
    //   Create the Mesh.
    //--------------------------------------------------
    AMP_INSIST( input_db->keyExists( "Mesh" ), "Key ''Mesh'' is missing!" );
    AMP::shared_ptr<AMP::Database> mesh_db = input_db->getDatabase( "Mesh" );
    AMP::shared_ptr<AMP::Mesh::MeshParameters> mgrParams(
        new AMP::Mesh::MeshParameters( mesh_db ) );
    mgrParams->setComm( AMP::AMP_MPI( AMP_COMM_WORLD ) );
    AMP::shared_ptr<AMP::Mesh::Mesh> meshAdapter = AMP::Mesh::Mesh::buildMesh( mgrParams );
    //--------------------------------------------------

    //--------------------------------------------------
    // Create a DOF manager for a nodal vector
    //--------------------------------------------------
    int DOFsPerNode          = 1;
    int DOFsPerElement       = 8;
    int nodalGhostWidth      = 1;
    int gaussPointGhostWidth = 1;
    bool split               = true;
    AMP::Discretization::DOFManager::shared_ptr nodalDofMap =
        AMP::Discretization::simpleDOFManager::create(
            meshAdapter, AMP::Mesh::GeomType::Vertex, nodalGhostWidth, DOFsPerNode, split );
    AMP::Discretization::DOFManager::shared_ptr gaussPointDofMap =
        AMP::Discretization::simpleDOFManager::create(
            meshAdapter, AMP::Mesh::GeomType::Volume, gaussPointGhostWidth, DOFsPerElement, split );
    //--------------------------------------------------

    AMP::LinearAlgebra::Vector::shared_ptr nullVec;

    double intguess = input_db->getDoubleWithDefault( "InitialGuess", 400 );

    //-----------------------------------------------
    //   CREATE THE NONLINEAR THERMAL OPERATOR 1 ----
    //-----------------------------------------------
    AMP_INSIST( input_db->keyExists( "NonlinearThermalOperator" ), "key missing!" );
    AMP::shared_ptr<AMP::Operator::ElementPhysicsModel> thermalTransportModel;
    AMP::shared_ptr<AMP::Operator::NonlinearBVPOperator> thermalNonlinearOperator =
        AMP::dynamic_pointer_cast<AMP::Operator::NonlinearBVPOperator>(
            AMP::Operator::OperatorBuilder::createOperator(
                meshAdapter, "NonlinearThermalOperator", input_db, thermalTransportModel ) );

    // initialize the input variable
    AMP::shared_ptr<AMP::Operator::DiffusionNonlinearFEOperator> thermalVolumeOperator =
        AMP::dynamic_pointer_cast<AMP::Operator::DiffusionNonlinearFEOperator>(
            thermalNonlinearOperator->getVolumeOperator() );


    auto globalSolVec =
        AMP::LinearAlgebra::createVector( nodalDofMap, thermalVolumeOperator->getOutputVariable() );
    auto globalRhsVec =
        AMP::LinearAlgebra::createVector( nodalDofMap, thermalVolumeOperator->getOutputVariable() );
    auto globalResVec =
        AMP::LinearAlgebra::createVector( nodalDofMap, thermalVolumeOperator->getOutputVariable() );

    globalSolVec->setToScalar( intguess );

    //-------------------------------------
    //   CREATE THE LINEAR THERMAL OPERATOR ----
    //-------------------------------------

    AMP::shared_ptr<AMP::Operator::ElementPhysicsModel> transportModel;
    AMP::shared_ptr<AMP::Operator::LinearBVPOperator> thermalLinearOperator =
        AMP::dynamic_pointer_cast<AMP::Operator::LinearBVPOperator>(
            AMP::Operator::OperatorBuilder::createOperator(
                meshAdapter, "LinearThermalOperator", input_db, transportModel ) );

    //-------------------------------------
    //  CREATE THE NEUTRONICS SOURCE  //
    //-------------------------------------
    AMP_INSIST( input_db->keyExists( "NeutronicsOperator" ),
                "Key ''NeutronicsOperator'' is missing!" );
    AMP::shared_ptr<AMP::Database> neutronicsOp_db = input_db->getDatabase( "NeutronicsOperator" );
    AMP::shared_ptr<AMP::Operator::NeutronicsRhsParameters> neutronicsParams(
        new AMP::Operator::NeutronicsRhsParameters( neutronicsOp_db ) );
    AMP::shared_ptr<AMP::Operator::NeutronicsRhs> neutronicsOperator(
        new AMP::Operator::NeutronicsRhs( neutronicsParams ) );

    AMP::LinearAlgebra::Variable::shared_ptr SpecificPowerVar =
        neutronicsOperator->getOutputVariable();
    AMP::LinearAlgebra::Vector::shared_ptr SpecificPowerVec =
        AMP::LinearAlgebra::createVector( gaussPointDofMap, SpecificPowerVar );

    neutronicsOperator->apply( nullVec, SpecificPowerVec );

    //----------------------------------------------------------
    //  Integrate Nuclear Rhs over Desnity * GeomType::Volume //
    //----------------------------------------------------------

    AMP_INSIST( input_db->keyExists( "VolumeIntegralOperator" ), "key missing!" );

    AMP::shared_ptr<AMP::Operator::ElementPhysicsModel> stransportModel;
    AMP::shared_ptr<AMP::Operator::VolumeIntegralOperator> sourceOperator =
        AMP::dynamic_pointer_cast<AMP::Operator::VolumeIntegralOperator>(
            AMP::Operator::OperatorBuilder::createOperator(
                meshAdapter, "VolumeIntegralOperator", input_db, stransportModel ) );

    // Create the power (heat source) vector.
    AMP::LinearAlgebra::Variable::shared_ptr PowerInWattsVar = sourceOperator->getOutputVariable();
    AMP::LinearAlgebra::Vector::shared_ptr PowerInWattsVec =
        AMP::LinearAlgebra::createVector( nodalDofMap, PowerInWattsVar );
    PowerInWattsVec->zero();

    // convert the vector of specific power to power for a given basis.
    sourceOperator->apply( SpecificPowerVec, PowerInWattsVec );

    //--------------------------------------
    AMP_INSIST( input_db->keyExists( "NonlinearSolver" ), "Key ''NonlinearSolver'' is missing!" );

    // AMP::shared_ptr<AMP::Database> nonlinearSolver_db1 =
    // input_db->getDatabase("NonlinearSolver");
    // AMP::shared_ptr<AMP::Database>    linearSolver_db1 =
    // nonlinearSolver_db1->getDatabase("LinearSolver");

    ut->passes( "set up to the iterations passes." );


    //-------------------------------------
    AMP::Operator::Operator::shared_ptr boundaryOp =
        thermalNonlinearOperator->getBoundaryOperator();

    //  AMP::shared_ptr<AMP::Operator::RobinVectorCorrection> robinBoundaryOp =
    //  AMP::dynamic_pointer_cast<AMP::Operator::RobinVectorCorrection>(
    //  thermalNonlinearOperator->getBoundaryOperator() );
    //  AMP::shared_ptr<AMP::Operator::NeumannVectorCorrectionParameters> correctionParameters =
    //  AMP::dynamic_pointer_cast<AMP::Operator::NeumannVectorCorrectionParameters>(robinBoundaryOp->getParameters())
    //  ;

    //  robinBoundaryOp->setVariableFlux( robinRHSVec );

    //------------------------------------------------------------------
    AMP::shared_ptr<AMP::Database> nonlinearSolver_db = input_db->getDatabase( "NonlinearSolver" );
    AMP::shared_ptr<AMP::Database> linearSolver_db =
        nonlinearSolver_db->getDatabase( "LinearSolver" );

    //----------------------------------------------------------------//
    // initialize the nonlinear solver
    AMP::shared_ptr<AMP::Solver::PetscSNESSolverParameters> nonlinearSolverParams(
        new AMP::Solver::PetscSNESSolverParameters( nonlinearSolver_db ) );

    // change the next line to get the correct communicator out
    nonlinearSolverParams->d_comm          = globalComm;
    nonlinearSolverParams->d_pOperator     = thermalNonlinearOperator;
    nonlinearSolverParams->d_pInitialGuess = globalSolVec;
    AMP::shared_ptr<AMP::Solver::PetscSNESSolver> nonlinearSolver(
        new AMP::Solver::PetscSNESSolver( nonlinearSolverParams ) );
    //-------------------------------------------------------------------------//
    // initialize the column preconditioner which is a diagonal block preconditioner
    AMP::shared_ptr<AMP::Database> columnPreconditioner_db =
        linearSolver_db->getDatabase( "Preconditioner" );

    AMP::shared_ptr<AMP::Database> thermalPreconditioner_db =
        columnPreconditioner_db->getDatabase( "pelletThermalPreconditioner" );
    AMP::shared_ptr<AMP::Solver::SolverStrategyParameters> thermalPreconditionerParams(
        new AMP::Solver::SolverStrategyParameters( thermalPreconditioner_db ) );
    thermalPreconditionerParams->d_pOperator = thermalLinearOperator;
    AMP::shared_ptr<AMP::Solver::TrilinosMLSolver> thermalPreconditioner(
        new AMP::Solver::TrilinosMLSolver( thermalPreconditionerParams ) );

    //--------------------------------------------------------------------//
    // register the preconditioner with the Jacobian free Krylov solver
    AMP::shared_ptr<AMP::Solver::PetscKrylovSolver> linearSolver =
        nonlinearSolver->getKrylovSolver();
    linearSolver->setPreconditioner( thermalPreconditioner );

    //-------------------------------------
    nonlinearSolver->setZeroInitialGuess( false );


    globalRhsVec->zero();

    globalRhsVec->copyVector( PowerInWattsVec );
    std::cout << "PowerInWattsVec norm  inside loop = " << globalRhsVec->L2Norm() << "\n";
    double expectedVal = 0.175811;
    if ( !AMP::Utilities::approx_equal( expectedVal, globalRhsVec->L2Norm(), 1e-5 ) ) {
        ut->failure( "the PowerInWattsVec norm has changed." );
    }

    //    robinBoundaryOp->reset(correctionParameters);

    thermalNonlinearOperator->modifyRHSvector( globalRhsVec );
    thermalNonlinearOperator->modifyInitialSolutionVector( globalSolVec );

    thermalNonlinearOperator->residual( globalRhsVec, globalSolVec, globalResVec );
    AMP::pout << "Initial Residual Norm for Step is: " << globalResVec->L2Norm() << std::endl;
    expectedVal = 4.84311;
    if ( !AMP::Utilities::approx_equal( expectedVal, globalResVec->L2Norm(), 1e-5 ) ) {
        ut->failure( "the Initial Residual Norm has changed." );
    }

    std::cout << " RHS Vec L2 Norm " << globalRhsVec->L2Norm() << std::endl;
    nonlinearSolver->solve( globalRhsVec, globalSolVec );

    std::cout << "Final Solution Norm: " << globalSolVec->L2Norm() << std::endl;
    expectedVal = 51541;
    if ( !AMP::Utilities::approx_equal( expectedVal, globalSolVec->L2Norm(), 1e-5 ) ) {
        ut->failure( "the Final Solution Norm has changed." );
    }

    thermalNonlinearOperator->residual( globalRhsVec, globalSolVec, globalResVec );
    AMP::pout << "Final   Residual Norm for Step is: " << globalResVec->L2Norm() << std::endl;
    expectedVal = 1. - 10;
    if ( !AMP::Utilities::approx_equal( expectedVal, globalResVec->L2Norm(), 10.0 ) ) {
        ut->failure( "the Final Residual Norm has changed." );
    }

    //---------------------------------------------------------------------------

#ifdef USE_EXT_SILO
    AMP::Utilities::Writer::shared_ptr siloWriter = AMP::Utilities::Writer::buildWriter( "Silo" );
    siloWriter->registerMesh( meshAdapter );

    siloWriter->registerVector(
        globalSolVec, meshAdapter, AMP::Mesh::GeomType::Vertex, "Temperature" );

    siloWriter->writeFile( input_file, 0 );
#endif

    if ( globalResVec->L2Norm() < 10e-6 ) {
        ut->passes( "Seggregated solve of Composite Operator using control loop of "
                    "Thermal+Robin->Map->Flow->Map ." );
    } else {
        ut->failure( "Seggregated solve of Composite Operator using control loop of "
                     "Thermal+Robin->Map->Flow->Map ." );
    }


    //-------------------------------------
    // The 3D-to-1D map is not working in parallel.
    //   -- See Bug 1219 and 1209.
    //} else {
    //  ut.expected_failure("parallel map3D-1D and map1D-3D fail in parallel, see bug #1219.");
    //}
    input_db.reset();

    ut->passes( exeName );
}

int main( int argc, char *argv[] )
{
    AMP::AMPManager::startup( argc, argv );
    AMP::UnitTest ut;

    flowTest( &ut, "testThermalRobinFlow-2" );

    ut.report();

    int num_failed = ut.NumFailGlobal();
    AMP::AMPManager::shutdown();
    return num_failed;
}
