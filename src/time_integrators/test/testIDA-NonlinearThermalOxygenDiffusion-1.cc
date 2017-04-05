#include "materials/Material.h"
#include "operators/NeutronicsRhs.h"
#include "utils/AMPManager.h"
#include "utils/AMP_MPI.h"
#include "utils/Database.h"
#include "utils/InputDatabase.h"
#include "utils/InputManager.h"
#include "utils/PIO.h"
#include "utils/UnitTest.h"
#include "utils/Utilities.h"
#include "utils/Utilities.h"
#include "utils/shared_ptr.h"
#include "vectors/Variable.h"
#include <string>

#include "operators/NeutronicsRhs.h"
#include "operators/NonlinearBVPOperator.h"
#include "operators/libmesh/VolumeIntegralOperator.h"
#include "utils/Writer.h"
#include "vectors/Vector.h"

#include "ampmesh/Mesh.h"
#include "discretization/DOF_Manager.h"
#include "discretization/simpleDOF_Manager.h"
#include "vectors/VectorBuilder.h"

#include "operators/OperatorBuilder.h"
#include "operators/diffusion/DiffusionLinearFEOperator.h"
#include "operators/diffusion/FickSoretNonlinearFEOperator.h"
#include "operators/libmesh/MassLinearFEOperator.h"

#include "operators/ColumnOperator.h"
#include "operators/LinearBVPOperator.h"
#include "operators/boundary/DirichletMatrixCorrection.h"
#include "operators/boundary/DirichletVectorCorrection.h"
#include "operators/mechanics/MechanicsLinearElement.h"
#include "operators/mechanics/MechanicsLinearFEOperator.h"
#include "solvers/ColumnSolver.h"
#include "solvers/trilinos/TrilinosMLSolver.h"

#include "time_integrators/ColumnTimeOperator.h"
#include "time_integrators/sundials/IDATimeIntegrator.h"
#include "time_integrators/sundials/IDATimeOperator.h"

//---------------------------------------------------------------------------//
// TESTS
//---------------------------------------------------------------------------//

#define __PI__ 3.14159265

#define __INIT_THERMAL_FN__( x, y, z, t )                                             \
    ( 750.0 +                                                                         \
      10000.0 * ( 0.5 + x ) * ( 0.5 - x ) * ( 0.5 + y ) * ( 0.5 - y ) * ( 0.5 + z ) * \
          ( 0.5 - z ) )
#define __INIT_OXYGEN_FN__( x, y, z, t ) \
    ( 0.01 + ( 0.5 + x ) * ( 0.5 - x ) * ( 0.5 + y ) * ( 0.5 - y ) * ( 0.5 + z ) * ( 0.5 - z ) )

void IDATimeIntegratorTest( AMP::UnitTest *ut )
{
    std::string input_file = "input_testIDA-NonlinearThermalOxygenDiffusion-1";
    std::string log_file   = "output_testIDA-NonlinearThermalOxygenDiffusion-1";

    AMP::PIO::logOnlyNodeZero( log_file );

    AMP::shared_ptr<AMP::InputDatabase> input_db( new AMP::InputDatabase( "input_db" ) );
    AMP::InputManager::getManager()->parseInputFile( input_file, input_db );
    input_db->printClassData( AMP::plog );

    AMP_INSIST( input_db->keyExists( "Mesh" ), "Key ''Mesh'' is missing!" );
    AMP::shared_ptr<AMP::Database> mesh_db = input_db->getDatabase( "Mesh" );
    AMP::shared_ptr<AMP::Mesh::MeshParameters> mgrParams(
        new AMP::Mesh::MeshParameters( mesh_db ) );
    mgrParams->setComm( AMP::AMP_MPI( AMP_COMM_WORLD ) );
    AMP::shared_ptr<AMP::Mesh::Mesh> meshAdapter = AMP::Mesh::Mesh::buildMesh( mgrParams );

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

    //----------------------------------------------------------------------------------------------------------------------------------------------//
    // create two nonlinear thermal BVP operators
    AMP_INSIST( input_db->keyExists( "NonlinearThermalOperator" ), "key missing!" );
    AMP_INSIST( input_db->keyExists( "NonlinearOxygenOperator" ), "key missing!" );

    AMP::shared_ptr<AMP::Operator::ElementPhysicsModel> thermalModel;
    AMP::shared_ptr<AMP::Operator::NonlinearBVPOperator> nonlinearThermalOperator =
        AMP::dynamic_pointer_cast<AMP::Operator::NonlinearBVPOperator>(
            AMP::Operator::OperatorBuilder::createOperator(
                meshAdapter, "NonlinearThermalOperator", input_db, thermalModel ) );

    AMP::shared_ptr<AMP::Operator::ElementPhysicsModel> oxygenModel;
    AMP::shared_ptr<AMP::Operator::NonlinearBVPOperator> nonlinearOxygenOperator =
        AMP::dynamic_pointer_cast<AMP::Operator::NonlinearBVPOperator>(
            AMP::Operator::OperatorBuilder::createOperator(
                meshAdapter, "NonlinearOxygenOperator", input_db, oxygenModel ) );

    //----------------------------------------------------------------------------------------------------------------------------------------------//
    // create a column rhs operator object with the nonlinear thermal in it for use in the nonlinear
    // problem definition
    AMP::shared_ptr<AMP::Operator::OperatorParameters> params;
    AMP::shared_ptr<AMP::Operator::ColumnOperator> columnNonlinearRhsOperator(
        new AMP::Operator::ColumnOperator( params ) );
    columnNonlinearRhsOperator->append( nonlinearThermalOperator );
    columnNonlinearRhsOperator->append( nonlinearOxygenOperator );
    // ---------------------------------------------------------------------------------------
    // create linear BVP operators
    AMP::shared_ptr<AMP::Operator::LinearBVPOperator> linearThermalOperator;
    AMP::shared_ptr<AMP::Operator::LinearBVPOperator> linearOxygenOperator;
    linearThermalOperator = AMP::dynamic_pointer_cast<AMP::Operator::LinearBVPOperator>(
        AMP::Operator::OperatorBuilder::createOperator(
            meshAdapter, "LinearThermalOperator", input_db, thermalModel ) );
    linearOxygenOperator = AMP::dynamic_pointer_cast<AMP::Operator::LinearBVPOperator>(
        AMP::Operator::OperatorBuilder::createOperator(
            meshAdapter, "LinearOxygenOperator", input_db, oxygenModel ) );
    //----------------------------------------------------------------------------------------------------------------------------------------------//
    // create a column rhs operator object with the linear thermal in it for use in the linear
    // problem definition
    AMP::shared_ptr<AMP::Operator::OperatorParameters> rhsparams;
    AMP::shared_ptr<AMP::Operator::ColumnOperator> columnLinearRhsOperator(
        new AMP::Operator::ColumnOperator( rhsparams ) );
    columnLinearRhsOperator->append( linearThermalOperator );
    columnLinearRhsOperator->append( linearOxygenOperator );

    // ---------------------------------------------------------------------------------------
    // create a mass linear BVP operator
    AMP::shared_ptr<AMP::Operator::ElementPhysicsModel> massThermalModel;
    AMP::shared_ptr<AMP::Operator::ElementPhysicsModel> massOxygenModel;
    AMP::shared_ptr<AMP::Operator::LinearBVPOperator> massThermalOp;
    AMP::shared_ptr<AMP::Operator::LinearBVPOperator> massOxygenOp;
    massThermalOp = AMP::dynamic_pointer_cast<AMP::Operator::LinearBVPOperator>(
        AMP::Operator::OperatorBuilder::createOperator(
            meshAdapter, "MassThermalOperator", input_db, massThermalModel ) );
    massOxygenOp = AMP::dynamic_pointer_cast<AMP::Operator::LinearBVPOperator>(
        AMP::Operator::OperatorBuilder::createOperator(
            meshAdapter, "MassOxygenOperator", input_db, massOxygenModel ) );

    // ---------------------------------------------------------------------------------------
    // create a column mass operator object for use in the nonlinear and linear problem definition
    AMP::shared_ptr<AMP::Operator::OperatorParameters> massParams;
    AMP::shared_ptr<AMP::Operator::ColumnOperator> columnMassOperator(
        new AMP::Operator::ColumnOperator( massParams ) );
    columnMassOperator->append( massThermalOp );
    columnMassOperator->append( massOxygenOp );

    // ---------------------------------------------------------------------------------------
    // create a  time operator for use in the preconditioner
    AMP::shared_ptr<AMP::InputDatabase> timeOperator_db(
        new AMP::InputDatabase( "TimeOperatorDatabase" ) );
    timeOperator_db->putDouble( "CurrentDt", 0.01 );
    timeOperator_db->putString( "name", "TimeOperator" );
    timeOperator_db->putBool( "bLinearMassOperator", true );
    timeOperator_db->putBool( "bLinearRhsOperator", false );
    timeOperator_db->putDouble( "ScalingFactor", 1.0 / 0.01 );

    AMP::shared_ptr<AMP::TimeIntegrator::TimeOperatorParameters> timeOperatorParameters(
        new AMP::TimeIntegrator::TimeOperatorParameters( timeOperator_db ) );
    timeOperatorParameters->d_pRhsOperator  = columnLinearRhsOperator;
    timeOperatorParameters->d_pMassOperator = columnMassOperator;
    timeOperatorParameters->d_Mesh          = meshAdapter;
    AMP::shared_ptr<AMP::TimeIntegrator::ColumnTimeOperator> columnLinearTimeOperator(
        new AMP::TimeIntegrator::ColumnTimeOperator( timeOperatorParameters ) );

    // ---------------------------------------------------------------------------------------
    // create vectors for initial conditions (IC) and time derivative at IC
    AMP::shared_ptr<AMP::Operator::DiffusionNonlinearFEOperator> thermalVolumeOperator =
        AMP::dynamic_pointer_cast<AMP::Operator::DiffusionNonlinearFEOperator>(
            nonlinearThermalOperator->getVolumeOperator() );
    AMP::shared_ptr<AMP::Operator::FickSoretNonlinearFEOperator> oxygenVolumeOperator =
        AMP::dynamic_pointer_cast<AMP::Operator::FickSoretNonlinearFEOperator>(
            nonlinearOxygenOperator->getVolumeOperator() );

    // note that the input variable for the time integrator and time operator will be a
    // multivariable
    AMP::shared_ptr<AMP::LinearAlgebra::MultiVariable> inputVar(
        new AMP::LinearAlgebra::MultiVariable( "inputVariable" ) );
    inputVar->add( thermalVolumeOperator->getOutputVariable() );
    inputVar->add( oxygenVolumeOperator->getOutputVariable() );

    AMP::shared_ptr<AMP::LinearAlgebra::MultiVariable> outputVar =
        AMP::dynamic_pointer_cast<AMP::LinearAlgebra::MultiVariable>(
            columnNonlinearRhsOperator->getOutputVariable() );

    AMP::LinearAlgebra::Vector::shared_ptr initialCondition =
        AMP::LinearAlgebra::createVector( nodalDofMap, outputVar );
    AMP::LinearAlgebra::Vector::shared_ptr initialConditionPrime =
        AMP::LinearAlgebra::createVector( nodalDofMap, outputVar );
    AMP::LinearAlgebra::Vector::shared_ptr f =
        AMP::LinearAlgebra::createVector( nodalDofMap, outputVar );

    f->zero();

    // ---------------------------------------------------------------------------------------
    //  create neutronics source
    AMP_INSIST( input_db->keyExists( "NeutronicsOperator" ),
                "Key ''NeutronicsOperator'' is missing!" );
    AMP::shared_ptr<AMP::Database> neutronicsOp_db = input_db->getDatabase( "NeutronicsOperator" );
    AMP::shared_ptr<AMP::Operator::NeutronicsRhsParameters> neutronicsParams(
        new AMP::Operator::NeutronicsRhsParameters( neutronicsOp_db ) );
    neutronicsParams->d_Mesh = meshAdapter;
    AMP::shared_ptr<AMP::Operator::NeutronicsRhs> neutronicsOperator(
        new AMP::Operator::NeutronicsRhs( neutronicsParams ) );

    AMP::LinearAlgebra::Variable::shared_ptr SpecificPowerVar =
        neutronicsOperator->getOutputVariable();
    AMP::LinearAlgebra::Vector::shared_ptr SpecificPowerVec =
        AMP::LinearAlgebra::createVector( gaussPointDofMap, SpecificPowerVar );

    // create the following shared pointers for ease of use
    AMP::LinearAlgebra::Vector::shared_ptr nullVec;

    neutronicsOperator->apply( nullVec, SpecificPowerVec );

    //  Integrate Nuclear Rhs over Density * GeomType::Volume //

    AMP_INSIST( input_db->keyExists( "VolumeIntegralOperator" ), "key missing!" );

    AMP::shared_ptr<AMP::Operator::ElementPhysicsModel> sourceTransportModel;
    AMP::shared_ptr<AMP::Operator::VolumeIntegralOperator> sourceOperator =
        AMP::dynamic_pointer_cast<AMP::Operator::VolumeIntegralOperator>(
            AMP::Operator::OperatorBuilder::createOperator(
                meshAdapter, "VolumeIntegralOperator", input_db, sourceTransportModel ) );

    // Create the power (heat source) vector.
    AMP::LinearAlgebra::Variable::shared_ptr powerInWattsVar = sourceOperator->getOutputVariable();
    AMP::LinearAlgebra::Vector::shared_ptr powerInWattsVec =
        AMP::LinearAlgebra::createVector( nodalDofMap, powerInWattsVar );
    powerInWattsVec->zero();

    // convert the vector of specific power to power for a given basis.
    sourceOperator->apply( SpecificPowerVec, powerInWattsVec );
    //----------------------------------------------------------------------------------------------------------------------------------------------//
    // set initial conditions, initialize created vectors

    int zeroGhostWidth           = 0;
    AMP::Mesh::MeshIterator node = meshAdapter->getIterator( AMP::Mesh::GeomType::Vertex, zeroGhostWidth );
    AMP::Mesh::MeshIterator end_node = node.end();

    AMP::LinearAlgebra::VS_Mesh vectorSelector1( meshAdapter );
    AMP::LinearAlgebra::VS_Mesh vectorSelector2( meshAdapter );

    AMP::shared_ptr<AMP::LinearAlgebra::Vector> thermalIC =
        initialCondition->select( vectorSelector1, ( outputVar->getVariable( 0 ) )->getName() );
    AMP::shared_ptr<AMP::LinearAlgebra::Vector> oxygenIC =
        initialCondition->select( vectorSelector2, ( outputVar->getVariable( 1 ) )->getName() );
    int counter = 0;
    for ( ; node != end_node; ++node ) {
        counter++;

        std::vector<size_t> gid;
        nodalDofMap->getDOFs( node->globalID(), gid );

        double px = ( node->coord() )[0];
        double py = ( node->coord() )[1];
        double pz = ( node->coord() )[2];

        double tval = __INIT_THERMAL_FN__( px, py, pz, 0 );
        double oval = __INIT_OXYGEN_FN__( px, py, pz, 0 );
        std::cout << "tval = " << tval << std::endl;
        std::cout << "oval = " << oval << std::endl;

        std::cout << "counter = " << counter << "bndGlobalIds.size() = " << gid.size() << std::endl;
        for ( auto &elem : gid ) {
            thermalIC->setValueByGlobalID( elem, tval );
            oxygenIC->setValueByGlobalID( elem, oval );
        } // end for i
    }     // end for node

    // ** please do not set the time derivative to be non-zero!!
    // ** as this causes trouble with the boundary - BP, 07/16/2010
    initialConditionPrime->zero();

    AMP::LinearAlgebra::Vector::shared_ptr thermalRhs =
        f->select( vectorSelector1, outputVar->getVariable( 0 )->getName() );
    // create a copy of the rhs which can be modified at each time step (maybe)
    thermalRhs->copyVector( powerInWattsVec );
    // modify the rhs to take into account boundary conditions
    nonlinearThermalOperator->modifyRHSvector( f );
    nonlinearThermalOperator->modifyInitialSolutionVector( initialCondition );
    nonlinearOxygenOperator->modifyRHSvector( f );
    nonlinearOxygenOperator->modifyInitialSolutionVector( initialCondition );

    // ---------------------------------------------------------------------------------------
    // create a preconditioner

    // get the ida database
    AMP_INSIST( input_db->keyExists( "IDATimeIntegrator" ),
                "Key ''IDATimeIntegrator'' is missing!" );
    AMP::shared_ptr<AMP::Database> ida_db = input_db->getDatabase( "IDATimeIntegrator" );
    // initialize the column preconditioner which is a diagonal block preconditioner
    AMP::shared_ptr<AMP::Database> columnPreconditioner_db =
        ida_db->getDatabase( "Preconditioner" );
    AMP::shared_ptr<AMP::Solver::SolverStrategyParameters> columnPreconditionerParams(
        new AMP::Solver::SolverStrategyParameters( columnPreconditioner_db ) );
    if ( columnPreconditionerParams.get() == nullptr ) {
        ut->failure( "Testing SolverStrategyParameters's constructor: FAIL" );
    } else {
        ut->passes( "Testing SolverStrategyParameters's constructor: PASS" );
    }

    columnPreconditionerParams->d_pOperator = columnLinearTimeOperator;
    AMP::shared_ptr<AMP::Solver::ColumnSolver> columnPreconditioner(
        new AMP::Solver::ColumnSolver( columnPreconditionerParams ) );

    AMP::shared_ptr<AMP::Database> thermalPreconditioner_db =
        columnPreconditioner_db->getDatabase( "thermalPreconditioner" );
    AMP::shared_ptr<AMP::Database> oxygenPreconditioner_db =
        columnPreconditioner_db->getDatabase( "oxygenPreconditioner" );
    AMP::shared_ptr<AMP::Solver::SolverStrategyParameters> thermalPreconditionerParams(
        new AMP::Solver::SolverStrategyParameters( thermalPreconditioner_db ) );
    AMP::shared_ptr<AMP::Solver::SolverStrategyParameters> oxygenPreconditionerParams(
        new AMP::Solver::SolverStrategyParameters( oxygenPreconditioner_db ) );
    thermalPreconditionerParams->d_pOperator = columnLinearTimeOperator->getOperator( 0 );
    oxygenPreconditionerParams->d_pOperator  = columnLinearTimeOperator->getOperator( 1 );
    AMP::shared_ptr<AMP::Solver::TrilinosMLSolver> linearThermalPreconditioner(
        new AMP::Solver::TrilinosMLSolver( thermalPreconditionerParams ) );
    AMP::shared_ptr<AMP::Solver::TrilinosMLSolver> linearOxygenPreconditioner(
        new AMP::Solver::TrilinosMLSolver( oxygenPreconditionerParams ) );

    columnPreconditioner->append( linearThermalPreconditioner );
    columnPreconditioner->append( linearOxygenPreconditioner );

    if ( columnPreconditioner.get() == nullptr ) {
        ut->failure( "Testing column preconditioner's constructor: FAIL" );
    } else {
        ut->passes( "Testing column preconditioner's constructor: PASS" );
    }

    // ---------------------------------------------------------------------------------------
    // create the IDA time integrator
    AMP::shared_ptr<AMP::TimeIntegrator::IDATimeIntegratorParameters> time_Params(
        new AMP::TimeIntegrator::IDATimeIntegratorParameters( ida_db ) );

    if ( ( time_Params.get() ) == nullptr ) {
        ut->failure( "Testing IDATimeIntegratorParameters' Constructor" );
    } else {
        ut->passes( "Testing IDATimeIntegratorParameters' Constructor" );
    }

    time_Params->d_pMassOperator   = columnMassOperator;
    time_Params->d_operator        = columnNonlinearRhsOperator;
    time_Params->d_pPreconditioner = columnPreconditioner;

    time_Params->d_ic_vector       = initialCondition;
    time_Params->d_ic_vector_prime = initialConditionPrime;

    time_Params->d_pSourceTerm = f;
    time_Params->d_object_name = "IDATimeIntegratorParameters";

    std::cout << "Before IDATimeIntegrator" << std::endl;
    AMP::shared_ptr<AMP::TimeIntegrator::IDATimeIntegrator> pIDATimeIntegrator(
        new AMP::TimeIntegrator::IDATimeIntegrator( time_Params ) );

    if ( pIDATimeIntegrator.get() == nullptr ) {
        ut->failure( "Testing IDATimeIntegrator's constructor" );
    } else {
        ut->passes( "Tested IDATimeIntegrator's constructor" );
    }

    // ---------------------------------------------------------------------------------------
    // step in time
    int retval          = 0;
    double current_time = 0;
    double maxT         = 0;
    double maxO         = 0;
    // double abs_error=0.0;
    double minT = 0;
    double minO = 0;
    // double rel_error=0.0;
    // double exact_sol=0.0;
    int j = 1;
    while ( pIDATimeIntegrator->getCurrentTime() < pIDATimeIntegrator->getFinalTime() ) {
        retval = pIDATimeIntegrator->advanceSolution( pIDATimeIntegrator->getCurrentDt(), 0 );
        // pIDATimeIntegrator->updateSolution();
        current_time = pIDATimeIntegrator->getCurrentTime();

        std::cout << j++ << "-th timestep" << std::endl;
        if ( retval == 0 ) {
            ut->passes( "Testing IDATimeIntegrator's advanceSolution. PASS!!" );
        } else {
            ut->failure( "Tested IDATimeIntegrator's advanceSolution. FAIL!!" );
        }

        AMP::shared_ptr<AMP::LinearAlgebra::Vector> currentSolution =
            pIDATimeIntegrator->getCurrentSolution();
        AMP::shared_ptr<AMP::LinearAlgebra::Vector> currentThermalSolution =
            currentSolution->subsetVectorForVariable( thermalVolumeOperator->getOutputVariable() );
        AMP::shared_ptr<AMP::LinearAlgebra::Vector> currentOxygenSolution =
            currentSolution->subsetVectorForVariable( oxygenVolumeOperator->getOutputVariable() );

        maxT = currentThermalSolution->max();
        minT = currentThermalSolution->min();

        maxO = currentOxygenSolution->max();
        minO = currentOxygenSolution->min();

        std::cout << "current_time = " << current_time << std::endl;
        std::cout << "max val of the current thermal solution = " << maxT << std::endl;
        std::cout << "min val of the current thermal solution = " << minT << std::endl;
        std::cout << "max val of the current oxygen solution = " << maxO << std::endl;
        std::cout << "min val of the current oxygen solution = " << minO << std::endl;
    }


    AMP::AMPManager::shutdown();

    if ( ut->NumFailLocal() == 0 ) {
        ut->passes( "testIDATimeIntegrator successful" );
    }
}


//---------------------------------------------------------------------------//

int main( int argc, char *argv[] )
{
    AMP::AMPManager::startup( argc, argv );
    AMP::UnitTest ut;

    IDATimeIntegratorTest( &ut );

    ut.report();

    int num_failed = ut.NumFailGlobal();
    AMP::AMPManager::shutdown();
    return num_failed;
}

//---------------------------------------------------------------------------//
//                        end of SundialsVectorTest.cc
//---------------------------------------------------------------------------//
