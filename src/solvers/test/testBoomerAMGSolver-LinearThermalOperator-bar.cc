#include "AMP/ampmesh/Mesh.h"
#include "AMP/discretization/DOF_Manager.h"
#include "AMP/discretization/simpleDOF_Manager.h"
#include "AMP/materials/Material.h"
#include "AMP/operators/ElementOperationFactory.h"
#include "AMP/operators/ElementPhysicsModelFactory.h"
#include "AMP/operators/LinearBVPOperator.h"
#include "AMP/operators/NeutronicsRhs.h"
#include "AMP/operators/OperatorBuilder.h"
#include "AMP/operators/boundary/DirichletMatrixCorrection.h"
#include "AMP/operators/diffusion/DiffusionLinearElement.h"
#include "AMP/operators/diffusion/DiffusionLinearFEOperator.h"
#include "AMP/operators/diffusion/DiffusionTransportModel.h"
#include "AMP/operators/libmesh/VolumeIntegralOperator.h"
#include "AMP/solvers/hypre/BoomerAMGSolver.h"
#include "AMP/utils/AMPManager.h"
#include "AMP/utils/Database.h"
#include "AMP/utils/InputDatabase.h"
#include "AMP/utils/InputManager.h"
#include "AMP/utils/PIO.h"
#include "AMP/utils/UnitTest.h"
#include "AMP/utils/Utilities.h"
#include "AMP/utils/Writer.h"
#include "AMP/utils/shared_ptr.h"
#include "AMP/vectors/Variable.h"
#include "AMP/vectors/Vector.h"
#include "AMP/vectors/VectorBuilder.h"

#include <fstream>
#include <limits>
#include <string>

#define ITFAILS ut.failure( __LINE__ );
#define UNIT_TEST( a ) \
    if ( !( a ) )      \
        ut.failure( __LINE__ );

void linearThermalTest( AMP::UnitTest *ut )
{
    // Input and output file names
    //  #include <string>
    std::string exeName( "testBoomerAMGSolver-LinearThermalOperator-bar" );
    std::string input_file = "input_" + exeName;
    std::string log_file   = "output_" + exeName;
    ////////////////////////////////////
    //    INITIALIZE THE PROBLEM      //
    ////////////////////////////////////

    // Construct a smart pointer to a new database.
    //  #include "utils/shared_ptr.h"
    //  #include "utils/InputDatabase.h"
    AMP::shared_ptr<AMP::InputDatabase> input_db( new AMP::InputDatabase( "input_db" ) );

    // Fill the database from the input file.
    //  #include "utils/InputManager.h"
    AMP::InputManager::getManager()->parseInputFile( input_file, input_db );
    input_db->printClassData( AMP::plog );


    // Print from all cores into the output files
    //   #include "utils/PIO.h"
    AMP::PIO::logAllNodes( log_file );

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
    ////////////////////////////////////
    //  CREATE THE NEUTRONICS SOURCE  //
    ////////////////////////////////////
    AMP_INSIST( input_db->keyExists( "NeutronicsOperator" ),
                "Key ''NeutronicsOperator'' is missing!" );
    AMP::shared_ptr<AMP::Operator::ElementPhysicsModel> unusedModel;
    AMP::shared_ptr<AMP::Operator::NeutronicsRhs> neutronicsOperator =
        AMP::dynamic_pointer_cast<AMP::Operator::NeutronicsRhs>(
            AMP::Operator::OperatorBuilder::createOperator(
                meshAdapter, "NeutronicsOperator", input_db, unusedModel ) );

    neutronicsOperator->setTimeStep( 0. );

    AMP::LinearAlgebra::Variable::shared_ptr SpecificPowerVar =
        neutronicsOperator->getOutputVariable();
    AMP::LinearAlgebra::Vector::shared_ptr SpecificPowerVec =
        AMP::LinearAlgebra::createVector( gaussPointDofMap, SpecificPowerVar );

    neutronicsOperator->apply( nullVec, SpecificPowerVec );

    /////////////////////////////////////////////////////
    //  Integrate Nuclear Source over Desnity * Volume //
    /////////////////////////////////////////////////////

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

    // convert the vector of specific power to power for a given basis.
    sourceOperator->apply( SpecificPowerVec, PowerInWattsVec );

    ////////////////////////////////////
    //   CREATE THE THERMAL OPERATOR  //
    ////////////////////////////////////
    AMP::shared_ptr<AMP::Operator::ElementPhysicsModel> transportModel;
    AMP::shared_ptr<AMP::Operator::LinearBVPOperator> diffusionOperator =
        AMP::dynamic_pointer_cast<AMP::Operator::LinearBVPOperator>(
            AMP::Operator::OperatorBuilder::createOperator(
                meshAdapter, "DiffusionBVPOperator", input_db, transportModel ) );


    AMP::LinearAlgebra::Vector::shared_ptr TemperatureInKelvinVec =
        AMP::LinearAlgebra::createVector( nodalDofMap, diffusionOperator->getInputVariable() );
    AMP::LinearAlgebra::Vector::shared_ptr RightHandSideVec =
        AMP::LinearAlgebra::createVector( nodalDofMap, diffusionOperator->getOutputVariable() );
    AMP::LinearAlgebra::Vector::shared_ptr ResidualVec =
        AMP::LinearAlgebra::createVector( nodalDofMap, diffusionOperator->getOutputVariable() );

    RightHandSideVec->copyVector( PowerInWattsVec );

    AMP::shared_ptr<AMP::Operator::BoundaryOperator> boundaryOp;
    boundaryOp = diffusionOperator->getBoundaryOperator();

    boundaryOp->addRHScorrection( RightHandSideVec );
    boundaryOp->setRHScorrection( RightHandSideVec );

    // make sure the database on theinput file exists for the linear solver
    AMP_INSIST( input_db->keyExists( "LinearSolver" ), "Key ''LinearSolver'' is missing!" );

    // Read the input file onto a database.
    AMP::shared_ptr<AMP::Database> mlSolver_db = input_db->getDatabase( "LinearSolver" );

    // Fill in the parameters fo the class with the info on the database.
    AMP::shared_ptr<AMP::Solver::SolverStrategyParameters> mlSolverParams(
        new AMP::Solver::SolverStrategyParameters( mlSolver_db ) );

    // Define the operature to be used by the Solver.
    mlSolverParams->d_pOperator = diffusionOperator;

    // Set initial guess
    TemperatureInKelvinVec->setToScalar( 1.0 );

    // Check the initial L2 norm of the solution
    double initSolNorm = TemperatureInKelvinVec->L2Norm();
    std::cout << "Initial Solution Norm: " << initSolNorm << std::endl;

    double rhsNorm = RightHandSideVec->L2Norm();
    std::cout << "RHS Norm: " << rhsNorm << std::endl;

    // Create the ML Solver
    auto mlSolver = std::make_shared<AMP::Solver::BoomerAMGSolver>( mlSolverParams );

    // Use a random initial guess?
    mlSolver->setZeroInitialGuess( false );

    // Solve the prblem.
    mlSolver->solve( RightHandSideVec, TemperatureInKelvinVec );

    // Compute the residual
    diffusionOperator->residual( RightHandSideVec, TemperatureInKelvinVec, ResidualVec );

    // Check the L2 norm of the final residual.
    double finalResidualNorm = ResidualVec->L2Norm();
    std::cout << "Final Residual Norm: " << finalResidualNorm << std::endl;

    if ( finalResidualNorm > 10.0 ) {
        ut->failure( "BoomerAMGSolver successfully solves a linear thermal problem with a nuclear "
                     "source term." );
    } else {
        ut->passes( "BoomerAMGSolver successfully solves a linear thermal problem with a nuclear "
                    "source term." );
    }

    // check the solution
    int zeroGhostWidth = 0;
    AMP::Mesh::MeshIterator iterator =
        meshAdapter->getIterator( AMP::Mesh::GeomType::Vertex, zeroGhostWidth );

    // The analytical solution is:  T = a + b*z + c*z*z
    //   c = -power/2
    //   b = -10*power
    //   a = 300 + 150*power

    double power = 1.;
    double c     = -power / 2.;
    double b     = -10. * power;
    double a     = 300. + 150. * power;
    bool passes  = 1;
    double cal, zee, sol, err;

    // Serial execution
    AMP::AMP_MPI globalComm( AMP_COMM_WORLD );
    for ( int i = 0; i < globalComm.getSize(); i++ ) {
        if ( globalComm.getRank() == i ) {
            std::string filename          = "data_" + exeName;
            int rank                      = globalComm.getRank();
            int nranks                    = globalComm.getSize();
            std::ios_base::openmode omode = std::ios_base::out;
            if ( rank > 0 )
                omode |= std::ios_base::app;
            std::ofstream file( filename.c_str(), omode );
            if ( rank == 0 ) {
                file << "(* x y z analytic calculated relative-error *)" << std::endl;
                file << "formula=" << a << " + " << b << "*z + " << c << "*z^2;" << std::endl;
                file << "results={" << std::endl;
            }
            file.precision( 14 );

            iterator        = iterator.begin();
            size_t numNodes = 0, iNode = 0;
            for ( ; iterator != iterator.end(); iterator++ )
                numNodes++;

            iterator   = iterator.end();
            double mse = 0.0;
            for ( ; iterator != iterator.end(); iterator++ ) {
                std::vector<size_t> gid;
                nodalDofMap->getDOFs( iterator->globalID(), gid );
                cal = TemperatureInKelvinVec->getValueByGlobalID( gid[0] );
                zee = ( iterator->coord() )[2];
                sol = a + b * zee + c * zee * zee;
                err =
                    fabs( cal - sol ) * 2. / ( cal + sol + std::numeric_limits<double>::epsilon() );
                double x, y, z;
                x = ( iterator->coord() )[0];
                y = ( iterator->coord() )[1];
                z = ( iterator->coord() )[2];
                mse += ( sol - cal ) * ( sol - cal );
                file << "{" << x << "," << y << "," << z << "," << sol << "," << cal << "," << err
                     << "}";
                if ( iNode < numNodes - 1 )
                    file << "," << std::endl;
                if ( fabs( cal - sol ) > cal * 1e-3 ) {
                    passes = 0;
                    ut->failure( "Error" );
                }
                iNode++;
            }

            if ( rank == nranks - 1 ) {
                file << "};" << std::endl;
                mse /= ( 1. * iNode );
                mse = sqrt( mse );
                file << "l2err = {" << iNode << "," << mse << "};\n";
            }
            file.close();
        }
        globalComm.barrier();
    }
    if ( passes )
        ut->passes( "The linear thermal solve is verified." );

// Plot the results
#ifdef USE_EXT_SILO
    AMP::Utilities::Writer::shared_ptr siloWriter = AMP::Utilities::Writer::buildWriter( "Silo" );
    siloWriter->registerMesh( meshAdapter );

    siloWriter->registerVector(
        PowerInWattsVec, meshAdapter, AMP::Mesh::GeomType::Vertex, "PowerInWatts" );
    siloWriter->registerVector(
        TemperatureInKelvinVec, meshAdapter, AMP::Mesh::GeomType::Vertex, "TemperatureInKelvin" );
    siloWriter->registerVector( ResidualVec, meshAdapter, AMP::Mesh::GeomType::Vertex, "Residual" );

    siloWriter->writeFile( input_file, 0 );
#endif

    input_db.reset();

    ut->passes( exeName );
}


int main( int argc, char *argv[] )
{
    AMP::AMPManager::startup( argc, argv );
    AMP::UnitTest ut;

    linearThermalTest( &ut );

    ut.report();

    int num_failed = ut.NumFailGlobal();
    AMP::AMPManager::shutdown();
    return num_failed;
}
