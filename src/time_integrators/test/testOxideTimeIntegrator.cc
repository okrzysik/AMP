#include "utils/AMPManager.h"
#include "utils/AMP_MPI.h"
#include "utils/InputManager.h"
#include "utils/UnitTest.h"

#include "ampmesh/Mesh.h"
#include "discretization/simpleDOF_Manager.h"
#include "utils/Writer.h"
#include "vectors/Vector.h"
#include "vectors/VectorBuilder.h"

#include "time_integrators/oxide/OxideTimeIntegrator.h"

void OxideTest( AMP::UnitTest *ut, std::string input_file )
{

    // Load the input file
    std::string log_file = input_file + ".log";
    AMP::PIO::logOnlyNodeZero( log_file );
    AMP::shared_ptr<AMP::InputDatabase> input_db( new AMP::InputDatabase( "input_db" ) );
    AMP::InputManager::getManager()->parseInputFile( input_file, input_db );
    input_db->printClassData( AMP::plog );

    // Get the Mesh database and create the mesh parameters
    AMP::shared_ptr<AMP::Database> database = input_db->getDatabase( "Mesh" );
    AMP::shared_ptr<AMP::Mesh::MeshParameters> params( new AMP::Mesh::MeshParameters( database ) );
    AMP::AMP_MPI globalComm( AMP_COMM_WORLD );
    params->setComm( globalComm );

    // Create the meshes from the input database
    AMP::Mesh::Mesh::shared_ptr manager = AMP::Mesh::Mesh::buildMesh( params );
    AMP::Mesh::Mesh::shared_ptr mesh    = manager->Subset( "clad" );
    globalComm.barrier();

    // Create the surface mesh that we will use to create the oxide layer
    AMP::Mesh::Mesh::shared_ptr surface =
        mesh->Subset( mesh->getBoundaryIDIterator( AMP::Mesh::Face, 4, 0 ) );
    surface->setName( "clad_surface" );

    // Create the temperature profile
    AMP::Discretization::DOFManager::shared_ptr DOF =
        AMP::Discretization::simpleDOFManager::create( mesh, AMP::Mesh::Vertex, 1, 1, true );
    AMP::LinearAlgebra::Variable::shared_ptr temp_var(
        new AMP::LinearAlgebra::Variable( "temperature" ) );
    AMP::LinearAlgebra::Vector::shared_ptr temp_vec =
        AMP::LinearAlgebra::createVector( DOF, temp_var, true );
    AMP::Mesh::MeshIterator iterator = mesh->getIterator( AMP::Mesh::Vertex );
    double T0                        = input_db->getDoubleWithDefault( "T0", 650 );
    std::vector<size_t> dofs;
    for ( size_t i = 0; i < iterator.size(); i++ ) {
        std::vector<double> coord = iterator->coord();
        DOF->getDOFs( iterator->globalID(), dofs );
        temp_vec->setValueByGlobalID( dofs[0], T0 + 100 * coord[2] );
        ++iterator;
    }
    temp_vec->makeConsistent( AMP::LinearAlgebra::Vector::CONSISTENT_SET );
    AMP_ASSERT( fabs( temp_vec->min() - T0 ) / T0 < 1e-9 );

    // Create the oxide time integrator
    AMP::shared_ptr<AMP::TimeIntegrator::OxideTimeIntegratorParameters> parameters(
        new AMP::TimeIntegrator::OxideTimeIntegratorParameters(
            AMP::shared_ptr<AMP::Database>() ) );
    parameters->d_mesh = surface;
    parameters->d_temp = temp_vec;
    parameters->depth  = 1e-3;
    AMP::TimeIntegrator::TimeIntegrator::shared_ptr timeIntegrator(
        new AMP::TimeIntegrator::OxideTimeIntegrator( parameters ) );
    AMP::LinearAlgebra::Vector::shared_ptr solution = timeIntegrator->getCurrentSolution();
    AMP::LinearAlgebra::Variable::shared_ptr oxide_var(
        new AMP::LinearAlgebra::Variable( "oxide" ) );
    AMP::LinearAlgebra::Variable::shared_ptr alpha_var(
        new AMP::LinearAlgebra::Variable( "alpha" ) );
    AMP::LinearAlgebra::Vector::shared_ptr oxide = solution->subsetVectorForVariable( oxide_var );
    AMP::LinearAlgebra::Vector::shared_ptr alpha = solution->subsetVectorForVariable( alpha_var );

// Register the data with the silo writer
#ifdef USE_EXT_SILO
    AMP::Utilities::Writer::shared_ptr siloWriter = AMP::Utilities::Writer::buildWriter( "Silo" );
    siloWriter->registerVector( temp_vec, mesh, AMP::Mesh::Vertex, "temperature" );
    siloWriter->registerVector( oxide, surface, AMP::Mesh::Vertex, "oxide_thickness" );
    siloWriter->registerVector( alpha, surface, AMP::Mesh::Vertex, "alpha_thickness" );
#endif

    // Run the time integration
    double integration_time   = 0.0;
    std::vector<double> times = input_db->getDoubleArray( "Time" );
    for ( size_t i = 0; i < times.size(); i++ ) {
        // Advance the solution
        double dT = times[i] - timeIntegrator->getCurrentTime();
        globalComm.barrier();
        double t0 = AMP::AMP_MPI::time();
        timeIntegrator->advanceSolution( dT, false );
        globalComm.barrier();
        integration_time += AMP::AMP_MPI::time() - t0;
#ifdef USE_EXT_SILO
        siloWriter->writeFile( input_file, i );
#endif
        // Check the solution
        if ( input_db->keyExists( "oxide" ) && input_db->keyExists( "alpha" ) ) {
            if ( times[i] == 0.0 )
                continue;
            std::vector<double> oxide_solution = input_db->getDoubleArray( "oxide" );
            std::vector<double> alpha_solution = input_db->getDoubleArray( "alpha" );
            double err_oxide                   = oxide->min() - oxide_solution[i];
            double err_alpha                   = alpha->min() - alpha_solution[i];
            if ( fabs( err_oxide / oxide_solution[i] ) < 1e-3 )
                ut->passes( "oxide solution matches" );
            else
                ut->failure( "oxide solution matches" );
            if ( fabs( err_alpha / alpha_solution[i] ) < 1e-3 )
                ut->passes( "alpha solution matches" );
            else
                ut->failure( "alpha solution matches" );
        }
    }
    AMP::pout << "Time required for integration: " << integration_time << std::endl;

    ut->passes( "Test runs to completion" );
}


int main( int argc, char *argv[] )
{

    AMP::AMPManager::startup( argc, argv );
    AMP::UnitTest ut;

    OxideTest( &ut, "input_testOxideTimeIntegrator-1" );

    ut.report();

    int num_failed = ut.NumFailGlobal();
    AMP::AMPManager::shutdown();
    return num_failed;
}
