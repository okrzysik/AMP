#include "utils/AMPManager.h"
#include "utils/UnitTest.h"
#include "utils/AMP_MPI.h"
#include "utils/InputManager.h"

#include "ampmesh/SiloIO.h"
#include "vectors/Vector.h"
#include "discretization/simpleDOF_Manager.h"
#include "vectors/VectorBuilder.h"

#include "time_integrators/oxide/OxideTimeIntegrator.h"

void OxideTest( AMP::UnitTest *ut, std::string input_file )
{
    
    // Load the input file    
    std::string log_file = input_file + ".log";  
    AMP::PIO::logOnlyNodeZero(log_file);
    boost::shared_ptr<AMP::InputDatabase> input_db(new AMP::InputDatabase("input_db"));
    AMP::InputManager::getManager()->parseInputFile(input_file, input_db);
    input_db->printClassData(AMP::plog);
    
    // Get the Mesh database and create the mesh parameters
    boost::shared_ptr<AMP::Database> database = input_db->getDatabase( "Mesh" );
    boost::shared_ptr<AMP::Mesh::MeshParameters> params(new AMP::Mesh::MeshParameters(database));
    AMP::AMP_MPI globalComm(AMP_COMM_WORLD);
    params->setComm(globalComm);

    // Create the meshes from the input database
    AMP::Mesh::Mesh::shared_ptr manager = AMP::Mesh::Mesh::buildMesh(params);
    AMP::Mesh::Mesh::shared_ptr mesh = manager->Subset("clad");
    globalComm.barrier();

    // Create the surface mesh that we will use to create the oxide layer
    AMP::Mesh::Mesh::shared_ptr surface = mesh->Subset( mesh->getBoundaryIDIterator( AMP::Mesh::Face, 4, 0 ) );
    surface->setName("clad_surface");

    // Create the temperature profile
    AMP::Discretization::DOFManager::shared_ptr DOF = AMP::Discretization::simpleDOFManager::create(mesh,AMP::Mesh::Vertex,1,1,true);
    AMP::LinearAlgebra::Variable::shared_ptr temp_var( new AMP::LinearAlgebra::Variable("temperature") );
    AMP::LinearAlgebra::Vector::shared_ptr temp_vec = AMP::LinearAlgebra::createVector( DOF, temp_var, true );
    temp_vec->setToScalar(650);

    // Create the oxide time integrator
    boost::shared_ptr<AMP::TimeIntegrator::OxideTimeIntegratorParameters> parameters( 
        new AMP::TimeIntegrator::OxideTimeIntegratorParameters(boost::shared_ptr<AMP::Database>()) );
    parameters->d_mesh = surface;
    parameters->d_temp = temp_vec;
    parameters->depth  = 1e-3;
    AMP::TimeIntegrator::TimeIntegrator::shared_ptr timeIntegrator( new AMP::TimeIntegrator::OxideTimeIntegrator( parameters ) );
    AMP::LinearAlgebra::Vector::shared_ptr solution = timeIntegrator->getCurrentSolution();
    AMP::LinearAlgebra::Variable::shared_ptr oxide_var( new AMP::LinearAlgebra::Variable("oxide") );
    AMP::LinearAlgebra::Variable::shared_ptr alpha_var( new AMP::LinearAlgebra::Variable("alpha") );
    AMP::LinearAlgebra::Vector::shared_ptr oxide = solution->subsetVectorForVariable( oxide_var );
    AMP::LinearAlgebra::Vector::shared_ptr alpha = solution->subsetVectorForVariable( alpha_var );
    
    // Register the data with the silo writer
    #ifdef USE_SILO
        AMP::Mesh::SiloIO::shared_ptr  siloWriter( new AMP::Mesh::SiloIO);
        siloWriter->registerVector( temp_vec, mesh, AMP::Mesh::Vertex, "temperature" );
        siloWriter->registerVector( oxide, surface, AMP::Mesh::Vertex, "oxide_thickness" );
        siloWriter->registerVector( alpha, surface, AMP::Mesh::Vertex, "alpha_thickness" );
    #endif
    
    // Run the time integration
    std::vector<double> times = input_db->getDoubleArray("Time");
    for (size_t i=0; i<times.size(); i++) {
        double dT = times[i] - timeIntegrator->getCurrentTime();
        timeIntegrator->advanceSolution( dT, false );
        #ifdef USE_SILO
            siloWriter->writeFile( input_file, i );
        #endif
    }
}


int main(int argc, char *argv[])
{
    
    AMP::AMPManager::startup(argc, argv);
    AMP::UnitTest ut;
   
    OxideTest( &ut, "input_testOxideTimeIntegrator-1" );

    ut.report();

    int num_failed = ut.NumFailGlobal();
    AMP::AMPManager::shutdown();
    return num_failed;
}


