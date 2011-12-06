#include "utils/Database.h"
#include "utils/InputDatabase.h"
#include "utils/InputManager.h"
#include "utils/AMP_MPI.h"
#include "utils/AMPManager.h"
#include "utils/UnitTest.h"
#include "utils/Utilities.h"
#include "utils/PIO.h"
#include "ampmesh/MeshManager.h"



void  runTest ( AMP::UnitTest *ut )
{

    // Get the Mesh database and create the mesh parameters
    AMP::AMP_MPI globalComm(AMP_COMM_WORLD);
    boost::shared_ptr<AMP::MemoryDatabase> input_db(new AMP::MemoryDatabase("input"));
    input_db->putInteger("NumberOfMeshes",1);
    input_db->putInteger("NumberOfMeshToMeshMaps",0);
    input_db->putInteger("DomainDecomposition",0);
    boost::shared_ptr<AMP::Database> mesh_db = input_db->putDatabase("Mesh_1");
    mesh_db->putString("Filename","pellet_1x.e");
    mesh_db->putString("MeshName","pellet_1");
    mesh_db->putDouble("x_offset",0.0);
    mesh_db->putDouble("y_offset",0.0);
    mesh_db->putDouble("z_offset",0.0);
    mesh_db->putInteger("NumberOfElements",1326);
    //mesh_db->putString("DatabaseName","db1");
    AMP::Mesh::MeshManagerParameters::shared_ptr  meshmgrParams ( new AMP::Mesh::MeshManagerParameters ( input_db ) );
    AMP::Mesh::MeshManager::shared_ptr  manager ( new AMP::Mesh::MeshManager ( meshmgrParams ) );
    AMP::Mesh::MeshManager::Adapter::shared_ptr mesh = manager->getMesh( "pellet_1" );
    globalComm.barrier();

    // Create a nodal variable
    std::string varName = "test";
    AMP::LinearAlgebra::Variable::shared_ptr variable = AMP::LinearAlgebra::Variable::shared_ptr( new AMP::Mesh::NodalScalarVariable( varName, mesh ) );

    // Create the vectors
    AMP::LinearAlgebra::Vector::shared_ptr  dummy;
    AMP::LinearAlgebra::Vector::shared_ptr  v1 = manager->createVector( variable );
    AMP::LinearAlgebra::Vector::shared_ptr  v2 = manager->createVector( variable );

    // Initialize the vectors
    v1->setToScalar(0.0);
    v2->setToScalar(0.0);

    // Time makeConsistentSet
    globalComm.barrier();
    double start_time = AMP::AMP_MPI::time();
    v1->makeConsistent ( AMP::LinearAlgebra::Vector::CONSISTENT_SET );
    globalComm.barrier();
    double end_time = AMP::AMP_MPI::time();
    std::cout << std::endl << "Time for makeConsistent: " << end_time-start_time << std::endl;

    // Print the number of ghost values in the communication list
    AMP::LinearAlgebra::CommunicationList::shared_ptr  communicationList = v1->getCommunicationList();
    std::vector<unsigned int> ghost_ids = communicationList->getGhostIDList();
    size_t N_ghosts = globalComm.sumReduce(ghost_ids.size());
    std::cout << std::endl << "There are " << N_ghosts << " global ghost values" << std::endl;

}


int main(int argc, char *argv[])
{
    AMP::AMPManager::startup(argc, argv);
    AMP::UnitTest ut;

    runTest ( &ut );
    ut.report();

    int num_failed = ut.NumFailGlobal();
    AMP::AMPManager::shutdown();
    return num_failed;
}   

