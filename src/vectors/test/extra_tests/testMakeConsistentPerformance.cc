#include "utils/Database.h"
#include "ampmesh/Mesh.h"
#include "discretization/DOF_Manager.h"
#include "discretization/simpleDOF_Manager.h"
#include "utils/AMPManager.h"
#include "utils/AMP_MPI.h"
#include "utils/InputDatabase.h"
#include "utils/InputManager.h"
#include "utils/PIO.h"
#include "utils/UnitTest.h"
#include "utils/Utilities.h"
#include "vectors/Variable.h"
#include "vectors/VectorBuilder.h"


void runTest( AMP::UnitTest *ut )
{

    // Get the Mesh database and create the mesh parameters
    AMP::AMP_MPI globalComm( AMP_COMM_WORLD );
    AMP::shared_ptr<AMP::MemoryDatabase> mesh_db( new AMP::MemoryDatabase( "Mesh" ) );
    mesh_db->putInteger( "dim", 3 );
    mesh_db->putString( "MeshName", "mesh1" );
    mesh_db->putString( "MeshType", "libMesh" );
    mesh_db->putString( "FileName", "pellet_1x.e" );
    AMP::shared_ptr<AMP::Mesh::MeshParameters> params( new AMP::Mesh::MeshParameters( mesh_db ) );
    params->setComm( globalComm );

    // Create the meshes from the input database
    AMP::shared_ptr<AMP::Mesh::Mesh> mesh = AMP::Mesh::Mesh::buildMesh( params );

    // Create a simple DOFManager
    int DOFsPerNode     = 1;
    std::string varName = "test";
    AMP::LinearAlgebra::Variable::shared_ptr nodalVariable(
        new AMP::LinearAlgebra::Variable( varName ) );
    AMP::Discretization::DOFManagerParameters::shared_ptr DOFparams(
        new AMP::Discretization::DOFManagerParameters( mesh ) );
    AMP::Discretization::DOFManager::shared_ptr DOFs =
        AMP::Discretization::simpleDOFManager::create( mesh, AMP::Mesh::GeomType::Vertex, 1, DOFsPerNode );

    // Create the vectors
    AMP::LinearAlgebra::Vector::shared_ptr dummy;
    AMP::LinearAlgebra::Vector::shared_ptr v1 = createVector( DOFs, nodalVariable );
    AMP::LinearAlgebra::Vector::shared_ptr v2 = createVector( DOFs, nodalVariable );

    // Initialize the vectors
    v1->setToScalar( 0.0 );
    v2->setToScalar( 0.0 );

    // Time makeConsistentSet
    globalComm.barrier();
    double start_time = AMP::AMP_MPI::time();
    v1->makeConsistent( AMP::LinearAlgebra::Vector::ScatterType::CONSISTENT_SET );
    globalComm.barrier();
    double end_time = AMP::AMP_MPI::time();
    std::cout << std::endl << "Time for makeConsistent: " << end_time - start_time << std::endl;

    // Print the number of ghost values in the communication list
    AMP::LinearAlgebra::CommunicationList::shared_ptr communicationList =
        v1->getCommunicationList();
    std::vector<size_t> ghost_ids = communicationList->getGhostIDList();
    size_t N_ghosts               = globalComm.sumReduce( ghost_ids.size() );
    size_t N_ghosts2 = globalComm.sumReduce( mesh->numGhostElements( AMP::Mesh::GeomType::Vertex, 1 ) );
    std::cout << std::endl << "There are " << N_ghosts << " global ghost values" << std::endl;
    std::cout << std::endl
              << "There are " << N_ghosts2 << " global ghost values in the iterator" << std::endl;

    ut->passes( "Test ran to completion" );
}


int main( int argc, char *argv[] )
{
    AMP::AMPManager::startup( argc, argv );
    AMP::UnitTest ut;

    runTest( &ut );
    ut.report();

    int num_failed = ut.NumFailGlobal();
    AMP::AMPManager::shutdown();
    return num_failed;
}
