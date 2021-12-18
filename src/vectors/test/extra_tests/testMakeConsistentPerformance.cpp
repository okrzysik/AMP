#include "AMP/IO/PIO.h"
#include "AMP/discretization/DOF_Manager.h"
#include "AMP/discretization/simpleDOF_Manager.h"
#include "AMP/mesh/Mesh.h"
#include "AMP/mesh/MeshParameters.h"
#include "AMP/utils/AMPManager.h"
#include "AMP/utils/AMP_MPI.h"
#include "AMP/utils/Database.h"
#include "AMP/utils/UnitTest.h"
#include "AMP/utils/Utilities.h"
#include "AMP/vectors/Variable.h"
#include "AMP/vectors/VectorBuilder.h"


static void runTest( AMP::UnitTest *ut )
{

    // Get the Mesh database and create the mesh parameters
    AMP::AMP_MPI globalComm( AMP_COMM_WORLD );
    auto mesh_db = std::make_shared<AMP::Database>( "Mesh" );
    // mesh_db->putScalar( "dim", 3 );
    // mesh_db->putScalar( "MeshName", "mesh1" );
    // mesh_db->putScalar( "MeshType", "libMesh" );
    // mesh_db->putScalar( "FileName", "pellet_1x.e" );
    mesh_db->putScalar( "dim", 3 );
    mesh_db->putScalar( "MeshName", "mesh1" );
    mesh_db->putScalar( "MeshType", "AMP" );
    mesh_db->putScalar( "Generator", "cylinder" );
    mesh_db->putVector<int>( "Size", { 20, 40 } );
    mesh_db->putVector<double>( "Range", { 0.004095, 0, 0.01016 } );
    auto params = std::make_shared<AMP::Mesh::MeshParameters>( mesh_db );
    params->setComm( globalComm );

    // Create the meshes from the input database
    auto mesh = AMP::Mesh::Mesh::buildMesh( params );

    // Create a simple DOFManager
    auto nodalVariable = std::make_shared<AMP::LinearAlgebra::Variable>( "test" );
    auto DOFs =
        AMP::Discretization::simpleDOFManager::create( mesh, AMP::Mesh::GeomType::Vertex, 1, 1 );

    // Create the vectors
    auto v1 = createVector( DOFs, nodalVariable );
    auto v2 = createVector( DOFs, nodalVariable );

    // Initialize the vectors
    v1->setToScalar( 0.0 );
    v2->setToScalar( 0.0 );

    // Time makeConsistentSet
    globalComm.barrier();
    double start_time = AMP::AMP_MPI::time();
    v1->makeConsistent( AMP::LinearAlgebra::VectorData::ScatterType::CONSISTENT_SET );
    globalComm.barrier();
    double end_time = AMP::AMP_MPI::time();

    // Get the number of ghost values in the communication list
    auto commList   = v1->getCommunicationList();
    auto ghost_ids  = commList->getGhostIDList();
    size_t N_global = v1->getGlobalSize();
    size_t N_ghosts = globalComm.sumReduce( ghost_ids.size() );
    size_t N_ghosts2 =
        globalComm.sumReduce( mesh->numGhostElements( AMP::Mesh::GeomType::Vertex, 1 ) );

    // Print the results
    if ( globalComm.getRank() == 0 ) {
        std::cout << "Time for makeConsistent: " << end_time - start_time << std::endl;
        std::cout << "There are " << N_global << " global values" << std::endl;
        std::cout << "There are " << N_ghosts << " global ghost values" << std::endl;
        std::cout << "There are " << N_ghosts2 << " global ghost values in the iterator"
                  << std::endl;
    }
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
