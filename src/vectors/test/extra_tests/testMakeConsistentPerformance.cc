#include "utils/Database.h"
#include "utils/InputDatabase.h"
#include "utils/InputManager.h"
#include "utils/AMP_MPI.h"
#include "utils/AMPManager.h"
#include "utils/UnitTest.h"
#include "utils/Utilities.h"
#include "utils/PIO.h"
#include "ampmesh/Mesh.h"
#include "discretization/DOF_Manager.h"
#include "discretization/simpleDOF_Manager.h"
#include "discretization/NodalVariable.h"
#include "operators/map/NodeToNodeMap.h"
#include "operators/map/AsyncMapColumnOperator.h"




void  runTest ( AMP::UnitTest *ut )
{

    // Get the Mesh database and create the mesh parameters
    AMP::AMP_MPI globalComm(AMP_COMM_WORLD);
    boost::shared_ptr<AMP::MemoryDatabase> mesh_db(new AMP::MemoryDatabase("Mesh"));
    mesh_db->putInteger("dim",3);
    mesh_db->putString("MeshName","mesh1");
    mesh_db->putString("MeshType","libMesh");
    mesh_db->putString("FileName","pellet_1x.e");
    boost::shared_ptr<AMP::Mesh::MeshParameters> params(new AMP::Mesh::MeshParameters(mesh_db));
    params->setComm(globalComm);

    // Create the meshes from the input database
    boost::shared_ptr<AMP::Mesh::Mesh> mesh = AMP::Mesh::Mesh::buildMesh(params);

    // Create a simple DOFManager
    int DOFsPerNode = 3;
    std::string varName = "test";
    AMP::LinearAlgebra::Variable::shared_ptr nodalVariable( new AMP::Discretization::NodalVariable(DOFsPerNode,varName) );
    AMP::Discretization::DOFManagerParameters::shared_ptr DOFparams( new AMP::Discretization::DOFManagerParameters(mesh) );
    boost::shared_ptr<AMP::Discretization::simpleDOFManager> DOFs( new AMP::Discretization::simpleDOFManager(mesh,AMP::Mesh::Vertex,1,DOFsPerNode) );

    // Create the vectors
    AMP::LinearAlgebra::Vector::shared_ptr  dummy;
    AMP::LinearAlgebra::Vector::shared_ptr v1 = DOFs->createVector( nodalVariable );
    AMP::LinearAlgebra::Vector::shared_ptr v2 = DOFs->createVector( nodalVariable );

    // Initialize the vectors
    v1->setToScalar(0.0);
    v2->setToScalar(0.0);

    // Time makeConsistentSet
    globalComm.barrier();
    double start_time = AMP::AMP_MPI::time();
    v1->makeConsistent ( AMP::LinearAlgebra::Vector::CONSISTENT_SET );
    globalComm.barrier();
    double end_time = AMP::AMP_MPI::time();
    std::cout << "Time for makeConsistent: " << end_time-start_time << std::endl;
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

