#include <utils/AMPManager.h>
#include <utils/UnitTest.h>
#include <utils/InputDatabase.h>
#include <ampmesh/mesh.h>
#include <discretization/simpleDOF_Manager.h>
#include <vectors/VectorBuilder.h>
#include <vectors/VectorSelector.h>

void myTest(AMP::UnitTest *ut)
{
    AMP::AMP_MPI globalComm(AMP_COMM_WORLD);

    // build mesh
    AMP::Database::shared_ptr meshDatabase(new AMP::InputDatabase("Mesh"));
    meshDatabase->putString("MeshName", "NoName");
    meshDatabase->putString("MeshType", "AMP");
    meshDatabase->putInteger("dim", 3);
    meshDatabase->putDouble("x_offset", 0.0);
    meshDatabase->putDouble("y_offset", 0.0);
    meshDatabase->putDouble("z_offset", 0.0);
    meshDatabase->putString("Generator", "cube");
    std::vector<int> size(3, 5);
    meshDatabase->putIntegerArray("Size", size);
    std::vector<double> range(6, 0.0);
    range[1] = range[3] = range[5] = 1.0;
    meshDatabase->putDoubleArray("Range", range);
    AMP::Mesh::MeshParameters::shared_ptr meshParams(new AMP::Mesh::MeshParameters(meshDatabase));
    meshParams->setComm(globalComm);
    AMP::Mesh::Mesh::shared_ptr mesh = AMP::Mesh::Mesh::buildMesh(meshParams);

    // make a simple vector 
    bool const split = true;
    int const ghostWidth = 0;
    int const dofsPerNode = 1;
    AMP::Discretization::DOFManager::shared_ptr dofManager = 
        AMP::Discretization::simpleDOFManager::create(mesh, AMP::Mesh::Vertex, ghostWidth, dofsPerNode);
    AMP::LinearAlgebra::Variable::shared_ptr variable(new AMP::LinearAlgebra::Variable("noname"));
    AMP::LinearAlgebra::Vector::shared_ptr vector = AMP::LinearAlgebra::createVector(dofManager, variable, split);

    // select from iterator 
    AMP::Mesh::MeshIterator meshIterator = mesh->getSurfaceIterator(AMP::Mesh::Vertex, ghostWidth);
    AMP::LinearAlgebra::Vector::const_shared_ptr boundaryVector
        = vector->constSelect(AMP::LinearAlgebra::VS_MeshIterator(meshIterator, globalComm), "dummy");
    AMP::pout<<"this works but look when iterator is set to end"<<std::endl;
    meshIterator = meshIterator.end();
    boundaryVector
        = vector->constSelect(AMP::LinearAlgebra::VS_MeshIterator(meshIterator, globalComm), "dummy");
}

int main(int argc, char *argv[])
{
    AMP::AMPManager::startup(argc, argv);
    AMP::UnitTest ut;
  
    try {
        myTest(&ut);
    } catch (std::exception &err) {
        AMP::pout << "ERROR: While testing "<<argv[0] << err.what() << std::endl;
        ut.failure("ERROR: While testing");
    } catch( ... ) {
        AMP::pout << "ERROR: While testing "<<argv[0] << "An unknown exception was thrown." << std::endl;
        ut.failure("ERROR: While testing");
    }
  
    ut.report();
    int num_failed = ut.NumFailGlobal();
  
    AMP::AMPManager::shutdown();
    return num_failed;
}  
