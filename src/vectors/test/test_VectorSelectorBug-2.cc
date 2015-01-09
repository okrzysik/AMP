#include <utils/AMPManager.h>
#include <utils/UnitTest.h>
#include <utils/InputDatabase.h>
#include <utils/InputManager.h>
#include <ampmesh/Mesh.h>
#include <discretization/simpleDOF_Manager.h>
#include <vectors/VectorBuilder.h>
#include <vectors/VectorSelector.h>

void myTest(AMP::UnitTest *ut)
{
    std::string const exeName = "test_VectorSelectorBug-2";
    AMP::AMP_MPI globalComm(AMP_COMM_WORLD);

    // Parse input file
    std::string const inputFile = "input_" + exeName;
    AMP::shared_ptr<AMP::InputDatabase> inputDatabase(new AMP::InputDatabase("input_db"));
    AMP::InputManager::getManager()->parseInputFile(inputFile, inputDatabase);
    
    // Read the mesh
    AMP::pout<<"--------------------\n";
    AMP::pout<<"    LOADING MESH    \n";
    AMP::pout<<"--------------------\n";
    AMP::Database::shared_ptr meshDatabase = inputDatabase->getDatabase("Mesh");
    AMP::Mesh::MeshParameters::shared_ptr meshParams(new AMP::Mesh::MeshParameters(meshDatabase));
    meshParams->setComm(globalComm);
    AMP::Mesh::Mesh::shared_ptr mesh = AMP::Mesh::Mesh::buildMesh(meshParams);
    
    int const boundaryID = 2;
    AMP::Mesh::Mesh::shared_ptr boundaryMesh = 
        mesh->Subset(mesh->getBoundaryIDIterator(AMP::Mesh::Volume, boundaryID));
    AMP::pout<<"mesh volumes = "<<mesh->numGlobalElements(AMP::Mesh::Volume)<<"\n";
    AMP::pout<<"boundary mesh volumes = "<<boundaryMesh->numGlobalElements(AMP::Mesh::Volume)<<"\n";
    
    // make a simple vector 
    bool const split = true;
    int const ghostWidth = 0;
    int const dofsPerNode = 1;
    AMP::Discretization::DOFManager::shared_ptr dofManager = 
        AMP::Discretization::simpleDOFManager::create(mesh, AMP::Mesh::Vertex, ghostWidth, dofsPerNode);
    AMP::LinearAlgebra::Variable::shared_ptr variable(new AMP::LinearAlgebra::Variable("noname"));
    AMP::LinearAlgebra::Vector::shared_ptr vector = AMP::LinearAlgebra::createVector(dofManager, variable, split);

    
    AMP::pout<<"mesh vector size "         <<vector->select(AMP::LinearAlgebra::VS_Mesh(mesh)        , "var")->getGlobalSize()<<"\n";
    AMP::pout<<"boundary mesh vector size "<<vector->select(AMP::LinearAlgebra::VS_Mesh(boundaryMesh), "var")->getGlobalSize()<<"\n";

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
