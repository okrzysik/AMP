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
    AMP::Database::shared_ptr meshDatabase = inputDatabase->getDatabase("Mesh");
    AMP::Mesh::MeshParameters::shared_ptr meshParams(new AMP::Mesh::MeshParameters(meshDatabase));
    meshParams->setComm(globalComm);
    AMP::Mesh::Mesh::shared_ptr mesh = AMP::Mesh::Mesh::buildMesh(meshParams);
    
    // Subset the mesh
    int const boundaryID = inputDatabase->getInteger("BoundaryID");
    AMP::Mesh::Mesh::shared_ptr boundaryMesh = 
        mesh->Subset(mesh->getBoundaryIDIterator(AMP::Mesh::Volume, boundaryID));

    // Check mesh sizes
    AMP::pout<<"mesh volumes = "<<mesh->numGlobalElements(AMP::Mesh::Volume)<<"\n";
    AMP::pout<<"boundary mesh volumes = "<<boundaryMesh->numGlobalElements(AMP::Mesh::Volume)<<"\n";
    AMP::pout<<"mesh vertices = "<<mesh->numGlobalElements(AMP::Mesh::Vertex)<<"\n";
    AMP::pout<<"boundary mesh vertices = "<<boundaryMesh->numGlobalElements(AMP::Mesh::Vertex)<<"\n";

    std::size_t const expectedSize = boundaryMesh->numGlobalElements(AMP::Mesh::Vertex);
    
    // Build a vector 
    bool const split = true;
    int const ghostWidth = 0;
    int const dofsPerNode = 1;
    AMP::Discretization::DOFManager::shared_ptr dofManager = 
        AMP::Discretization::simpleDOFManager::create(mesh, AMP::Mesh::Vertex, ghostWidth, dofsPerNode);
    AMP::LinearAlgebra::Variable::shared_ptr variable(new AMP::LinearAlgebra::Variable("var"));
    AMP::LinearAlgebra::Vector::shared_ptr vector = 
        AMP::LinearAlgebra::createVector(dofManager, variable, split);

    // Subset the vector
    AMP::LinearAlgebra::Vector::shared_ptr boundaryVector = 
        vector->select(AMP::LinearAlgebra::VS_Mesh(boundaryMesh), "var");

    // Check vector sizes
    AMP::pout<<"vector size "         <<vector        ->getGlobalSize()<<"\n";
    AMP::pout<<"boundary vector size "<<boundaryVector->getGlobalSize()<<"\n";

    std::size_t const actualSize = boundaryVector->getGlobalSize();

    AMP_ASSERT(expectedSize == actualSize);

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
