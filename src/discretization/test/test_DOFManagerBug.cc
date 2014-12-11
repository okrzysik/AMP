
#include <utils/UnitTest.h>
#include <utils/Utilities.h>
#include <utils/shared_ptr.h>
#include <utils/Database.h>
#include <utils/InputDatabase.h>
#include <utils/InputManager.h>
#include <utils/AMP_MPI.h>
#include <utils/AMPManager.h>
#include <utils/PIO.h>
#include <discretization/simpleDOF_Manager.h>
#include <ampmesh/Mesh.h>
#include <vectors/VectorBuilder.h>

void myTest(AMP::UnitTest *ut)
{
    std::string exeName("DOFManagerBug");
    std::string log_file = "output_" + exeName;
    std::string msgPrefix;
    AMP::PIO::logOnlyNodeZero(log_file);

    boost::shared_ptr<AMP::InputDatabase> input_db(new AMP::InputDatabase("input_db"));
    AMP::AMP_MPI globalComm(AMP_COMM_WORLD);

    std::string input_file = "input_" + exeName;
    AMP::InputManager::getManager()->parseInputFile(input_file, input_db);
    input_db->printClassData(AMP::plog);

    // load the mesh
    boost::shared_ptr<AMP::Database> meshDatabase = input_db->getDatabase("Mesh");
    boost::shared_ptr<AMP::Mesh::MeshParameters> meshParams(new AMP::Mesh::MeshParameters(meshDatabase));
    meshParams->setComm(AMP::AMP_MPI(AMP_COMM_WORLD));
    AMP::Mesh::Mesh::shared_ptr mesh = AMP::Mesh::Mesh::buildMesh(meshParams);

    bool const split = true;
    int const ghostWidth = 1;
    int const dofsPerNode = 1;
    AMP::Discretization::DOFManager::shared_ptr	dofManager = 
        AMP::Discretization::simpleDOFManager::create(mesh, AMP::Mesh::Vertex, ghostWidth, dofsPerNode);

    AMP::LinearAlgebra::Variable::shared_ptr variable(new AMP::LinearAlgebra::Variable("var"));
    AMP::LinearAlgebra::Vector::shared_ptr ampVector = 
        AMP::LinearAlgebra::createVector(dofManager, variable, split);

    AMP::Mesh::MeshIterator meshIterator = mesh->getIterator(AMP::Mesh::Vertex, ghostWidth);

    AMP_ASSERT( meshIterator.size() == ampVector->getDOFManager()->getIterator().size() );


    ut->passes( exeName );
}


int main(int argc, char *argv[])
{
    AMP::AMPManagerProperties startup_properties;
    startup_properties.use_MPI_Abort = false;
    AMP::AMPManager::startup(argc,argv,startup_properties);
    AMP::UnitTest ut;

    myTest(&ut);

    ut.report();

    int num_failed = ut.NumFailGlobal();
    AMP::AMPManager::shutdown();
    return num_failed;
}   


