
#include <utils/UnitTest.h>
#include <utils/Utilities.h>
#include <utils/shared_ptr.h>
#include <utils/Database.h>
#include <utils/InputDatabase.h>
#include <utils/InputManager.h>
#include <utils/AMP_MPI.h>
#include <utils/AMPManager.h>
#include <utils/PIO.h>

#include <vectors/VectorBuilder.h>

#include <discretization/simpleDOF_Manager.h>

#include <ampmesh/Mesh.h>

#include <operators/map/dtk/DTKAMPVectorHelpers.h>
#include <operators/map/dtk/DTKAMPMeshManager.h>

#include <DTK_ConsistentInterpolationOperator.hpp>

#include <iostream>
#include <string>
#include <cstdlib>

double testFunction1( const std::vector<double>& coords )
{
    return coords[0] * coords[1] * coords[2];
}

void myTest(AMP::UnitTest *ut)
{
    std::string exeName("testDTKInterpolation");
    std::string log_file = "output_" + exeName;
    std::string msgPrefix;
    AMP::PIO::logOnlyNodeZero(log_file);

    // load the source mesh
    AMP::pout<<"Loading the source mesh"<<std::endl;
    boost::shared_ptr<AMP::InputDatabase> input_db(new AMP::InputDatabase("input_db"));
    AMP::AMP_MPI globalComm(AMP_COMM_WORLD);

    std::string input_file = "input_" + exeName;
    AMP::InputManager::getManager()->parseInputFile(input_file, input_db);
    input_db->printClassData(AMP::plog);

    boost::shared_ptr<AMP::Database> sourceMeshDatabase = input_db->getDatabase("SourceMesh");
    boost::shared_ptr<AMP::Mesh::MeshParameters> sourceMeshParams(new AMP::Mesh::MeshParameters(sourceMeshDatabase));
    sourceMeshParams->setComm(AMP::AMP_MPI(AMP_COMM_WORLD));
    AMP::Mesh::Mesh::shared_ptr sourceMesh = AMP::Mesh::Mesh::buildMesh(sourceMeshParams);
    std::size_t const numVerticesOnSourceMesh = sourceMesh->numGlobalElements(AMP::Mesh::Vertex);
    std::size_t const numElementsOnSourceMesh = sourceMesh->numGlobalElements(AMP::Mesh::Volume);
    AMP::pout<<"source mesh contains "<<numVerticesOnSourceMesh<<" vertices\n";
    AMP::pout<<"source mesh contains "<<numElementsOnSourceMesh<<" elements\n";

    // build source vector
    AMP::pout<<"Building the source vector"<<std::endl;
    bool const split = true;
    int const ghostWidth = 0;
    int const dofsPerNode = 1;
    AMP::Discretization::DOFManager::shared_ptr sourceDofManager = AMP::Discretization::simpleDOFManager::create(sourceMesh, AMP::Mesh::Vertex, ghostWidth, dofsPerNode);
    AMP::LinearAlgebra::Variable::shared_ptr variable(new AMP::LinearAlgebra::Variable("dummy"));
    AMP::LinearAlgebra::Vector::shared_ptr sourceVector = AMP::LinearAlgebra::createVector(sourceDofManager, variable, split);
#ifdef USE_EXT_SILO
    AMP::Utilities::Writer::shared_ptr siloWriter = AMP::Utilities::Writer::buildWriter("silo");
    siloWriter->setDecomposition(1);
    siloWriter->registerVector(sourceVector, sourceMesh, AMP::Mesh::Vertex, "vector");
    siloWriter->writeFile("source", 0);
#endif

    AMP::Mesh::MeshIterator sourceMeshIterator = sourceMesh->getIterator();
    for ( sourceMeshIterator = sourceMeshIterator.begin();
          sourceMeshIterator != sourceMeshIterator.end();
          ++sourceMeshIterator ) {
        souceDofManager->getDOFs(sourceMeshIterator->globalID(), dofIndices);
        
        AMP_ASSERT( static_cast<std::size_t>(value) == (2 * dofIndices[0]) );
    }

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


