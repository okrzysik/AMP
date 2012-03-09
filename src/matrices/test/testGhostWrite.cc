
#include "utils/AMPManager.h"
#include "utils/UnitTest.h"
#include "utils/Utilities.h"
#include <iostream>
#include <string>

#include "boost/shared_ptr.hpp"

#include "utils/Database.h"
#include "utils/InputDatabase.h"
#include "utils/InputManager.h"
#include "utils/AMP_MPI.h"
#include "utils/AMPManager.h"
#include "utils/PIO.h"

#include "ampmesh/Mesh.h"

#include "test_Matrix.h"

void myTest(AMP::UnitTest *ut)
{
    std::string exeName("testGhostWrite");
    std::string log_file = "output_" + exeName;
    AMP::PIO::logOnlyNodeZero(log_file);

    // Create the mesh parameter object
    boost::shared_ptr<AMP::MemoryDatabase> database(new AMP::MemoryDatabase("Mesh"));
    database->putString("MeshType","libMesh");
    database->putString("MeshName","mesh");
    database->putString("FileName","cube8.e");
    database->putInteger("dim",3);
    boost::shared_ptr<AMP::Mesh::MeshParameters> params(new AMP::Mesh::MeshParameters(database));
    params->setComm(AMP::AMP_MPI(AMP_COMM_WORLD));

    // Create the mesh
    AMP::Mesh::Mesh::shared_ptr mesh = AMP::Mesh::Mesh::buildMesh( params );

    // Create the matrix
    AMP::Discretization::DOFManager::shared_ptr  DOFs = AMP::Discretization::simpleDOFManager::create(mesh,AMP::Mesh::Vertex,1,1);
    AMP::LinearAlgebra::Variable::shared_ptr var1( new AMP::LinearAlgebra::Variable("a") );
    AMP::LinearAlgebra::Variable::shared_ptr var2( new AMP::LinearAlgebra::Variable("b") );
    AMP::LinearAlgebra::Vector::shared_ptr  vec1 = AMP::LinearAlgebra::createVector( DOFs, var1 );
    AMP::LinearAlgebra::Vector::shared_ptr  vec2 = AMP::LinearAlgebra::createVector( DOFs, var2 );
    AMP::LinearAlgebra::Matrix::shared_ptr  matrix = AMP::LinearAlgebra::createMatrix ( vec1, vec2 );

    // Loop through the owned and ghost nodes
    AMP::Mesh::MeshIterator el = mesh->getIterator(AMP::Mesh::Volume,0);
    AMP::Mesh::MeshIterator end_el = el.end();
    for( ; el != end_el; ++el) {
        // Get the DOFs for all nodes
        std::vector<size_t> dofs;
        std::vector<size_t> dofIndices;
        std::vector<AMP::Mesh::MeshElement> elements = el->getElements(AMP::Mesh::Vertex);
        for (size_t i=0; i<elements.size(); i++) {
            DOFs->getDOFs( elements[i].globalID(), dofs );
            for (size_t j=0; j<dofs.size(); j++) {
                dofIndices.push_back(dofs[j]);
            }
        }
        for (size_t j=0; j<dofIndices.size(); j++) {
            for (size_t i=0; i<dofIndices.size(); i++) {
                matrix->setValueByGlobalID ( dofIndices[j], dofIndices[i], 1.0 );
            }//end for i
        }//end for j
    }//end for el

    matrix->makeConsistent();
}


int main(int argc, char *argv[])
{
  AMP::AMPManager::startup(argc, argv);
  AMP::UnitTest ut;

  try {
    myTest(&ut);
  } catch (std::exception &err) {
    std::cout << "ERROR: While testing "<<argv[0] << err.what() << std::endl;
    ut.failure("ERROR: While testing");
  } catch( ... ) {
    std::cout << "ERROR: While testing "<<argv[0] << "An unknown exception was thrown." << std::endl;
    ut.failure("ERROR: While testing");
  }

  ut.report();

  int num_failed = ut.NumFailGlobal();
  AMP::AMPManager::shutdown();
  return num_failed;
}   


