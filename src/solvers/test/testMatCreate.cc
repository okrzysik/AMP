#include <iostream>
#include <string>

#include "utils/InputManager.h"
#include "utils/AMPManager.h"
#include "utils/UnitTest.h"
#include "utils/Utilities.h"
#include "utils/ReadTestMesh.h"

#include "ampmesh/Mesh.h"
#include "discretization/DOF_Manager.h"
#include "discretization/simpleDOF_Manager.h"
#include "vectors/Variable.h"
#include "vectors/Vector.h"
#include "vectors/VectorBuilder.h"
#include "matrices/MatrixBuilder.h"


void myTest(AMP::UnitTest *ut, std::string input_file) {

    std::string log_file = "output_testMatCreate";
    AMP::PIO::logOnlyNodeZero(log_file);

    // Read the input file
    boost::shared_ptr<AMP::InputDatabase>  input_db ( new AMP::InputDatabase ( "input_db" ) );
    AMP::InputManager::getManager()->parseInputFile ( input_file , input_db );

    // Get the Mesh database and create the mesh parameters
    boost::shared_ptr<AMP::Database> database = input_db->getDatabase( "Mesh" );
    boost::shared_ptr<AMP::Mesh::MeshParameters> params(new AMP::Mesh::MeshParameters(database));
    params->setComm(AMP::AMP_MPI(AMP_COMM_WORLD));

    // Create the meshes from the input database
    AMP::Mesh::Mesh::shared_ptr mesh = AMP::Mesh::Mesh::buildMesh(params);

    // Create the DOF manager
    AMP::Discretization::DOFManager::shared_ptr scalarDOFs = AMP::Discretization::simpleDOFManager::create(mesh,AMP::Mesh::Vertex,1,1);
    AMP::Discretization::DOFManager::shared_ptr vectorDOFs = AMP::Discretization::simpleDOFManager::create(mesh,AMP::Mesh::Vertex,1,3);

    // Create the vectors
    AMP::LinearAlgebra::Variable::shared_ptr inVar (new AMP::LinearAlgebra::Variable("inputVar"));
    AMP::LinearAlgebra::Variable::shared_ptr outVar (new AMP::LinearAlgebra::Variable("outputVar"));
    AMP::LinearAlgebra::Vector::shared_ptr inVec = AMP::LinearAlgebra::createVector( vectorDOFs, inVar );
    AMP::LinearAlgebra::Vector::shared_ptr outVec = AMP::LinearAlgebra::createVector( scalarDOFs, outVar );

    // Create the matrix
    AMP::LinearAlgebra::Matrix::shared_ptr mat1 = AMP::LinearAlgebra::createMatrix ( inVec, outVec );
    if(mat1.get()!=NULL){
        ut->passes("Able to create a non-square matrices");
    } else {
        ut->failure("Unable to create a non-square matrices");
    }

}


int main(int argc, char *argv[])
{
  AMP::AMPManager::startup(argc, argv);
  AMP::UnitTest ut;

  try {
    myTest(&ut, "input_testLinearFlow-1");
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



