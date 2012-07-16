
#include "utils/AMPManager.h"
#include "utils/UnitTest.h"
#include "utils/Utilities.h"
#include "utils/Database.h"
#include "utils/InputDatabase.h"
#include "utils/InputManager.h"
#include "utils/AMP_MPI.h"
#include "utils/PIO.h"

#include "discretization/DOF_Manager.h"
#include "discretization/simpleDOF_Manager.h"
#include "vectors/Variable.h"
#include "vectors/Vector.h"
#include "vectors/VectorBuilder.h"

#include "externVars.h"

#include "ampmesh/Mesh.h"
#include "ampmesh/dendro/DendroSearch.h"

#include "operators/contact/NodeToSegmentConstraintsOperator.h"

#include <fstream>
#include <boost/lexical_cast.hpp>


void myTest(AMP::UnitTest *ut, std::string exeName) {
  std::string input_file = "input_" + exeName;
  std::string log_file = "output_" + exeName; 

  AMP::PIO::logOnlyNodeZero(log_file);
  AMP::AMP_MPI globalComm(AMP_COMM_WORLD);

#ifdef USE_SILO
//  AMP::Mesh::SiloIO::shared_ptr siloWriter(new AMP::Mesh::SiloIO);
#endif

  int npes = globalComm.getSize();
  int rank = globalComm.getRank();

  // Load the input file
  globalComm.barrier();
  double inpReadBeginTime = MPI_Wtime();
  boost::shared_ptr<AMP::InputDatabase> input_db(new AMP::InputDatabase("input_db"));
  AMP::InputManager::getManager()->parseInputFile(input_file, input_db);
  input_db->printClassData(AMP::plog);
  globalComm.barrier();
  double inpReadEndTime = MPI_Wtime();
  if(!rank) {
    std::cout<<"Finished parsing the input file in "<<(inpReadEndTime - inpReadBeginTime)<<" seconds."<<std::endl;
  }

  // Load the mesh
  globalComm.barrier();
  double meshBeginTime = MPI_Wtime();
  AMP_INSIST(input_db->keyExists("Mesh"), "Key ''Mesh'' is missing!");
  boost::shared_ptr<AMP::Database> mesh_db = input_db->getDatabase("Mesh");
  boost::shared_ptr<AMP::Mesh::MeshParameters> meshParams(new AMP::Mesh::MeshParameters(mesh_db));
  meshParams->setComm(globalComm);
  AMP::Mesh::Mesh::shared_ptr meshAdapter = AMP::Mesh::Mesh::buildMesh(meshParams);
  globalComm.barrier();
  double meshEndTime = MPI_Wtime();
  if(!rank) {
    std::cout<<"Finished reading the mesh in "<<(meshEndTime - meshBeginTime)<<" seconds."<<std::endl;
  }

  // build the contact operator
  boost::shared_ptr<AMP::Operator::NodeToSegmentConstraintsOperatorParameters> 
      nodeToSegmentConstraintsOperatorParams( new AMP::Operator::NodeToSegmentConstraintsOperatorParameters(mesh_db) );

  std::vector<AMP::Mesh::MeshID> meshIDs = meshAdapter->getBaseMeshIDs();
  nodeToSegmentConstraintsOperatorParams->d_MasterMeshID = meshIDs[0];
  nodeToSegmentConstraintsOperatorParams->d_SlaveMeshID = meshIDs[1];
  nodeToSegmentConstraintsOperatorParams->d_MasterBoundaryID = 1;
  nodeToSegmentConstraintsOperatorParams->d_SlaveBoundaryID = 2;
  
  int dofsPerNode = 3;
  int nodalGhostWidth = 1;
  bool split = true;
  AMP::Discretization::DOFManager::shared_ptr dofManager = AMP::Discretization::simpleDOFManager::create(meshAdapter,
      AMP::Mesh::Vertex, nodalGhostWidth, dofsPerNode, split);
  nodeToSegmentConstraintsOperatorParams->d_DOFsPerNode = dofsPerNode;
  nodeToSegmentConstraintsOperatorParams->d_DOFManager = dofManager;

  nodeToSegmentConstraintsOperatorParams->d_GlobalComm = globalComm;
  nodeToSegmentConstraintsOperatorParams->d_Mesh = meshAdapter;

  boost::shared_ptr<AMP::Operator::NodeToSegmentConstraintsOperator> 
      nodeToSegmentConstraintsOperator( new AMP::Operator::NodeToSegmentConstraintsOperator(nodeToSegmentConstraintsOperatorParams) );

  nodeToSegmentConstraintsOperator->reset(nodeToSegmentConstraintsOperatorParams);

  AMP::LinearAlgebra::Variable::shared_ptr dummyVariable(new AMP::LinearAlgebra::Variable("Dummy"));
  AMP::LinearAlgebra::Vector::shared_ptr dummyInVector = createVector(dofManager, dummyVariable, split);
  AMP::LinearAlgebra::Vector::shared_ptr dummyOutVector = createVector(dofManager, dummyVariable, split);

  nodeToSegmentConstraintsOperator->apply(dummyInVector, dummyInVector, dummyOutVector);
  nodeToSegmentConstraintsOperator->applyTranspose(dummyInVector, dummyInVector, dummyOutVector);

  

#ifdef USE_SILO
//  AMP::Linear
//  siloWriter->registerVector(mechNlSolVec, meshAdapter, AMP::Mesh::Vertex, "Solution");
//  char outFileName[256];
//  sprintf(outFileName, "LoadPrescribed-DeformedPlateWithHole-LinearElasticity_%d", step);
//  siloWriter->writeFile(outFileName, 0);
#endif

  ut->passes(exeName);
}


int main(int argc, char *argv[])
{
  AMP::AMPManager::startup(argc, argv);
  AMP::UnitTest ut;

  std::string exeName = "testNodeToSegmentConstraintsOperator";

  try {
    myTest(&ut, exeName);
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



