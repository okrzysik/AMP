
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


double dummyFunction(const std::vector<double> &xyz, const int dof) {
  AMP_ASSERT(xyz.size() == 3);
  double x = xyz[0], y = xyz[1], z = xyz[2];
  //  return 7.0;
  //return (1.0 + 6.0 * x) * (2.0 - 5.0 * y) * (3.0 + 4.0 * z);
  return (1.0 + 6.0 * x) + (2.0 - 5.0 * y) + (3.0 + 4.0 * z);
}

void myTest(AMP::UnitTest *ut, std::string exeName) {
  std::string input_file = "input_" + exeName;
  std::string log_file = "output_" + exeName; 

  ot::RegisterEvents();

  AMP::PIO::logOnlyNodeZero(log_file);
  AMP::AMP_MPI globalComm(AMP_COMM_WORLD);

  int rank = globalComm.getRank();
  int npes = globalComm.getSize();

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
  meshParams->setComm(AMP::AMP_MPI(AMP_COMM_WORLD));
  AMP::Mesh::Mesh::shared_ptr meshAdapter = AMP::Mesh::Mesh::buildMesh(meshParams);
  globalComm.barrier();
  double meshEndTime = MPI_Wtime();
  if(!rank) {
    std::cout<<"Finished reading the mesh in "<<(meshEndTime - meshBeginTime)<<" seconds."<<std::endl;
  }

  // Create a vector field
  int DOFsPerNode = 1;
  int nodalGhostWidth = 1;
  bool split = true;
  AMP::Discretization::DOFManager::shared_ptr DOFs = AMP::Discretization::simpleDOFManager::create(meshAdapter,
      AMP::Mesh::Vertex, nodalGhostWidth, DOFsPerNode, split);
  AMP::LinearAlgebra::Variable::shared_ptr dummyVariable(new AMP::LinearAlgebra::Variable("Dummy"));
  AMP::LinearAlgebra::Vector::shared_ptr dummyVector = createVector(DOFs, dummyVariable, split);

  AMP::Mesh::MeshIterator node = meshAdapter->getIterator(AMP::Mesh::Vertex, 0);
  AMP::Mesh::MeshIterator end_node = node.end();
  for ( ; node != end_node; ++node) {
    std::vector<size_t> globalID;
    DOFs->getDOFs(node->globalID(), globalID); 
    for(int d = 0; d < globalID.size(); ++d) {
      dummyVector->setLocalValueByGlobalID(globalID[d], dummyFunction(node->coord(), d));
    }//end d
  }

  double minCoords[3];
  double maxCoords[3];
  std::vector<double> box = meshAdapter->getBoundingBox();
  for(int i = 0; i < meshAdapter->getDim(); ++i) {
    minCoords[i] = box[2*i+0];
    maxCoords[i] = box[2*i+1];
  }

  int totalNumPts = input_db->getInteger("TotalNumberOfPoints");
  int avgNumPts = totalNumPts/npes;
  int extraNumPts = totalNumPts%npes;

  int numLocalPts = avgNumPts;
  if(rank < extraNumPts) {
    numLocalPts++;
  }

  // Generate Random points in [min, max]
  const unsigned int seed = (0x1234567 + (24135*rank));
  srand48(seed);

  std::vector<double> pts(3*numLocalPts);
  for(int i = 0; i < numLocalPts; ++i) {
    double x = ((maxCoords[0] - minCoords[0])*drand48()) + minCoords[0];
    double y = ((maxCoords[1] - minCoords[1])*drand48()) + minCoords[1];
    double z = ((maxCoords[2] - minCoords[2])*drand48()) + minCoords[2];
    pts[3*i] = x;
    pts[(3*i) + 1] = y;
    pts[(3*i) + 2] = z;
  }//end i
  if(!rank) {
    std::cout<<"Finished generating "<<totalNumPts <<" random points for search!"<<std::endl;
  }

  DendroSearch dendroSearch(globalComm, meshAdapter);
  std::vector<double> interpolatedData; 
  std::vector<bool> interpolationWasDone;
  dendroSearch.searchAndInterpolate(dummyVector, DOFsPerNode, pts, interpolatedData, interpolationWasDone);
  AMP_ASSERT(interpolatedData.size() == (DOFsPerNode*numLocalPts));
  AMP_ASSERT(interpolationWasDone.size() == numLocalPts);

  int localNotFound = 0;
  if(numLocalPts > 0) {
    localNotFound = static_cast<int>(std::count(interpolationWasDone.begin(), interpolationWasDone.end(), false));
  }
  int globalNotFound;
  MPI_Allreduce(&localNotFound, &globalNotFound, 1, MPI_INT, MPI_SUM, globalComm.getCommunicator());
  if(!rank) {
    std::cout<<globalNotFound<<" points (total) were not found"<<std::endl;
  }

  std::vector<double> interpolationError((DOFsPerNode*numLocalPts), 0.0);
  for(unsigned int i = 0; i < numLocalPts; ++i) {
    if(interpolationWasDone[i]) {
      for(int d = 0; d < DOFsPerNode; ++d) {
        interpolationError[(i*DOFsPerNode) + d] = fabs(interpolatedData[(i*DOFsPerNode) + d] -
            dummyFunction(std::vector<double>(&(pts[3*i]), &(pts[3*i]) + 3), d));
      }//end d
    }
  } // end for i
  double localMaxError = 0;
  if(numLocalPts > 0) {
    localMaxError = *(std::max_element(interpolationError.begin(), interpolationError.end()));
  }
  if(!rank) {
    std::cout<<"Finished computing the local squared norm of the interpolation error."<<std::endl;
  }
  globalComm.barrier();

  double globalMaxError = globalComm.maxReduce<double>(localMaxError);
  if(!rank) {
    std::cout<<"Global max error is "<<std::setprecision(15)<<globalMaxError<<std::endl;
  }

  AMP_ASSERT(globalMaxError < 1.0e-12);

  std::vector<AMP::Mesh::MeshElementID> faceVerticesGlobalIDs;
  std::vector<double> shiftGlobalCoords, projectionLocalCoordsOnFace;
  std::vector<int> flags;
  dendroSearch.projectOnBoundaryID(4, faceVerticesGlobalIDs, shiftGlobalCoords, projectionLocalCoordsOnFace, flags);
  globalComm.barrier();
  if(!rank) {
    std::cout<<"PAR ICI"<<std::endl;
    for (size_t i = 0; i < 100; ++i) { std::cout<<std::flush; }
  }
  globalComm.barrier();
  unsigned int localPtsNotFound = std::count(flags.begin(), flags.end(), DendroSearch::NotFound);
  unsigned int localPtsFoundNotOnBoundary = std::count(flags.begin(), flags.end(), DendroSearch::FoundNotOnBoundary);
  unsigned int localPtsFoundOnBoundary = std::count(flags.begin(), flags.end(), DendroSearch::FoundOnBoundary);
  unsigned int globalPtsNotFound = globalComm.sumReduce(localPtsNotFound);
  unsigned int globalPtsFoundNotOnBoundary = globalComm.sumReduce(localPtsFoundNotOnBoundary);
  unsigned int globalPtsFoundOnBoundary = globalComm.sumReduce(localPtsFoundOnBoundary);
  if(!rank) {
    std::cout<<"Global number of points not found is "<<globalPtsNotFound<<std::endl;
    std::cout<<"Global number of points found not on boundary is "<<globalPtsFoundNotOnBoundary<<std::endl;
    std::cout<<"Global number of points found on boundary is "<<globalPtsFoundOnBoundary<<std::endl;
    std::cout<<"Total number of points is "<<globalPtsNotFound+globalPtsFoundNotOnBoundary+globalPtsFoundOnBoundary<<std::endl;
  }

  ut->passes(exeName);
}


int main(int argc, char *argv[])
{
  AMP::AMPManager::startup(argc, argv);
  AMP::UnitTest ut;

  std::string exeName = "testDendroInterpolation";

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



