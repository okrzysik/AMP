
#include "utils/AMPManager.h"
#include "utils/UnitTest.h"
#include "utils/Utilities.h"
#include "utils/Database.h"
#include "utils/InputDatabase.h"
#include "utils/InputManager.h"
#include "utils/AMP_MPI.h"

#include "discretization/DOF_Manager.h"
#include "discretization/simpleDOF_Manager.h"
#include "vectors/Variable.h"
#include "vectors/Vector.h"
#include "vectors/VectorBuilder.h"

#include "externVars.h"

#include "ampmesh/Mesh.h"
#include "ampmesh/dendro/DendroSearch.h"

#include <fstream>
#include <iomanip>
#include <numeric>
#include <boost/lexical_cast.hpp>

double dummyFunction(const std::vector<double> &xyz, const int dof) {
  AMP_ASSERT(xyz.size() == 3);
  double x = xyz[0], y = xyz[1], z = xyz[2];
  //  return 7.0;
  //return (1.0 + 6.0 * x) * (2.0 - 5.0 * y) * (3.0 + 4.0 * z);
  return (1.0 + 6.0 * x) + (2.0 - 5.0 * y) + (3.0 + 4.0 * z);
}

void genGaussPts(int rank, size_t numOctPtsPerProc, std::vector<double> & pts);
void genUniformPts(int rank, size_t numOctPtsPerProc, std::vector<double> & pts);

void rescalePts(std::vector<double> & pts);

double gaussian(double mean, double std_deviation);

void run(const std::string & meshFileName, 
         size_t numRandomPts, 
         void (*randomPtsGenerator)(int, size_t, std::vector<double>&),
         std::vector<double> & timingMeasurements) {

  AMP::AMP_MPI globalComm(AMP_COMM_WORLD);

  int rank = globalComm.getRank();
  int npes = globalComm.getSize();

  // Load the mesh
  globalComm.barrier();
  double meshBeginTime = MPI_Wtime();
  boost::shared_ptr<AMP::InputDatabase> mesh_db(new AMP::InputDatabase("input_db"));
  mesh_db->putString("MeshName", "PelletMesh");
  mesh_db->putString("MeshType", "libMesh");
  mesh_db->putString("FileName", meshFileName);
  mesh_db->putInteger("dim", 3);
  boost::shared_ptr<AMP::Mesh::MeshParameters> meshParams(new AMP::Mesh::MeshParameters(mesh_db));
  meshParams->setComm(AMP::AMP_MPI(AMP_COMM_WORLD));
  AMP::Mesh::Mesh::shared_ptr meshAdapter = AMP::Mesh::Mesh::buildMesh(meshParams);
  globalComm.barrier();
  double meshEndTime = MPI_Wtime();
  if(!rank) {
    std::cout<<"Finished reading the mesh in "<<(meshEndTime - meshBeginTime)<<" seconds."<<std::endl;
  }

  // Build the vector field
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

  // Generate the random points
  int totalNumPts = numRandomPts;

  int avgNumPts = totalNumPts/npes;
  int extraNumPts = totalNumPts%npes;

  int numLocalPts = avgNumPts;
  if(rank < extraNumPts) {
    numLocalPts++;
  }

  std::vector<double> pts;
  randomPtsGenerator(rank, numLocalPts, pts);

  double minCoords[3];
  double maxCoords[3];
  std::vector<double> box = meshAdapter->getBoundingBox();
  for(unsigned char i = 0; i < meshAdapter->getDim(); ++i) {
    minCoords[i] = box[2*i+0];
    maxCoords[i] = box[2*i+1];
  } // end for i

  for(size_t i = 0; i < numLocalPts; ++i) {
    pts[3*i+0] = (maxCoords[0] - minCoords[0])*pts[3*i+0] + minCoords[0];
    pts[3*i+1] = (maxCoords[1] - minCoords[1])*pts[3*i+1] + minCoords[1];
    pts[3*i+2] = (maxCoords[2] - minCoords[2])*pts[3*i+2] + minCoords[2];
  } // end for i
  if(!rank) {
    std::cout<<"Finished generating "<<totalNumPts <<" random points for search!"<<std::endl;
  }

  // Perform the search
  DendroSearch dendroSearch(meshAdapter, true);
  dendroSearch.search(globalComm, pts);

  // Interpolate
  std::vector<double> interpolatedData; 
  std::vector<bool> interpolationWasDone;
  dendroSearch.interpolate(globalComm, dummyVector, DOFsPerNode, interpolatedData, interpolationWasDone);

  std::vector<double> interpolationError((DOFsPerNode*numLocalPts), 0.0);
  for(unsigned int i = 0; i < numLocalPts; ++i) {
    if(interpolationWasDone[i]) {
      for(int d = 0; d < DOFsPerNode; ++d) {
        interpolationError[(i*DOFsPerNode) + d] = fabs(interpolatedData[(i*DOFsPerNode) + d] -
            dummyFunction(std::vector<double>(&(pts[3*i]), &(pts[3*i]) + 3), d));
      }//end d
    }
  } // end for i
  double localMaxError = 0.0;
  if(numLocalPts > 0) {
    localMaxError = *(std::max_element(interpolationError.begin(), interpolationError.end()));
  }
  double globalMaxError = globalComm.maxReduce<double>(localMaxError);
  if(!rank) {
    std::cout<<"Intepolation global max error is "<<std::setprecision(15)<<globalMaxError<<std::endl;
  }

  AMP_ASSERT(globalMaxError < 1.0e-12);

  // Project on boundary
  std::vector<AMP::Mesh::MeshElementID> faceVerticesGlobalIDs;
  std::vector<double> shiftGlobalCoords, projectionLocalCoordsOnFace;
  std::vector<int> flags;
  dendroSearch.projectOnBoundaryID(globalComm, 4, faceVerticesGlobalIDs, shiftGlobalCoords, projectionLocalCoordsOnFace, flags);

  size_t localPtsNotFound = std::count(flags.begin(), flags.end(), DendroSearch::NotFound);
  size_t localPtsFoundNotOnBoundary = std::count(flags.begin(), flags.end(), DendroSearch::FoundNotOnBoundary);
  size_t localPtsFoundOnBoundary = std::count(flags.begin(), flags.end(), DendroSearch::FoundOnBoundary);
  size_t globalPtsNotFound = globalComm.sumReduce(localPtsNotFound);
  size_t globalPtsFoundNotOnBoundary = globalComm.sumReduce(localPtsFoundNotOnBoundary);
  size_t globalPtsFoundOnBoundary = globalComm.sumReduce(localPtsFoundOnBoundary);
  if(!rank) {
    std::cout<<"Global number of points not found is "<<globalPtsNotFound<<std::endl;
    std::cout<<"Global number of points found not on boundary is "<<globalPtsFoundNotOnBoundary<<std::endl;
    std::cout<<"Global number of points found on boundary is "<<globalPtsFoundOnBoundary<<std::endl;
    std::cout<<"Total number of points is "<<globalPtsNotFound+globalPtsFoundNotOnBoundary+globalPtsFoundOnBoundary<<std::endl;
  }

  DendroSearch::TimingType timingTypes[5] = { 
    DendroSearch::Setup, 
    DendroSearch::CoarseSearch, 
    DendroSearch::FineSearch,
    DendroSearch::Interpolation,
    DendroSearch::ProjectionOnBoundaryID
  };
  timingMeasurements.resize(5);
  dendroSearch.reportTiming(5, timingTypes, &(timingMeasurements[0]));

  timingMeasurements.push_back(timingMeasurements[1]+timingMeasurements[2]);
  timingMeasurements.push_back(timingMeasurements[0]+timingMeasurements[1]+timingMeasurements[2]);
//  localTimingMeasurements[4] = std::accumulate(localTimingMeasurements, localTimingMeasurements+4, 0.0);
//  
//  globalComm.maxReduce(localTimingMeasurements, globalTimingMeasurements, 5);
//
//  if(!rank) {
//    std::cout<<"Setup time = "<<globalTimingMeasurements[0]<<"\n";
//    std::cout<<"Coarse Search time = "<<globalTimingMeasurements[1]<<"\n";
//    std::cout<<"Fine Search time = "<<globalTimingMeasurements[2]<<"\n";
//    std::cout<<"Interpolation time = "<<globalTimingMeasurements[3]<<"\n";
//    std::cout<<"Total time = "<<globalTimingMeasurements[4]<<"\n";
//  }
}

void myTest(AMP::UnitTest *ut, std::string exeName) {
  AMP::AMP_MPI globalComm(AMP_COMM_WORLD);

  int rank = globalComm.getRank();
  int npes = globalComm.getSize();

  size_t n_k = 3;
  std::string meshFileNames[] = { "pellet_4x.e", "pellet_2x.e", "pellet_4x.e" };
  size_t meshNumElements[] = { 3705, 26146, 183210 };

  size_t n_j = 1;
  size_t numRandomPts[] = { 10000, 20000, 40000, 80000, 160000 };
  
  size_t n_i = 1;
  void (*randomPtsGenerators[])(int, size_t, std::vector<double>&) = { &genUniformPts, &genGaussPts };
  std::string prefixes[] = { "uniform", "gaussian" };

  std::string suffix = boost::lexical_cast<std::string>(npes);


  for (size_t i = 0; i < n_i; ++i) {
    std::fstream fout;
    if (!rank) { 
      std::string fileName = prefixes[i] + "_" + suffix;
      fout.open(fileName.c_str(), std::fstream::out);
      fout<<n_j<<"  "<<n_k<<"\n";
      for (size_t j = 0; j < n_j; ++j) {
        fout<<numRandomPts[j]<<"  ";
      } // end for j
      fout<<"\n";
      for (size_t k = 0; k < n_k; ++k) { 
        fout<<meshNumElements[k]<<"  ";
      } // end for k
      fout<<"\n";
    } // end if
    for (size_t j = 0; j < n_j; ++j) {
      for (size_t k = 0; k < n_k; ++k) {
        std::vector<double> localTimingMeasurements;
        run(meshFileNames[k], numRandomPts[j], randomPtsGenerators[i], localTimingMeasurements);
        size_t numTimingMeasurements = localTimingMeasurements.size();
        AMP_ASSERT(numTimingMeasurements > 0);
        std::vector<double> globalTimingMeasurements(numTimingMeasurements, -1.0); 
        globalComm.maxReduce(&(localTimingMeasurements[0]), &(globalTimingMeasurements[0]), numTimingMeasurements);
        if (!rank) {
          fout<<globalTimingMeasurements[0]<<"  " // setup
              <<globalTimingMeasurements[1]<<"  " // coarse
              <<globalTimingMeasurements[2]<<"  " // fine
              <<globalTimingMeasurements[5]<<"  "<<std::flush; // coarse+fine
        } // end if
      } // end for k
      if (!rank) {
        fout<<"\n";
      } // end if
    } // end for j
    if (!rank) {
      fout.close();
    } // end if
  } // end for i

  ut->passes(exeName);
}


int main(int argc, char *argv[])
{
  AMP::AMPManager::startup(argc, argv);
  AMP::UnitTest ut;

  std::string exeName = "testDendroSearch";

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

void genUniformPts(int rank, size_t numPtsPerProc, std::vector<double> & pts) {
  const int seed = (0x1234567 + (24135*rank));
  srand48(seed);

  pts.resize(3*numPtsPerProc);
  for(size_t i = 0; i < 3*numPtsPerProc; ++i) {
    pts[i] = drand48();
  } // end for i
}

void genGaussPts(int rank, size_t numPtsPerProc, std::vector<double> & pts)
{
  const int seed = (0x12345678  + (76543*rank));
  srand48(seed);

  pts.resize(3*numPtsPerProc);
  for(size_t i = 0; i < (3*numPtsPerProc); i++) {
    pts[i]= gaussian(0.5, 0.16);
  } // end for i
}

double gaussian(double mean, double std_deviation) {
  static double t1 = 0, t2=0;
  double x1, x2, x3, r;

  using namespace std;

  // reuse previous calculations
  if(t1) {
    const double tmp = t1;
    t1 = 0;
    return mean + std_deviation * tmp;
  }
  if(t2) {
    const double tmp = t2;
    t2 = 0;
    return mean + std_deviation * tmp;
  }

  // pick randomly a point inside the unit disk
  do {
    x1 = 2 * drand48() - 1;
    x2 = 2 * drand48() - 1;
    x3 = 2 * drand48() - 1;
    r = x1 * x1 + x2 * x2 + x3*x3;
  } while(r >= 1);

  // Box-Muller transform
  r = sqrt(-2.0 * log(r) / r);

  // save for next call
  t1 = (r * x2);
  t2 = (r * x3);

  return mean + (std_deviation * r * x1);
}//end gaussian

void rescalePts(std::vector<double> & pts) 
{
  double minX = pts[0];
  double maxX = pts[0];

  double minY = pts[1];
  double maxY = pts[1];

  double minZ = pts[2];
  double maxZ = pts[2];

  for(unsigned int i = 0; i < pts.size(); i += 3) {
    double xPt = pts[i];
    double yPt = pts[i + 1];
    double zPt = pts[i + 2];

    if(xPt < minX) {
      minX = xPt;
    }

    if(xPt > maxX) {
      maxX = xPt;
    }

    if(yPt < minY) {
      minY = yPt;
    }

    if(yPt > maxY) {
      maxY = yPt;
    }

    if(zPt < minZ) {
      minZ = zPt;
    }

    if(zPt > maxZ) {
      maxZ = zPt;
    }
  }//end for i

  double xRange = (maxX - minX);
  double yRange = (maxY - minY);
  double zRange = (maxZ - minZ);

  for(unsigned int i = 0; i < pts.size();  i += 3) {
    double xPt = pts[i];
    double yPt = pts[i + 1];
    double zPt = pts[i + 2];

    pts[i] = 0.05 + (0.925*(xPt - minX)/xRange);
    pts[i + 1] = 0.05 + (0.925*(yPt - minY)/yRange);
    pts[i + 2] = 0.05 + (0.925*(zPt - minZ)/zRange);
  }//end for i
}

