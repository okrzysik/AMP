
#include "mpi.h"

#include "utils/AMPManager.h"
#include "utils/UnitTest.h"
#include "utils/Utilities.h"
#include "utils/Database.h"
#include "utils/InputDatabase.h"
#include "utils/InputManager.h"
#include "utils/AMP_MPI.h"
#include "utils/PIO.h"

#include "ampmesh/Mesh.h"

#include <iostream>
#include <string>
#include <cmath>
#include <cstdlib>

#include <vector>

#include "sys/sys.h"
#include "par/parUtils.h"
#include "binOps/binUtils.h"
#include "oct/TreeNode.h"
#include "oct/octUtils.h"
#include "oct/nodeAndValues.h"
#include "oct/nodeAndRanks.h"
#include "externVars.h"
#include "dendro.h"

void computeMinAndMaxCoords( double* minCoords, double* maxCoords, double* ScalingFactor, 
    AMP::Mesh::Mesh::shared_ptr meshAdapter, AMP::AMP_MPI globalComm ) {
  AMP::Mesh::MeshIterator nd = meshAdapter->getIterator(AMP::Mesh::Vertex, 0);
  AMP::Mesh::MeshIterator end_nd = nd.end();

  AMP_ASSERT(nd != end_nd);
  {
    std::vector<double> pt = nd->coord();
    for(int i = 0; i < pt.size(); ++i) {
      minCoords[i] = pt[i];
      maxCoords[i] = pt[i];
    }//end i
    ++nd;
  }
  for( ; nd != end_nd; ++nd) {
    std::vector<double> pt = nd->coord();
    for(int i = 0; i < pt.size(); ++i) {
      if(minCoords[i] > pt[i]) {
        minCoords[i] = pt[i];
      }
      if(maxCoords[i] < pt[i]) {
        maxCoords[i] = pt[i];
      }
    }//end i
  }//end nd

  globalComm.minReduce<double>(minCoords, 3, NULL);
  globalComm.maxReduce<double>(maxCoords, 3, NULL);

  for(int i = 0; i < 3; ++i) {
    ScalingFactor[i] = 1.0/(1.0e-10 + maxCoords[i] - minCoords[i]);
  }//end i
}

void createLocalMeshElementArray(std::vector<AMP::Mesh::MeshElement>& localElemArr, 
    AMP::Mesh::Mesh::shared_ptr meshAdapter) {
  localElemArr.clear();
  AMP::Mesh::MeshIterator el = meshAdapter->getIterator(AMP::Mesh::Volume, 0);
  AMP::Mesh::MeshIterator end_el = el.end();
  AMP_ASSERT(el != end_el);
  for(; el != end_el; ++el) {
    localElemArr.push_back(*el);
  }//end el
}

void setupDSforSearchType(unsigned int & BoxLevel, std::vector<ot::TreeNode>& nodeList, std::vector<unsigned int>& numIndicesList,
    std::vector<ot::TreeNode>& mins, std::vector<unsigned int>& rankList,
    std::vector<int>& elemIdList, std::vector<AMP::Mesh::MeshElement>& localElemArr,
    double* minCoords, double* maxCoords, double* ScalingFactor, 
    AMP::Mesh::Mesh::shared_ptr meshAdapter, AMP::AMP_MPI globalComm) {
  int rank = globalComm.getRank();
  int npes = globalComm.getSize();

  computeMinAndMaxCoords(minCoords, maxCoords, ScalingFactor, meshAdapter, globalComm);

  createLocalMeshElementArray(localElemArr, meshAdapter);

  const unsigned int MaxDepth = 30;

  unsigned int totalNumElems = meshAdapter->numGlobalElements(AMP::Mesh::Volume);

  double avgHboxInv = std::pow(totalNumElems, (1.0/3.0));
  assert(avgHboxInv > 1.0);
  BoxLevel = binOp::fastLog2(static_cast<unsigned int>(std::ceil(avgHboxInv)));
  assert(BoxLevel < MaxDepth);

  if(!rank) {
    std::cout<<"BoxLevel = "<<BoxLevel<<std::endl;
  }

  const unsigned int ITPMD = (1u << MaxDepth);
  const double DTPMD = static_cast<double>(ITPMD);
  const double hBox = 1.0/(static_cast<double>(1u<<BoxLevel));

  std::vector< ot::NodeAndValues<int, 1> > nodeAndElemIdList;

  assert(!(localElemArr.empty()));

  for(int eId = 0; eId < localElemArr.size(); ++eId) {
    std::vector<AMP::Mesh::MeshElement> currNodes = localElemArr[eId].getElements(AMP::Mesh::Vertex);
    int minId[3];
    int maxId[3];
    for(size_t i = 0; i < currNodes.size(); ++i) {
      std::vector<double> pt = currNodes[i].coord();
      double scaledPt[3];
      for(int j = 0; j < 3; ++j) {
        scaledPt[j] = ((pt[j] - minCoords[j])*ScalingFactor[j]);
        int id = static_cast<int>(scaledPt[j]/hBox);
        if(i == 0) {
          minId[j] = id;
          maxId[j] = id;
        } else {
          if(minId[j] > id) {
            minId[j] = id;
          }
          if(maxId[j] < id) {
            maxId[j] = id;
          }
        }
      }//end j
    }//end i
    //Performance Improvement: We can skip the boxes that lie
    //completely outside the element.
    for(int k = minId[2]; k <= maxId[2]; ++k) {
      for(int j = minId[1]; j <= maxId[1]; ++j) {
        for(int i = minId[0]; i <= maxId[0]; ++i) {
          unsigned int bX = i*(1u<<(MaxDepth - BoxLevel));
          unsigned int bY = j*(1u<<(MaxDepth - BoxLevel));
          unsigned int bZ = k*(1u<<(MaxDepth - BoxLevel));
          ot::TreeNode box(bX, bY, bZ, BoxLevel, 3, MaxDepth);
          box.setWeight(rank);
          ot::NodeAndValues<int, 1> obj;
          obj.node = box;
          obj.values[0] = eId;
          nodeAndElemIdList.push_back(obj);
        }//end i 
      }//end j 
    }//end k 
  }//end eId

  std::vector< ot::NodeAndValues<int, 1> > tmpList;
  par::sampleSort< ot::NodeAndValues<int, 1> >(
      nodeAndElemIdList, tmpList, (globalComm.getCommunicator()));
  swap(nodeAndElemIdList, tmpList);
  tmpList.clear();

  int numLocalOcts = nodeAndElemIdList.size();
  int numGlobalOcts = globalComm.sumReduce<int>(numLocalOcts);
  if(!rank) {
    std::cout<<"Total num initial octants = "<<numGlobalOcts <<std::endl;
  }

  assert(rankList.empty());
  assert(elemIdList.empty());

  for(int i = 0; i < numLocalOcts; ++i) {
    ot::TreeNode currNode = nodeAndElemIdList[i].node;
    rankList.push_back(currNode.getWeight());
    elemIdList.push_back(nodeAndElemIdList[i].values[0]);
  }//end i

  assert(numLocalOcts > 0);
  assert(nodeList.empty());
  assert(numIndicesList.empty());

  //Local Merge
  {
    ot::TreeNode currNode = nodeAndElemIdList[0].node;
    nodeList.push_back(currNode);
    numIndicesList.push_back(1);
  }
  for(size_t i = 1; i < nodeAndElemIdList.size(); ++i) {
    ot::TreeNode currNode = nodeAndElemIdList[i].node;
    if( nodeList[nodeList.size() - 1] == currNode ) {
      numIndicesList[nodeList.size() - 1]++;
    } else {
      nodeList.push_back(currNode);
      numIndicesList.push_back(1);
    }
  }//end i
  nodeAndElemIdList.clear();

  assert(nodeList.size() >= 2);

  //PARALLEL MERGE HERE

  if(!(nodeList.empty())) {
    nodeList[0].setWeight(0);
  }
  for(int i = 1; i < nodeList.size(); ++i) {
    nodeList[i].setWeight(nodeList[i - 1].getWeight() + numIndicesList[i - 1]);
  }//end i

  ot::TreeNode firstNode;
  if(!(nodeList.empty())) {
    firstNode = nodeList[0];
    firstNode.setWeight(rank);
  }
  mins.resize(npes);
  MPI_Allgather(&firstNode, 1, par::Mpi_datatype<ot::TreeNode>::value(), 
      &(mins[0]), 1, par::Mpi_datatype<ot::TreeNode>::value(), globalComm.getCommunicator() );

  std::vector<ot::TreeNode> tmpMins;
  for(int i = 0; i < npes; ++i) {
    if(mins[i].getDim() > 0) {
      tmpMins.push_back(mins[i]);
    }
  }//end i
  swap(mins, tmpMins);
  tmpMins.clear();

  int minNumIndices = numIndicesList[0];
  int maxNumIndices = numIndicesList[0];
  for(int i = 1; i < numIndicesList.size(); ++i) {
    if(minNumIndices > numIndicesList[i]) {
      minNumIndices = numIndicesList[i];
    }
    if(maxNumIndices < numIndicesList[i]) {
      maxNumIndices = numIndicesList[i];
    }
  }//end i

  int globalMinNumIndices = globalComm.minReduce<int>(minNumIndices);
  int globalMaxNumIndices = globalComm.maxReduce<int>(maxNumIndices);

  numLocalOcts = nodeList.size();
  numGlobalOcts = globalComm.sumReduce<int>(numLocalOcts);

  if(!rank) {
    std::cout<<"Total num final octants = "<<numGlobalOcts <<std::endl;
    std::cout<<"Global Min Num Indices = "<<globalMinNumIndices <<std::endl;
    std::cout<<"Global Max Num Indices = "<<globalMaxNumIndices <<std::endl;
  }
}

void myTest(AMP::UnitTest *ut, std::string exeName) {
  std::string input_file = "input_" + exeName;
  std::string log_file = "output_" + exeName; 

  ot::RegisterEvents();

  AMP::PIO::logOnlyNodeZero(log_file);
  AMP::AMP_MPI globalComm(AMP_COMM_WORLD);

  int rank = globalComm.getRank();
  int npes = globalComm.getSize();

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

  globalComm.barrier();

  double setupBeginTime = MPI_Wtime();

  double minCoords[3];
  double maxCoords[3];
  double ScalingFactor[3];

  std::vector<ot::TreeNode> nodeList;
  std::vector<unsigned int> numIndicesList;
  std::vector<unsigned int> rankList;
  std::vector<int> elemIdList;
  std::vector<AMP::Mesh::MeshElement> localElemArr;
  std::vector<ot::TreeNode> mins;
  unsigned int BoxLevel;

  setupDSforSearchType(BoxLevel, nodeList, numIndicesList, mins, rankList, elemIdList, localElemArr,
      minCoords, maxCoords, ScalingFactor, meshAdapter, globalComm );

  globalComm.barrier();

  double setupEndTime = MPI_Wtime();

  if(!rank) {
    std::cout<<"Finished setting up DS for search in "<<(setupEndTime - setupBeginTime)<<" seconds."<<std::endl;
  }

  int totalNumPts = input_db->getInteger("TotalNumberOfPoints");
  int avgNumPts = totalNumPts/npes;
  int extraNumPts = totalNumPts%npes;

  int numLocalPts = avgNumPts;
  if(rank < extraNumPts) {
    numLocalPts++;
  }

  //Generate Random points in [min, max]
  const unsigned int seed = (0x1234567 + (24135*rank));
  srand48(seed);

  std::vector<double> pts;
  for(int i = 0; i < numLocalPts; ++i) {
    double x = ((maxCoords[0] - minCoords[0])*drand48()) + minCoords[0];
    double y = ((maxCoords[1] - minCoords[1])*drand48()) + minCoords[1];
    double z = ((maxCoords[2] - minCoords[2])*drand48()) + minCoords[2];
    pts.push_back(x);
    pts.push_back(y);
    pts.push_back(z);
  }//end i

  if(!rank) {
    std::cout<<"Finished generating random points for search!"<<std::endl;
  }

  globalComm.barrier();

  double searchBeginTime = MPI_Wtime();

  const unsigned int MaxDepth = 30;
  const unsigned int ITPMD = (1u << MaxDepth);
  const double DTPMD = static_cast<double>(ITPMD);

  std::vector<ot::NodeAndValues<double, 4> > ptsWrapper;
  for(int i = 0; i < numLocalPts; ++i) {
    double x = pts[3*i];
    double y = pts[(3*i) + 1];
    double z = pts[(3*i) + 2];
    double scaledX = ((x - minCoords[0])*ScalingFactor[0]);
    double scaledY = ((y - minCoords[1])*ScalingFactor[1]);
    double scaledZ = ((z - minCoords[2])*ScalingFactor[2]);
    unsigned int pX = static_cast<unsigned int>(scaledX*DTPMD);
    unsigned int pY = static_cast<unsigned int>(scaledY*DTPMD);
    unsigned int pZ = static_cast<unsigned int>(scaledZ*DTPMD);

    ot::NodeAndValues<double, 4> tmpObj;
    tmpObj.node =  ot::TreeNode(pX, pY, pZ, MaxDepth, 3, MaxDepth);
    tmpObj.node.setWeight(rank);
    tmpObj.values[0] = x;
    tmpObj.values[1] = y;
    tmpObj.values[2] = z;
    tmpObj.values[3] = i;

    ptsWrapper.push_back(tmpObj);
  }//end i

  //Performance Question: Should PtsWrapper be sorted or not?
  //If PtsWrapper is sorted (even just a local sort), we can skip the
  //binary searches and use binning instead.  Binning is amortized constant
  //time and using binary searches would be logarithmic. This is just a matter
  //of constants since sorting is also logarithmic.

  int* sendCnts = new int[npes];
  for(int i = 0; i < npes; ++i) {
    sendCnts[i] = 0;
  }//end i

  std::vector<int> part(numLocalPts, -1);
  for(int i = 0; i < numLocalPts; ++i) {
    unsigned int retIdx;
    bool found = seq::maxLowerBound<ot::TreeNode>(mins, (ptsWrapper[i].node), retIdx, NULL, NULL);
    if(found) {
      part[i] = mins[retIdx].getWeight();
      sendCnts[part[i]]++;
    }
  }//end i

  int* recvCnts = new int[npes];
  MPI_Alltoall(sendCnts, 1, MPI_INT, recvCnts, 1, MPI_INT, (globalComm.getCommunicator()));

  int* sendDisps = new int[npes]; 
  int* recvDisps = new int[npes]; 
  sendDisps[0] = 0;
  recvDisps[0] = 0;
  for(int i = 1; i < npes; ++i) {
    sendDisps[i] = sendDisps[i - 1] + sendCnts[i - 1];
    recvDisps[i] = recvDisps[i - 1] + recvCnts[i - 1];
  }//end i

  std::vector<ot::NodeAndValues<double, 4> > sendList(sendDisps[npes - 1] + sendCnts[npes - 1]);

  for(int i = 0; i < npes; ++i) {
    sendCnts[i] = 0;
  }//end i

  for(int i = 0; i < numLocalPts; ++i) {
    if(part[i] >= 0) {
      sendList[sendDisps[part[i]] + sendCnts[part[i]]] = ptsWrapper[i];
      sendCnts[part[i]]++;
    }
  }//end i

  std::vector<ot::NodeAndValues<double, 4> > recvList(recvDisps[npes - 1] + recvCnts[npes - 1]);
  ot::NodeAndValues<double, 4>* sendListPtr = NULL;
  ot::NodeAndValues<double, 4>* recvListPtr = NULL;
  if(!(sendList.empty())) {
    sendListPtr = &(sendList[0]);
  }
  if(!(recvList.empty())) {
    recvListPtr = &(recvList[0]);
  }
  MPI_Alltoallv( sendListPtr, sendCnts, sendDisps, par::Mpi_datatype<ot::NodeAndValues<double, 4> >::value(),
      recvListPtr, recvCnts, recvDisps, par::Mpi_datatype<ot::NodeAndValues<double, 4> >::value(), 
      (globalComm.getCommunicator()) );

  for(int i = 0; i < npes; ++i) {
    sendCnts[i] = 0;
  }//end i

  std::vector<int> ptToOctMap((recvList.size()), -1);
  for(int i = 0; i < recvList.size(); ++i) {
    unsigned int retIdx;
    seq::maxLowerBound<ot::TreeNode>(nodeList, (recvList[i].node), retIdx, NULL, NULL);
    if( nodeList[retIdx].isAncestor(recvList[i].node) ) {
      ptToOctMap[i] = retIdx;
      int stIdx = nodeList[retIdx].getWeight();
      for(int j = 0; j < numIndicesList[retIdx]; ++j) {
        sendCnts[rankList[stIdx + j]]++;
      }//end j
    }
  }//end i

  MPI_Alltoall(sendCnts, 1, MPI_INT, recvCnts, 1, MPI_INT, (globalComm.getCommunicator()));

  sendDisps[0] = 0;
  recvDisps[0] = 0;
  for(int i = 1; i < npes; ++i) {
    sendDisps[i] = sendDisps[i - 1] + sendCnts[i - 1];
    recvDisps[i] = recvDisps[i - 1] + recvCnts[i - 1];
  }//end i

  std::vector<int> sendElemIdList(sendDisps[npes - 1] + sendCnts[npes - 1]);
  std::vector<double> sendPtsList(5*(sendDisps[npes - 1] + sendCnts[npes - 1]));

  for(int i = 0; i < npes; ++i) {
    sendCnts[i] = 0;
  }//end i

  for(int i = 0; i < ptToOctMap.size(); ++i) {
    if(ptToOctMap[i] >= 0) {
      int stIdx = nodeList[ptToOctMap[i]].getWeight();
      for(int j = 0; j < numIndicesList[ptToOctMap[i]]; ++j) {
        sendElemIdList[sendDisps[rankList[stIdx + j]] + sendCnts[rankList[stIdx + j]]] = elemIdList[stIdx + j];
        sendPtsList[(5*(sendDisps[rankList[stIdx + j]] + sendCnts[rankList[stIdx + j]]))] = recvList[i].values[0];
        sendPtsList[(5*(sendDisps[rankList[stIdx + j]] + sendCnts[rankList[stIdx + j]])) + 1] = recvList[i].values[1];
        sendPtsList[(5*(sendDisps[rankList[stIdx + j]] + sendCnts[rankList[stIdx + j]])) + 2] = recvList[i].values[2];
        sendPtsList[(5*(sendDisps[rankList[stIdx + j]] + sendCnts[rankList[stIdx + j]])) + 3] = recvList[i].values[3];
        sendPtsList[(5*(sendDisps[rankList[stIdx + j]] + sendCnts[rankList[stIdx + j]])) + 4] = recvList[i].node.getWeight();
        sendCnts[rankList[stIdx + j]]++;
      }//end j
    }
  }//end i

  std::vector<int> recvElemIdList(recvDisps[npes - 1] + recvCnts[npes - 1]);
  int* sendElemPtr = NULL;
  int* recvElemPtr = NULL;
  if(!(sendElemIdList.empty())) {
    sendElemPtr = &(sendElemIdList[0]);
  }
  if(!(recvElemIdList.empty())) {
    recvElemPtr = &(recvElemIdList[0]);
  }
  MPI_Alltoallv( sendElemPtr, sendCnts, sendDisps, MPI_INT,
      recvElemPtr, recvCnts, recvDisps, MPI_INT, (globalComm.getCommunicator()) );

  for(int i = 0; i < npes; ++i) {
    sendCnts[i] *= 5;
    sendDisps[i] *= 5;
    recvCnts[i] *= 5;
    recvDisps[i] *= 5;
  }//end i

  std::vector<double> recvPtsList(recvDisps[npes - 1] + recvCnts[npes - 1]);
  double* sendPtsPtr = NULL;
  double* recvPtsPtr = NULL;
  if(!(sendPtsList.empty())) {
    sendPtsPtr = &(sendPtsList[0]);
  }
  if(!(recvPtsList.empty())) {
    recvPtsPtr = &(recvPtsList[0]);
  }
  MPI_Alltoallv( sendPtsPtr, sendCnts, sendDisps, MPI_DOUBLE,
      recvPtsPtr, recvCnts, recvDisps, MPI_DOUBLE, (globalComm.getCommunicator()) );

  delete [] sendCnts;
  delete [] sendDisps;
  delete [] recvCnts;
  delete [] recvDisps;

  int numRecvPts = recvElemIdList.size();
  std::vector<bool> results(numRecvPts, false);
  for(int i = 0; i < numRecvPts; ++i) {
    assert(recvElemIdList[i] >= 0);
    assert(recvElemIdList[i] < localElemArr.size());
    AMP::Mesh::MeshElement el = localElemArr[recvElemIdList[i]];
    //results[i] = el.containsPoint(recvPtsList[5*i], recvPtsList[(5*i) + 1], recvPtsList[(5*i) + 2]);
  }//end i

  globalComm.barrier();

  double searchEndTime = MPI_Wtime();

  if(!rank) {
    std::cout<<"Finished searching "<<totalNumPts<<" points in "<<
      (searchEndTime - searchBeginTime)<<" seconds using "<<npes<<"processors."<<std::endl;
  }

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



