
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

#include "mpi.h"
#include <vector>

#include "sys/sys.h"
#include "par/parUtils.h"
#include "oct/TreeNode.h"
#include "oct/octUtils.h"
#include "oct/nodeAndValues.h"
#include "oct/nodeAndRanks.h"
#include "externVars.h"
#include "dendro.h"

void myTest(AMP::UnitTest *ut, std::string exeName) {
  std::string input_file = "input_" + exeName;
  std::string log_file = "output_" + exeName; 

  ot::RegisterEvents();

  AMP::PIO::logOnlyNodeZero(log_file);
  AMP::AMP_MPI globalComm(AMP_COMM_WORLD);

  boost::shared_ptr<AMP::InputDatabase> input_db(new AMP::InputDatabase("input_db"));
  AMP::InputManager::getManager()->parseInputFile(input_file, input_db);
  input_db->printClassData(AMP::plog);

  AMP_INSIST(input_db->keyExists("Mesh"), "Key ''Mesh'' is missing!");
  boost::shared_ptr<AMP::Database> mesh_db = input_db->getDatabase("Mesh");
  boost::shared_ptr<AMP::Mesh::MeshParameters> meshParams(new AMP::Mesh::MeshParameters(mesh_db));
  meshParams->setComm(AMP::AMP_MPI(AMP_COMM_WORLD));
  AMP::Mesh::Mesh::shared_ptr meshAdapter = AMP::Mesh::Mesh::buildMesh(meshParams);

  int rank = globalComm.getRank();

  AMP::Mesh::MeshIterator nd = meshAdapter->getIterator(AMP::Mesh::Vertex, 0);
  AMP::Mesh::MeshIterator end_nd = nd.end();

  double minCoords[3];
  double maxCoords[3];

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

  double ScalingFactor[3];
  for(int i = 0; i < 3; ++i) {
    ScalingFactor[i] = 1.0/(1.0e-10 + maxCoords[i] - minCoords[i]);
  }//end i

  const unsigned int MaxDepth = 30;
  const unsigned int ITPMD = (1u << MaxDepth);
  const double DTPMD = static_cast<double>(ITPMD);

  std::vector< ot::NodeAndValues<AMP::Mesh::MeshElementID, 1> > nodeAndElemIdList;

  AMP::Mesh::MeshIterator el = meshAdapter->getIterator(AMP::Mesh::Volume, 0);
  AMP::Mesh::MeshIterator end_el = el.end();
  AMP_ASSERT(el != end_el);
  for( ; el != end_el; ++el) {
    std::vector<AMP::Mesh::MeshElement> currNodes = el->getElements(AMP::Mesh::Vertex);
    std::vector<ot::TreeNode> ptOcts;
    for(size_t i = 0; i < currNodes.size(); ++i) {
      std::vector<double> pt = currNodes[i].coord();
      double scaledX = ((pt[0] - minCoords[0])*ScalingFactor[0]);
      double scaledY = ((pt[1] - minCoords[1])*ScalingFactor[1]);
      double scaledZ = ((pt[2] - minCoords[2])*ScalingFactor[2]);
      unsigned int pX = static_cast<unsigned int>(scaledX*DTPMD);
      unsigned int pY = static_cast<unsigned int>(scaledY*DTPMD);
      unsigned int pZ = static_cast<unsigned int>(scaledZ*DTPMD);
      ptOcts.push_back( ot::TreeNode(pX, pY, pZ, MaxDepth, 3, MaxDepth) );
    }//end i
    ot::TreeNode nca = ptOcts[0];
    for(size_t i = 1; i < ptOcts.size(); ++i) {
      nca = ot::getNCA(nca, ptOcts[i]);
    }//end i
    nca.setWeight(rank);
    ot::NodeAndValues<AMP::Mesh::MeshElementID, 1> obj;
    obj.node = nca;
    obj.values[0] = el->globalID();
    nodeAndElemIdList.push_back(obj);
  }//end for el

  std::vector< ot::NodeAndValues<AMP::Mesh::MeshElementID, 1> > tmpList;
  par::sampleSort< ot::NodeAndValues<AMP::Mesh::MeshElementID, 1> >(
      nodeAndElemIdList, tmpList, (globalComm.getCommunicator()));
  swap(nodeAndElemIdList, tmpList);
  tmpList.clear();

  AMP_ASSERT(!(nodeAndElemIdList.empty()));

  std::vector<unsigned int> rankList;
  std::vector<AMP::Mesh::MeshElementID> elemIdList;
  std::vector<ot::NodeAndRanks> nodeAndIndexList;
  //Serial Merge
  {
    ot::TreeNode currNode = nodeAndElemIdList[0].node;
    ot::NodeAndRanks tmpObj;
    tmpObj.node = currNode;
    tmpObj.node.setWeight(1);
    tmpObj.ranks.push_back(rankList.size());
    nodeAndIndexList.push_back(tmpObj);
    rankList.push_back(currNode.getWeight());
    elemIdList.push_back(nodeAndElemIdList[0].values[0]);
  }
  for(size_t i = 1; i < nodeAndElemIdList.size(); ++i) {
    ot::TreeNode currNode = nodeAndElemIdList[i].node;
    if( (nodeAndIndexList[nodeAndIndexList.size() - 1].node == currNode) ||
        (nodeAndIndexList[nodeAndIndexList.size() - 1].node.isAncestor(currNode)) ) {
      nodeAndIndexList[nodeAndIndexList.size() - 1].ranks.push_back(rankList.size());
    } else {
      ot::NodeAndRanks tmpObj;
      tmpObj.node = currNode;
      tmpObj.node.setWeight(1);
      tmpObj.ranks.push_back(rankList.size());
      nodeAndIndexList.push_back(tmpObj);
    }
    rankList.push_back(currNode.getWeight());
    elemIdList.push_back(nodeAndElemIdList[i].values[0]);
  }//end i

  int npes = globalComm.getSize();

  std::vector< ot::TreeNode > firstAndLastList(2*npes);
  firstAndLastList[(2*rank)] = nodeAndIndexList[0].node;
  firstAndLastList[(2*rank) + 1] = nodeAndIndexList[nodeAndIndexList.size() - 1].node;

  MPI_Allgather( &(firstAndLastList[(2*rank)]), 2, par::Mpi_datatype<ot::TreeNode>::value(), 
      &(firstAndLastList[0]), 2, par::Mpi_datatype<ot::TreeNode>::value(), globalComm.getCommunicator() );

  int numToSend = 0;
  for(int i = (rank + 1); i < npes; ++i) {
    if( (firstAndLastList[(2*rank) + 1] == firstAndLastList[2*i]) ||
        (firstAndLastList[(2*rank) + 1].isAncestor(firstAndLastList[2*i])) ) {
      ++numToSend;
    } else {
      break;
    }
  }//end i

  int rankOfSending = -1; 
  for(int i = 0; i < rank; ++i) {
    if( (firstAndLastList[(2*i) + 1] == firstAndLastList[2*rank]) ||
        (firstAndLastList[(2*i) + 1].isAncestor(firstAndLastList[2*rank])) ) {
      rankOfSending = i;
      break;
    }
  }//end i

  //If global max(numToSend) is small then we can use Point2Point can be used instead of the
  //following All2allv.   
  //We can do a Reduce on numToSend to see if this is the case.

  int *sendCnts = new int[npes];
  int *recvCnts = new int[npes];
  int *sendDisps = new int[npes];
  int *recvDisps = new int[npes];

  for(int i = 0; i < npes; ++i) {
    sendCnts[i] = 0;
    recvCnts[i] = 0;
  }//end i
  if(rankOfSending >= 0) {
    recvCnts[rankOfSending] = 1;
  }
  for(int i = 1; i <= numToSend; ++i) {
    sendCnts[rank + i] = 1;
  }//end i

  for(int i = 0; i < npes; ++i) {
    sendDisps[i] = 0;
    recvDisps[i] = 0;
  }//end i

  ot::TreeNode remoteElem;
  MPI_Alltoallv( (&(firstAndLastList[(2*rank) + 1])), sendCnts, sendDisps, par::Mpi_datatype<ot::TreeNode>::value(),
      (&remoteElem), recvCnts, recvDisps, par::Mpi_datatype<ot::TreeNode>::value(), globalComm.getCommunicator() );

  firstAndLastList.clear();



  delete [] sendCnts;
  delete [] recvCnts;
  delete [] sendDisps;
  delete [] recvDisps;

  //TO DO: Form Mins List (Partition Info) for nodeAndIndexList.


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



