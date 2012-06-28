
#include "ampmesh/dendro/DendroSearch.h"

DendroSearch::DendroSearch(AMP::AMP_MPI comm, AMP::Mesh::Mesh::shared_ptr mesh) 
  : d_globalComm(comm), d_meshAdapter(mesh) {
    d_verbose = true;
    d_rank = d_globalComm.getRank();
    d_npes = d_globalComm.getSize();
    d_minCoords.resize(3);
    d_maxCoords.resize(3);
    d_scalingFactor.resize(3);
    setupDendro();
  }

void DendroSearch::interpolate(AMP::LinearAlgebra::Vector::shared_ptr vectorField, const unsigned int dofsPerNode,
    const std::vector<double> & pts, std::vector<double> & results, std::vector<bool> & foundPt) {
  search(pts);
  interpolate(vectorField, dofsPerNode, results, foundPt); 
}

void DendroSearch::setupDendro() {
  double setupBeginTime, setupEndTime;
  if(d_verbose) {
    d_globalComm.barrier();
    setupBeginTime = MPI_Wtime();
  }

  setupDSforSearch();

  if(d_verbose) {
    d_globalComm.barrier();
    setupEndTime = MPI_Wtime();
    if(!d_rank) {
      std::cout<<"Finished setting up DS for search in "<<(setupEndTime - setupBeginTime)<<" seconds."<<std::endl;
    }
  }
}

void DendroSearch::setupDSforSearch() {
  std::vector<double> box = d_meshAdapter->getBoundingBox();
  for(int i = 0; i < d_meshAdapter->getDim(); ++i) {
    d_minCoords[i] = box[2*i+0];
    d_maxCoords[i] = box[2*i+1];
    d_scalingFactor[i] = 1.0/(1.0e-10 + d_maxCoords[i] - d_minCoords[i]);
  }

  createLocalMeshElementArray();

  const unsigned int MaxDepth = 30;

  unsigned int totalNumElems = d_meshAdapter->numGlobalElements(AMP::Mesh::Volume);

  double avgHboxInv = std::pow(totalNumElems, (1.0/3.0));
  AMP_ASSERT(avgHboxInv > 1.0);
  d_boxLevel = binOp::fastLog2(static_cast<unsigned int>(std::ceil(avgHboxInv)));
  AMP_ASSERT(d_boxLevel < MaxDepth);

  if(d_verbose) {
    if(!d_rank) {
      std::cout<<"BoxLevel = "<<d_boxLevel<<std::endl;
    }
  }

  const double hBox = 1.0/(static_cast<double>(1u << d_boxLevel));

  std::vector< ot::NodeAndValues<int, 1> > nodeAndElemIdList;

  AMP_ASSERT(!(d_localElemArr.empty()));

  for(int eId = 0; eId < d_localElemArr.size(); ++eId) {
    std::vector<AMP::Mesh::MeshElement> currNodes = d_localElemArr[eId].getElements(AMP::Mesh::Vertex);
    int minId[3];
    int maxId[3];
    for(size_t i = 0; i < currNodes.size(); ++i) {
      std::vector<double> pt = currNodes[i].coord();
      double scaledPt[3];
      for(int j = 0; j < 3; ++j) {
        scaledPt[j] = ((pt[j] - d_minCoords[j])*d_scalingFactor[j]);
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
          unsigned int bX = i*(1u << (MaxDepth - d_boxLevel));
          unsigned int bY = j*(1u << (MaxDepth - d_boxLevel));
          unsigned int bZ = k*(1u << (MaxDepth - d_boxLevel));
          ot::TreeNode box(bX, bY, bZ, d_boxLevel, 3, MaxDepth);
          box.setWeight(d_rank);
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
      nodeAndElemIdList, tmpList, (d_globalComm.getCommunicator()));
  swap(nodeAndElemIdList, tmpList);
  tmpList.clear();

  int numLocalOcts = nodeAndElemIdList.size();

  if(d_verbose) {
    int numGlobalOcts = d_globalComm.sumReduce<int>(numLocalOcts);
    if(!d_rank) {
      std::cout<<"Total num initial octants = "<<numGlobalOcts <<std::endl;
    }
  }

  AMP_ASSERT(d_rankList.empty());
  AMP_ASSERT(d_elemIdList.empty());

  for(int i = 0; i < numLocalOcts; ++i) {
    ot::TreeNode currNode = nodeAndElemIdList[i].node;
    d_rankList.push_back(currNode.getWeight());
    d_elemIdList.push_back(nodeAndElemIdList[i].values[0]);
  }//end i

  AMP_ASSERT(numLocalOcts > 0);
  AMP_ASSERT(d_nodeList.empty());

  //Local Merge
  {
    ot::TreeNode currNode = nodeAndElemIdList[0].node;
    currNode.setWeight(1);
    d_nodeList.push_back(currNode);
  }
  for(size_t i = 1; i < nodeAndElemIdList.size(); ++i) {
    ot::TreeNode currNode = nodeAndElemIdList[i].node;
    if( d_nodeList[d_nodeList.size() - 1] == currNode ) {
      d_nodeList[d_nodeList.size() - 1].addWeight(1);
    } else {
      currNode.setWeight(1);
      d_nodeList.push_back(currNode);
    }
  }//end i
  nodeAndElemIdList.clear();

  int localFlag = 0;
  if( (d_rank > 0) && (d_rank < (d_npes - 1)) && ((d_nodeList.size()) == 1) ) {
    localFlag = 1;
  }

  int globalFlag = d_globalComm.sumReduce(localFlag);

  int prevRank = d_rank - 1;
  int nextRank = d_rank + 1;

  if(globalFlag > 0) {
    int gatherSendBuf = 0;
    if( (d_rank > 0) && (d_rank < (d_npes - 1)) && (d_nodeList.size() == 1) ) {
      gatherSendBuf = d_rankList.size();
    }

    int* gatherList = new int[d_npes];

    d_globalComm.allGather(gatherSendBuf, gatherList);

    if(d_rank > 0) {
      while(gatherList[prevRank] > 0) {
        --prevRank;
      }//end while
    }

    if(d_rank < (d_npes - 1)) {
      while(gatherList[nextRank] > 0) {
        ++nextRank;
      }//end while
    }

    int* sendBoxCnts = new int[d_npes];
    int* recvBoxCnts = new int[d_npes];

    int* sendSourceCnts = new int[d_npes];
    int* recvSourceCnts = new int[d_npes];

    for(int i = 0; i < d_npes; ++i) {
      sendBoxCnts[i] = 0;
      recvBoxCnts[i] = 0;
      sendSourceCnts[i] = 0;
      recvSourceCnts[i] = 0;
    }//end i

    if(gatherSendBuf > 0) {
      sendBoxCnts[prevRank] = 1;
      sendSourceCnts[prevRank] = gatherSendBuf;
    }
    for(int i = d_rank + 1; i < nextRank; ++i) {
      recvBoxCnts[i] = 1;
      recvSourceCnts[i] = gatherList[i];
    }//end i

    delete [] gatherList;

    int* sendBoxDisps = new int[d_npes];
    int* recvBoxDisps = new int[d_npes];
    sendBoxDisps[0] = 0;
    recvBoxDisps[0] = 0;
    for(int i = 1; i < d_npes; ++i) {
      sendBoxDisps[i] = sendBoxDisps[i - 1] + sendBoxCnts[i - 1];
      recvBoxDisps[i] = recvBoxDisps[i - 1] + recvBoxCnts[i - 1];
    }//end i

    std::vector<ot::TreeNode> tmpBoxList(recvBoxDisps[d_npes - 1] + recvBoxCnts[d_npes - 1]);

    ot::TreeNode* recvBoxBuf = NULL;
    if(!(tmpBoxList.empty())) {
      recvBoxBuf = (&(tmpBoxList[0]));
    }

    d_globalComm.allToAll( (&(d_nodeList[0])), sendBoxCnts, sendBoxDisps, recvBoxBuf, recvBoxCnts, recvBoxDisps, true);

    if(gatherSendBuf > 0) {
      d_nodeList.clear();
    } else {
      for(int i = 0; i < tmpBoxList.size(); ++i) {
        if(tmpBoxList[i] == d_nodeList[d_nodeList.size() - 1]) {
          d_nodeList[d_nodeList.size() - 1].addWeight(tmpBoxList[i].getWeight());
        } else {
          d_nodeList.push_back(tmpBoxList[i]);
        }
      }//end i
    }

    delete [] sendBoxCnts;
    delete [] recvBoxCnts;
    delete [] sendBoxDisps;
    delete [] recvBoxDisps;

    int* sendSourceDisps = new int[d_npes];
    int* recvSourceDisps = new int[d_npes];
    sendSourceDisps[0] = 0;
    recvSourceDisps[0] = 0;
    for(int i = 1; i < d_npes; ++i) {
      sendSourceDisps[i] = sendSourceDisps[i - 1] + sendSourceCnts[i - 1];
      recvSourceDisps[i] = recvSourceDisps[i - 1] + recvSourceCnts[i - 1];
    }//end i

    std::vector<int> tmpRankList(recvSourceDisps[d_npes - 1] + recvSourceCnts[d_npes - 1]);
    std::vector<int> tmpElemIdList(recvSourceDisps[d_npes - 1] + recvSourceCnts[d_npes - 1]);

    int* recvRankBuf = NULL;
    int* recvElemIdBuf = NULL;
    if(!(tmpRankList.empty())) {
      recvRankBuf = (&(tmpRankList[0]));
      recvElemIdBuf = (&(tmpElemIdList[0]));
    }

    d_globalComm.allToAll( (&(d_rankList[0])), sendSourceCnts, sendSourceDisps, recvRankBuf, recvSourceCnts, recvSourceDisps, true);
    d_globalComm.allToAll( (&(d_elemIdList[0])), sendSourceCnts, sendSourceDisps, recvElemIdBuf, recvSourceCnts, recvSourceDisps, true);

    if(gatherSendBuf > 0) {
      d_rankList.clear();
      d_elemIdList.clear();
    } else {
      if(!(tmpRankList.empty())) {
        d_rankList.insert(d_rankList.end(), tmpRankList.begin(), tmpRankList.end());
        d_elemIdList.insert(d_elemIdList.end(), tmpElemIdList.begin(), tmpElemIdList.end());
      }
    }

    delete [] sendSourceCnts;
    delete [] recvSourceCnts;
    delete [] sendSourceDisps;
    delete [] recvSourceDisps;
  }

  if(!(d_nodeList.empty())) {
    AMP_ASSERT(d_nodeList.size() >= 2);

    ot::TreeNode prevBox;
    ot::TreeNode nextBox;
    ot::TreeNode firstBox = d_nodeList[0];
    ot::TreeNode lastBox = d_nodeList[d_nodeList.size() - 1];
    MPI_Request recvPrevReq;
    MPI_Request recvNextReq;
    MPI_Request sendFirstReq;
    MPI_Request sendLastReq;
    if(d_rank > 0) {
      MPI_Irecv(&prevBox, 1, par::Mpi_datatype<ot::TreeNode>::value(),
          prevRank, 1, (d_globalComm.getCommunicator()), &recvPrevReq);
      MPI_Isend(&firstBox, 1, par::Mpi_datatype<ot::TreeNode>::value(),
          prevRank, 2, (d_globalComm.getCommunicator()), &sendFirstReq);
    }
    if(d_rank < (d_npes - 1)) {
      MPI_Irecv(&nextBox, 1, par::Mpi_datatype<ot::TreeNode>::value(),
          nextRank, 2, (d_globalComm.getCommunicator()), &recvNextReq);
      MPI_Isend(&lastBox, 1, par::Mpi_datatype<ot::TreeNode>::value(),
          nextRank, 1, (d_globalComm.getCommunicator()), &sendLastReq);
    }

    if(d_rank > 0) {
      MPI_Status status;
      MPI_Wait(&recvPrevReq, &status);
      MPI_Wait(&sendFirstReq, &status);
    }
    if(d_rank < (d_npes - 1)) {
      MPI_Status status;
      MPI_Wait(&recvNextReq, &status);
      MPI_Wait(&sendLastReq, &status);
    }

    bool removeFirst = false;
    bool addToLast = false;
    if(d_rank > 0) {
      if(prevBox == firstBox) {
        removeFirst = true;
      }
    }
    if(d_rank < (d_npes - 1)) {
      if(nextBox == lastBox) {
        addToLast = true;
      }
    }

    MPI_Request recvRankReq;
    MPI_Request recvElemIdReq;
    if(addToLast) {
      int numPts = d_rankList.size();
      d_rankList.resize(numPts + (nextBox.getWeight()));
      d_elemIdList.resize(numPts + (nextBox.getWeight()));
      d_nodeList[d_nodeList.size() - 1].addWeight(nextBox.getWeight());
      MPI_Irecv((&(d_rankList[numPts])), ((nextBox.getWeight())), MPI_INT, nextRank,
          3, (d_globalComm.getCommunicator()), &recvRankReq);
      MPI_Irecv((&(d_elemIdList[numPts])), ((nextBox.getWeight())), MPI_INT, nextRank,
          4, (d_globalComm.getCommunicator()), &recvElemIdReq);
    }
    if(removeFirst) {
      MPI_Send((&(d_rankList[0])), ((firstBox.getWeight())), MPI_INT, prevRank, 3, (d_globalComm.getCommunicator()));
      MPI_Send((&(d_elemIdList[0])), ((firstBox.getWeight())), MPI_INT, prevRank, 4, (d_globalComm.getCommunicator()));
      d_nodeList.erase(d_nodeList.begin());
    }
    if(addToLast) {
      MPI_Status status;
      MPI_Wait(&recvRankReq, &status);
      MPI_Wait(&recvElemIdReq, &status);
    }
    if(removeFirst) {
      d_rankList.erase(d_rankList.begin(), d_rankList.begin() + ((firstBox.getWeight())));
      d_elemIdList.erase(d_elemIdList.begin(), d_elemIdList.begin() + ((firstBox.getWeight())));
    }

    d_stIdxList.resize(d_nodeList.size());

    d_stIdxList[0] = 0;
    for(int i = 1; i < d_nodeList.size(); ++i) {
      d_stIdxList[i] = d_stIdxList[i - 1] + d_nodeList[i - 1].getWeight();
    }//end i
  }

  ot::TreeNode firstNode;
  if(!(d_nodeList.empty())) {
    firstNode = d_nodeList[0];
    firstNode.setWeight(d_rank);
  }
  d_mins.resize(d_npes);
  d_globalComm.allGather(firstNode, &(d_mins[0]));

  std::vector<ot::TreeNode> tmpMins;
  for(int i = 0; i < d_npes; ++i) {
    if(d_mins[i].getDim() > 0) {
      tmpMins.push_back(d_mins[i]);
    }
  }//end i
  swap(d_mins, tmpMins);
  tmpMins.clear();

  if(d_verbose) {
    int minFineListLen = d_nodeList[0].getWeight();
    int maxFineListLen = d_nodeList[0].getWeight();
    for(int i = 1; i < d_nodeList.size(); ++i) {
      if(minFineListLen > d_nodeList[i].getWeight()) {
        minFineListLen = d_nodeList[i].getWeight();
      }
      if(maxFineListLen < d_nodeList[i].getWeight()) {
        maxFineListLen = d_nodeList[i].getWeight();
      }
    }//end i

    int globalMinFineListLen = d_globalComm.minReduce(minFineListLen);
    int globalMaxFineListLen = d_globalComm.maxReduce(maxFineListLen);

    numLocalOcts = d_nodeList.size();

    int numGlobalOcts = d_globalComm.sumReduce(numLocalOcts);

    if(!d_rank) {
      std::cout<<"Total num final octants = "<<numGlobalOcts <<std::endl;
      std::cout<<"Global Min Fine List Length = "<<globalMinFineListLen <<std::endl;
      std::cout<<"Global Max Fine List Length = "<<globalMaxFineListLen <<std::endl;
    }
  }
}

void DendroSearch::search(const std::vector<double> & pts) {
  d_numLocalPts = (pts.size())/3;

  double searchBeginTime, searchStep1Time, searchStep2Time, searchStep3Time, searchStep4Time, searchStep5Time, searchStep6Time;
  if(d_verbose) {
    d_globalComm.barrier();
    searchBeginTime = MPI_Wtime();
  }

  const unsigned int MaxDepth = 30;
  const unsigned int ITPMD = (1u << MaxDepth);
  const double DTPMD = static_cast<double>(ITPMD);

  std::vector<ot::NodeAndValues<double, 4> > ptsWrapper(d_numLocalPts);
  for(int i = 0; i < d_numLocalPts; ++i) {
    double x = pts[3*i];
    double y = pts[(3*i) + 1];
    double z = pts[(3*i) + 2];
    double scaledX = ((x - d_minCoords[0])*d_scalingFactor[0]);
    double scaledY = ((y - d_minCoords[1])*d_scalingFactor[1]);
    double scaledZ = ((z - d_minCoords[2])*d_scalingFactor[2]);
    unsigned int pX = static_cast<unsigned int>(scaledX*DTPMD);
    unsigned int pY = static_cast<unsigned int>(scaledY*DTPMD);
    unsigned int pZ = static_cast<unsigned int>(scaledZ*DTPMD);

    ptsWrapper[i].node =  ot::TreeNode(pX, pY, pZ, MaxDepth, 3, MaxDepth);
    ptsWrapper[i].node.setWeight(d_rank);
    ptsWrapper[i].values[0] = x;
    ptsWrapper[i].values[1] = y;
    ptsWrapper[i].values[2] = z;
    ptsWrapper[i].values[3] = i;     
  }//end i

  if(d_verbose) {
    d_globalComm.barrier();
    searchStep1Time = MPI_Wtime();
    if(!d_rank) {
      std::cout<<"Time for step-1 of search: "<<(searchStep1Time - searchBeginTime)<<" seconds."<<std::endl;
    }
  }

  //Performance Question: Should PtsWrapper be sorted or not?
  //If PtsWrapper is sorted (even just a local sort), we can skip the
  //binary searches and use binning instead.  Binning is amortized constant
  //time and using binary searches would be logarithmic. This is just a matter
  //of constants since sorting is also logarithmic.

  d_sendCnts.resize(d_npes);
  d_recvCnts.resize(d_npes);
  d_sendDisps.resize(d_npes);
  d_recvDisps.resize(d_npes);
  std::fill(d_sendCnts.begin(), d_sendCnts.end(), 0);

  std::vector<int> part(d_numLocalPts, -1);

  for(int i = 0; i < d_numLocalPts; ++i) {
    unsigned int retIdx;
    bool found = seq::maxLowerBound<ot::TreeNode>(d_mins, (ptsWrapper[i].node), retIdx, NULL, NULL);
    if(found) {
      part[i] = d_mins[retIdx].getWeight();
      ++(d_sendCnts[part[i]]);
    }
  }//end i

  if(d_verbose) {
    d_globalComm.barrier();
    searchStep2Time = MPI_Wtime();
    if(!d_rank) {
      std::cout<<"Time for step-2 of search: "<<(searchStep2Time - searchStep1Time)<<" seconds."<<std::endl;
    }
  }

  d_globalComm.allToAll(1, &(d_sendCnts[0]), &(d_recvCnts[0]));

  d_sendDisps[0] = 0;
  d_recvDisps[0] = 0;
  for(int i = 1; i < d_npes; ++i) {
    d_sendDisps[i] = d_sendDisps[i - 1] + d_sendCnts[i - 1];
    d_recvDisps[i] = d_recvDisps[i - 1] + d_recvCnts[i - 1];
  }//end i

  std::vector<ot::NodeAndValues<double, 4> > sendList(d_sendDisps[d_npes - 1] + d_sendCnts[d_npes - 1]);

  std::fill(d_sendCnts.begin(), d_sendCnts.end(), 0);

  for(int i = 0; i < d_numLocalPts; ++i) {
    if(part[i] >= 0) {
      sendList[d_sendDisps[part[i]] + d_sendCnts[part[i]]] = ptsWrapper[i];
      ++(d_sendCnts[part[i]]);
    }
  }//end i
  ptsWrapper.clear();

  std::vector<ot::NodeAndValues<double, 4> > recvList(d_recvDisps[d_npes - 1] + d_recvCnts[d_npes - 1]);
  d_globalComm.allToAll((!(sendList.empty()) ? &(sendList[0]) : NULL), &(d_sendCnts[0]), &(d_sendDisps[0]),
      (!(recvList.empty()) ? &(recvList[0]) : NULL), &(d_recvCnts[0]), &(d_recvDisps[0]), true);
  sendList.clear();

  if(d_verbose) {
    d_globalComm.barrier();
    searchStep3Time = MPI_Wtime();
    if(!d_rank) {
      std::cout<<"Time for step-3 of search: "<<(searchStep3Time - searchStep2Time)<<" seconds."<<std::endl;
    }
  }

  std::fill(d_sendCnts.begin(), d_sendCnts.end(), 0);

  std::vector<int> ptToOctMap((recvList.size()), -1);
  for(int i = 0; i < recvList.size(); ++i) {
    unsigned int retIdx;
    seq::maxLowerBound<ot::TreeNode>(d_nodeList, (recvList[i].node), retIdx, NULL, NULL);
    if( d_nodeList[retIdx].isAncestor(recvList[i].node) ) {
      ptToOctMap[i] = retIdx;
      int stIdx = d_stIdxList[retIdx];
      for(int j = 0; j < d_nodeList[retIdx].getWeight(); ++j) {
        d_sendCnts[d_rankList[stIdx + j]]++;
      }//end j
    }
  }//end i

  if(d_verbose) {
    d_globalComm.barrier();
    searchStep4Time = MPI_Wtime();
    if(!d_rank) {
      std::cout<<"Time for step-4 of search: "<<(searchStep4Time - searchStep3Time)<<" seconds."<<std::endl;
    }
  }

  d_globalComm.allToAll(1, &(d_sendCnts[0]), &(d_recvCnts[0]));

  d_sendDisps[0] = 0;
  d_recvDisps[0] = 0;
  for(int i = 1; i < d_npes; ++i) {
    d_sendDisps[i] = d_sendDisps[i - 1] + d_sendCnts[i - 1];
    d_recvDisps[i] = d_recvDisps[i - 1] + d_recvCnts[i - 1];
  }//end i

  std::vector<double> sendPtsList(6*(d_sendDisps[d_npes - 1] + d_sendCnts[d_npes - 1]));

  std::fill(d_sendCnts.begin(), d_sendCnts.end(), 0);

  for(int i = 0; i < ptToOctMap.size(); ++i) {
    if(ptToOctMap[i] >= 0) {
      int stIdx = d_stIdxList[ptToOctMap[i]];
      for(int j = 0; j < d_nodeList[ptToOctMap[i]].getWeight(); ++j) {
        int recvRank = d_rankList[stIdx + j];
        int currIdx = 6*(d_sendDisps[recvRank] + d_sendCnts[recvRank]);
        //Local Id of this element on the processor that owns this element
        sendPtsList[currIdx] = d_elemIdList[stIdx + j];
        //Pt's x coordinate
        sendPtsList[currIdx + 1] = recvList[i].values[0];
        //Pt's y coordinate
        sendPtsList[currIdx + 2] = recvList[i].values[1];
        //Pt's z coordinate
        sendPtsList[currIdx + 3] = recvList[i].values[2];
        //Local Id of Pt on the processor that owns this Pt
        sendPtsList[currIdx + 4] = recvList[i].values[3];
        //rank of processor that owns Pt
        sendPtsList[currIdx + 5] = recvList[i].node.getWeight();
        ++(d_sendCnts[recvRank]);
      }//end j
    }
  }//end i
  recvList.clear();

  for(int i = 0; i < d_npes; ++i) {
    d_sendCnts[i] *= 6;
    d_sendDisps[i] *= 6;
    d_recvCnts[i] *= 6;
    d_recvDisps[i] *= 6;
  }//end i

  std::vector<double> recvPtsList(d_recvDisps[d_npes - 1] + d_recvCnts[d_npes - 1]);
  d_globalComm.allToAll((!(sendPtsList.empty()) ? &(sendPtsList[0]) : NULL), &(d_sendCnts[0]), &(d_sendDisps[0]), 
      (!(recvPtsList.empty()) ? &(recvPtsList[0]) : NULL), &(d_recvCnts[0]), &(d_recvDisps[0]), true);
  sendPtsList.clear();

  if(d_verbose) {
    d_globalComm.barrier();
    searchStep5Time = MPI_Wtime();
    if(!d_rank) {
      std::cout<<"Time for step-5 of search: "<<(searchStep5Time - searchStep4Time)<<" seconds."<<std::endl;
    }
  }

  unsigned int n_volume_elements = d_localElemArr.size();
  std::vector<hex8_element_t> volume_elements;
  volume_elements.reserve(n_volume_elements);
  for (unsigned int i = 0; i < n_volume_elements; ++i) {
    AMP::Mesh::MeshElement* amp_element = &(d_localElemArr[i]);
    std::vector<AMP::Mesh::MeshElement> amp_vector_support_points = amp_element->getElements(AMP::Mesh::Vertex);
    AMP_ASSERT(amp_vector_support_points.size() == 8);
    std::vector<double> support_points(24);
    for (unsigned int j = 0; j < 8; ++j) {
      std::vector<double> point_coord = amp_vector_support_points[j].coord();
      support_points[3*j+0] = point_coord[0];
      support_points[3*j+1] = point_coord[1];
      support_points[3*j+2] = point_coord[2];
    } // end j
    volume_elements.push_back(hex8_element_t(support_points));
  } // end for i

  std::fill(d_sendCnts.begin(), d_sendCnts.end(), 0);

  int numRecvPts = recvPtsList.size()/6;
  std::vector<double> tmpPtLocalCoord(3);

  d_foundPts.reserve(6*numRecvPts);
  unsigned int numFoundPts = 0;
  bool coordinates_are_local = true;
  for(int i = 0; i < numRecvPts; ++i) {
    double const * tmpPtGlobalCoordPtr = &(recvPtsList[6*i])+1;
    unsigned int eId = static_cast<unsigned int>(recvPtsList[6*i]);
    unsigned int procId = static_cast<unsigned int>(recvPtsList[6*i+5]);

    if (volume_elements[eId].within_bounding_box(tmpPtGlobalCoordPtr)) {
      if (volume_elements[eId].within_bounding_polyhedron(tmpPtGlobalCoordPtr)) {
        volume_elements[eId].map_global_to_local(tmpPtGlobalCoordPtr, &(tmpPtLocalCoord[0]));
        if (volume_elements[eId].contains_point(&(tmpPtLocalCoord[0]), coordinates_are_local)) {
          d_foundPts.push_back(recvPtsList[6*i]);
          for (unsigned int d = 0; d < 3; ++d) { d_foundPts.push_back(tmpPtLocalCoord[d]); }
          d_foundPts.push_back(recvPtsList[6*i+4]);
          d_foundPts.push_back(recvPtsList[6*i+5]);
          ++numFoundPts;
          ++(d_sendCnts[procId]);
        } // end if
      } // end if
    } // end if
  }//end i
  recvPtsList.clear();

  d_globalComm.allToAll(1, &(d_sendCnts[0]), &(d_recvCnts[0]));

  d_sendDisps[0] = 0;
  d_recvDisps[0] = 0;
  for(int i = 1; i < d_npes; ++i) {
    d_sendDisps[i] = d_sendDisps[i - 1] + d_sendCnts[i - 1];
    d_recvDisps[i] = d_recvDisps[i - 1] + d_recvCnts[i - 1];
  }//end i

  if(d_verbose) {
    d_globalComm.barrier();
    searchStep6Time = MPI_Wtime();
    if(!d_rank) {
      std::cout<<"Time for step-6 of search: "<<(searchStep6Time - searchStep5Time)<<" seconds."<<std::endl;
    }
  }
}

void DendroSearch::interpolate(AMP::LinearAlgebra::Vector::shared_ptr vectorField, const unsigned int dofsPerNode,
    std::vector<double> & results, std::vector<bool> & foundPt) {
  double interpolateBeginTime, interpolateStep1Time, interpolateStep2Time;
  if(d_verbose) {
    d_globalComm.barrier();
    interpolateBeginTime = MPI_Wtime();
  }

  vectorField->makeConsistent(  AMP::LinearAlgebra::Vector::CONSISTENT_SET );
  AMP::Discretization::DOFManager::shared_ptr dofManager = vectorField->getDOFManager();

  for(unsigned int i = 0; i < d_npes; ++i) {
    d_sendCnts[i] *= (dofsPerNode + 1);
    d_recvCnts[i] *= (dofsPerNode + 1);
    d_sendDisps[i] *= (dofsPerNode + 1);
    d_recvDisps[i] *= (dofsPerNode + 1);
  }// end for i

  std::vector<double> sendResults(d_sendDisps[d_npes - 1] + d_sendCnts[d_npes - 1]);

  std::vector<int> tmpSendCnts(d_npes, 0);

  std::vector<double> basis_functions_values(8);
  for(int i = 0; i < d_foundPts.size(); i += 6) {
    AMP::Mesh::MeshElement* amp_element = &(d_localElemArr[static_cast<unsigned int>(d_foundPts[i])]);
    std::vector<AMP::Mesh::MeshElement> amp_vector_support_points = amp_element->getElements(AMP::Mesh::Vertex);
    get_basis_functions_values(&(d_foundPts[i + 1]), &(basis_functions_values[0]));

    std::vector<double> value(dofsPerNode, 0.0);
    for (unsigned int j = 0; j < 8; ++j) {
      std::vector<size_t> globalID;
      dofManager->getDOFs(amp_vector_support_points[j].globalID(), globalID);
      AMP_ASSERT(globalID.size() == dofsPerNode);
      for(int d = 0; d < dofsPerNode; ++d) {
        double vecVal = vectorField->getValueByGlobalID(globalID[d]);
        value[d] += (vecVal * basis_functions_values[j]);
      }//end d
    } // end j
    unsigned int ptProcId = static_cast<unsigned int>(d_foundPts[i + 5]);
    sendResults[d_sendDisps[ptProcId] + tmpSendCnts[ptProcId]] = d_foundPts[i + 4];
    ++(tmpSendCnts[ptProcId]);
    for(int d = 0; d < dofsPerNode; ++d) {
      sendResults[d_sendDisps[ptProcId] + tmpSendCnts[ptProcId]] = value[d];
      ++(tmpSendCnts[ptProcId]);
    }//end d   
  }//end i
  tmpSendCnts.clear();

  if(d_verbose) {
    d_globalComm.barrier();
    interpolateStep1Time = MPI_Wtime();
    if(!d_rank) {
      std::cout<<"Time for step-1 of interpolate: "<<(interpolateStep1Time - interpolateBeginTime)<<" seconds."<<std::endl;
    }
  }

  std::vector<double> recvResults(d_recvDisps[d_npes - 1] + d_recvCnts[d_npes - 1]);

  d_globalComm.allToAll((!(sendResults.empty()) ? &(sendResults[0]) : NULL), &(d_sendCnts[0]), &(d_sendDisps[0]),
      (!(recvResults.empty()) ? &(recvResults[0]) : NULL), &(d_recvCnts[0]), &(d_recvDisps[0]), true);
  sendResults.clear();

  //Points that are not found will have a result = 0. 
  results.resize(dofsPerNode*d_numLocalPts);
  for (unsigned int i = 0; i < results.size(); ++i) {
    results[i] = 0.0;
  } // end for i

  foundPt.resize(d_numLocalPts);
  for (unsigned int i = 0; i < foundPt.size(); ++i) {
    foundPt[i] = false;
  } // end for i

  for(size_t i = 0; i < recvResults.size(); i += (dofsPerNode + 1)) {
    unsigned int locId = static_cast<unsigned int>(recvResults[i]);
    foundPt[locId] = true;
    for(int d = 0; d < dofsPerNode; ++d) {
      results[(locId*dofsPerNode) + d] = recvResults[i + d + 1];
    }//end d
  }//end i

  for(unsigned int i = 0; i < d_npes; ++i) {
    d_sendCnts[i] /= (dofsPerNode + 1);
    d_recvCnts[i] /= (dofsPerNode + 1);
    d_sendDisps[i] /= (dofsPerNode + 1);
    d_recvDisps[i] /= (dofsPerNode + 1);
  } // end for i

  if(d_verbose) {
    d_globalComm.barrier();
    interpolateStep2Time = MPI_Wtime();
    if(!d_rank) {
      std::cout<<"Time for step-2 of interpolate: "<<(interpolateStep2Time - interpolateStep1Time)<<" seconds."<<std::endl;
    }
  }
}

void DendroSearch::createLocalMeshElementArray() {
  d_localElemArr.clear();
  AMP::Mesh::MeshIterator el = d_meshAdapter->getIterator(AMP::Mesh::Volume, 0);
  AMP::Mesh::MeshIterator end_el = el.end();
  AMP_ASSERT(el != end_el);
  for(; el != end_el; ++el) {
    d_localElemArr.push_back(*el);
  }//end el
}



