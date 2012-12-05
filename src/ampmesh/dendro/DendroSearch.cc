
#include "ampmesh/dendro/DendroSearch.h"

#include <numeric>
#include <fstream>
#include <iomanip>
#include <boost/lexical_cast.hpp>
#include <set>
#include "ampmesh/latex_visualization_tools.h"

namespace AMP {
  namespace Mesh {

    DendroSearch::DendroSearch(AMP::Mesh::Mesh::shared_ptr mesh, bool verbose, std::ostream & oStream)
      : d_meshAdapter(mesh),
      d_verbose(verbose),
      d_oStream(oStream),
      d_timingMeasurements(std::vector<double>(numTimingTypes, -1.0)),
      d_tolerance(1.0e-12) {
        d_minCoords.resize(3);
        d_scalingFactor.resize(3);
        setupDSforSearch();
      }

    void DendroSearch::searchAndInterpolate(AMP::AMP_MPI comm, AMP::LinearAlgebra::Vector::const_shared_ptr vectorField, 
        const unsigned int dofsPerNode, const std::vector<double> & pts, std::vector<double> & results, std::vector<bool> & foundPt) {
      search(comm, pts);
      interpolate(comm, vectorField, dofsPerNode, results, foundPt); 
    }

    void DendroSearch::setTolerance(double tolerance) {
      d_tolerance = tolerance;
    }

    void DendroSearch::projectOnBoundaryID(AMP::AMP_MPI comm, const int boundaryID, std::vector<AMP::Mesh::MeshElementID> & faceVerticesGlobalIDs, 
        std::vector<double> & shiftGlobalCoords, std::vector<double> & projectionLocalCoordsOnFace, std::vector<int> & flags) {
      const int rank = comm.getRank();
      const int npes = comm.getSize();

      double projectBeginTime, projectStep1Time=0., projectStep2Time;
      if(d_verbose) {
        comm.barrier();
      }
      projectBeginTime = MPI_Wtime();

      std::vector<ProjectOnBoundaryData> sendData(d_sendDisps[npes-1] + d_sendCnts[npes-1]);

      std::vector<int> tmpSendCnts(npes, 0);

      for (unsigned int i = 0; i < d_foundPts.size(); i += 6) {
        ProjectOnBoundaryData tmpData;
        const double * pointLocalCoords_ptr = &(d_foundPts[i+1]);
        const size_t pointLocalID = static_cast<size_t>(d_foundPts[i+4]);
        const size_t pointOwnerRank = static_cast<size_t>(d_foundPts[i+5]);
        const size_t elementLocalID = static_cast<size_t>(d_foundPts[i]);
        tmpData.d_PointLocalID = pointLocalID;
        if (d_localElems[elementLocalID].isOnBoundary(boundaryID)) { // point was found and element is on boundary
          std::vector<AMP::Mesh::MeshElement> meshElementFaces = d_localElems[elementLocalID].getElements(AMP::Mesh::Face);
          AMP_CHECK_ASSERT( meshElementFaces.size() == 6 );
          for (size_t f = 0; f < 6; ++f) {
            if (meshElementFaces[f].isOnBoundary(boundaryID)) {
              tmpData.d_SearchStatus = FoundOnBoundary;
              std::vector<AMP::Mesh::MeshElement> faceVertices = meshElementFaces[f].getElements(AMP::Mesh::Vertex);
              AMP_CHECK_ASSERT( faceVertices.size() == 4 );
              for (size_t v = 0; v < 4; ++v) {
                tmpData.d_FaceVerticesIDs[v] = faceVertices[v].globalID();
              } // end for v
              d_volume_elements[elementLocalID]->project_on_face(f, pointLocalCoords_ptr,
                  &(tmpData.d_ProjectionLocalCoordsOnFace[0]), &(tmpData.d_ShiftGlobalCoords[0]));
              break; // we assume only one face will be on the boundary
            } // end if
          } // end for f
        } else { // point was found but element is not on boundary
          tmpData.d_SearchStatus = FoundNotOnBoundary;
        } // end if
        sendData[d_sendDisps[pointOwnerRank] + tmpSendCnts[pointOwnerRank]] = tmpData;
        ++tmpSendCnts[pointOwnerRank];
      } //end i
      AMP_CHECK_ASSERT( std::equal(tmpSendCnts.begin(), tmpSendCnts.end(), d_sendCnts.begin()) );
      tmpSendCnts.clear();

      if(d_verbose) {
        comm.barrier();
        projectStep1Time = MPI_Wtime();
        if(!rank) {
          std::cout<<"Time for step-1 of project on boundary: "<<(projectStep1Time - projectBeginTime)<<" seconds."<<std::endl;
        }
      }

      std::vector<ProjectOnBoundaryData> recvData(d_recvDisps[npes-1] + d_recvCnts[npes-1]);

      comm.allToAll((!(sendData.empty()) ? &(sendData[0]) : NULL), &(d_sendCnts[0]), &(d_sendDisps[0]),
          (!(recvData.empty()) ? &(recvData[0]) : NULL), &(d_recvCnts[0]), &(d_recvDisps[0]), true);
      sendData.clear();

      faceVerticesGlobalIDs.resize(4*d_numLocalPts);
      std::fill(faceVerticesGlobalIDs.begin(), faceVerticesGlobalIDs.end(), AMP::Mesh::MeshElementID());

      projectionLocalCoordsOnFace.resize(2*d_numLocalPts);
      std::fill(projectionLocalCoordsOnFace.begin(), projectionLocalCoordsOnFace.end(), 0.0);

      shiftGlobalCoords.resize(3*d_numLocalPts);
      std::fill(shiftGlobalCoords.begin(), shiftGlobalCoords.end(), 0.0);

      flags.resize(d_numLocalPts);
      std::fill(flags.begin(), flags.end(), NotFound);

      for (int i = 0; i < npes; ++i) {
        for (int j = 0; j < d_recvCnts[i]; ++j) {
          const ProjectOnBoundaryData tmpData = recvData[d_recvDisps[i] + j];
          const size_t pointLocalID = tmpData.d_PointLocalID;
          if (tmpData.d_SearchStatus > flags[pointLocalID]) { // FoundOnBoundary overwrites FoundNotOnBoundary 
            flags[pointLocalID] = tmpData.d_SearchStatus;
            if (flags[pointLocalID] == FoundOnBoundary) {
              for (size_t d = 0; d < 2; ++d) {
                projectionLocalCoordsOnFace[2*pointLocalID+d] = tmpData.d_ProjectionLocalCoordsOnFace[d];
              } // end for d
              for (size_t d = 0; d < 3; ++d) {
                shiftGlobalCoords[3*pointLocalID+d] = tmpData.d_ShiftGlobalCoords[d];
              } // end for d
              for (size_t v = 0; v < 4; ++v) {
                faceVerticesGlobalIDs[4*pointLocalID+v] = tmpData.d_FaceVerticesIDs[v]; 
              } // end for v 
            } // end if
          } // end if
        } // end for j
      } // end for i

      recvData.clear();

      if(d_verbose) {
        comm.barrier();
      }
      projectStep2Time = MPI_Wtime();
      d_timingMeasurements[ProjectionOnBoundaryID] = projectStep2Time - projectBeginTime;
      if(d_verbose) {
        if(!rank) {
          std::cout<<"Time for step-2 of project on boundary: "<<(projectStep2Time - projectStep1Time)<<" seconds."<<std::endl;
        }
      }
    }

    void DendroSearch::setupDSforSearch() {
      if(d_meshAdapter == NULL) {
        return;
      }

      const unsigned int MaxDepth = 30;

      AMP::AMP_MPI meshComm = d_meshAdapter->getComm();
      const int rank = meshComm.getRank();
      const int npes = meshComm.getSize();

      double setupBeginTime, setupEndTime;
      if(d_verbose) {
        meshComm.barrier();
      }
      setupBeginTime = MPI_Wtime();

      std::vector<double> box = d_meshAdapter->getBoundingBox();
      for(int i = 0; i < d_meshAdapter->getDim(); ++i) {
        d_minCoords[i] = box[(2*i) + 0];
        double maxCoord = box[(2*i) + 1];
        d_scalingFactor[i] = 1.0/(1.0e-10 + maxCoord - d_minCoords[i]);
      }//end i

      size_t globalNumElems = d_meshAdapter->numGlobalElements(AMP::Mesh::Volume);
      if(d_verbose) {
        if(!rank) {
          d_oStream<<"Total number of mesh elements = "<<globalNumElems<<std::endl;
        }
      }

      double avgHboxInv = std::pow(globalNumElems, (1.0/3.0));
      AMP_CHECK_ASSERT(avgHboxInv > 1.0);
      d_boxLevel = binOp::fastLog2(static_cast<unsigned int>(std::ceil(avgHboxInv)));
      AMP_CHECK_ASSERT(d_boxLevel < MaxDepth);
      if(d_boxLevel > 5) {
        d_boxLevel = 5;
      }
      const double hBox = 1.0/(static_cast<double>(1u << d_boxLevel));

      if(d_verbose) {
        if(!rank) {
          d_oStream<<"BoxLevel = "<<d_boxLevel<<std::endl;
        }
      }

      unsigned int twoPowFactor = (1u << (MaxDepth - d_boxLevel));

      size_t localNumElems = d_meshAdapter->numLocalElements(AMP::Mesh::Volume);
      AMP_CHECK_ASSERT(localNumElems > 0);

      meshComm.barrier();
      double loop1TimeBegin = MPI_Wtime();

      std::vector<ot::TreeNode> tmpNodeList;
      std::vector<std::vector<int> > tmpElemIdList;

      d_volume_elements.clear();
      d_localElems.clear();
      d_volume_elements.reserve(localNumElems);
      d_localElems.reserve(localNumElems);
      AMP::Mesh::MeshIterator el = d_meshAdapter->getIterator(AMP::Mesh::Volume, 0);
      for(size_t eId = 0; eId < localNumElems; ++eId, ++el) {
        std::vector<int> eIdSingleton(1, eId);
        d_localElems.push_back(*el);
        std::vector<AMP::Mesh::MeshElement> vertices = el->getElements(AMP::Mesh::Vertex);
        double support_points[24];
        int minId[3];
        int maxId[3];
        for (size_t j = 0; j < vertices.size(); ++j) {
          std::vector<double> pt = vertices[j].coord();
          double scaledPt[3];
          for(int k = 0; k < 3; ++k) {
            support_points[(3*j) + k] = pt[k];
            scaledPt[k] = ((pt[k] - d_minCoords[k])*d_scalingFactor[k]);
            int id = static_cast<int>(scaledPt[k]/hBox);
            if(j == 0) {
              minId[k] = id;
              maxId[k] = id;
            } else {
              if(minId[k] > id) {
                minId[k] = id;
              }
              if(maxId[k] < id) {
                maxId[k] = id;
              }
            }
          }//end k
        }//end j
        d_volume_elements.push_back(new hex8_element_t(support_points));
        //PERFORMANCE IMPROVEMENT: We can skip the boxes that lie
        //completely outside the element.
        for(int k = minId[2]; k <= maxId[2]; ++k) {
          for(int j = minId[1]; j <= maxId[1]; ++j) {
            for(int i = minId[0]; i <= maxId[0]; ++i) {
              unsigned int bX = i*twoPowFactor;
              unsigned int bY = j*twoPowFactor;
              unsigned int bZ = k*twoPowFactor;
              ot::TreeNode box(bX, bY, bZ, d_boxLevel, 3, MaxDepth);
              unsigned int retIdx;
              bool found = seq::maxLowerBound<ot::TreeNode>(tmpNodeList, box, retIdx, NULL, NULL);
              if(found) {
                if(tmpNodeList[retIdx] == box) {
                  tmpElemIdList[retIdx].push_back(eId);
                } else {
                  tmpNodeList.insert(tmpNodeList.begin() + retIdx + 1, box);
                  tmpElemIdList.insert(tmpElemIdList.begin() + retIdx + 1, eIdSingleton);
                }
              } else {
                tmpNodeList.insert(tmpNodeList.begin(), box);
                tmpElemIdList.insert(tmpElemIdList.begin(), eIdSingleton);
              }
            }//end i 
          }//end j 
        }//end k 
      }//end eId

      meshComm.barrier();
      double loop1TimeEnd = MPI_Wtime();
      if(!rank) {
        d_oStream<<"SetupDS: Loop1-time = "<<(loop1TimeEnd - loop1TimeBegin)<<std::endl; 
      }

      d_nodeList.clear();
      d_rankList.clear();
      d_elemIdList.clear();

      if(npes == 1) {
        swap(d_nodeList, tmpNodeList);

        for(size_t i = 0; i < d_nodeList.size(); ++i) {
          d_nodeList[i].setWeight(tmpElemIdList[i].size());
          d_elemIdList.insert(d_elemIdList.end(), tmpElemIdList[i].begin(), tmpElemIdList[i].end());
          tmpElemIdList[i].clear();
        }//end i
        tmpElemIdList.clear();

        d_rankList.resize(d_elemIdList.size(), 0);
      } else {
        int numInitialLocalOcts = tmpNodeList.size();
        int numInitialGlobalOcts = meshComm.sumReduce<int>(numInitialLocalOcts);
        AMP_CHECK_ASSERT(numInitialGlobalOcts > 0);
        if(d_verbose) {
          if(!rank) {
            d_oStream<<"Total num initial octants = "<<numInitialGlobalOcts <<std::endl;
          }
        }

        if(numInitialGlobalOcts <= npes) {
          std::vector<ot::TreeNode> globalNodeList(numInitialGlobalOcts);

          ot::TreeNode* sendOctPtr = NULL;
          if(numInitialLocalOcts > 0) {
            sendOctPtr = &(tmpNodeList[0]);
          }

          meshComm.allGather<ot::TreeNode>(sendOctPtr, numInitialLocalOcts,
              &(globalNodeList[0]), NULL, NULL, false);

          seq::makeVectorUnique(globalNodeList, false);

          std::vector<int> sendEidList;
          std::vector<int> sendEidCnts(npes, 0);
          for(size_t i = 0, j = 0; i < numInitialLocalOcts; ++i, ++j) {
            while(tmpNodeList[i] < globalNodeList[j]) {
              ++j;
            }
            sendEidCnts[j] = tmpElemIdList[i].size();
            sendEidList.insert(sendEidList.end(), tmpElemIdList[i].begin(), tmpEledIdList[i].end());
            tmpElemIdList[i].clear();
          }//end i
          tmpElemIdList.clear();
          tmpNodeList.clear();

          std::vector<int> recvEidCnts(npes);
          meshComm.allToAll<int>(npes, &(sendEidCnts[0]), &(recvEidCnts[0]));

          std::vector<int> sendEidDisps(npes);
          std::vector<int> recvEidDisps(npes);
          sendEidDisps[0] = 0;
          recvEidDisps[0] = 0;
          for(int i = 1; i < npes; ++i) {
            sendEidDisps[i] = sendEidDisps[i - 1] + sendEidCnts[i - 1];
            recvEidDisps[i] = recvEidDisps[i - 1] + recvEidCnts[i - 1];
          }//end i

          std::vector<int> recvEidList(recvEidDisps[npes - 1] + recvEidCnts[npes - 1]);
          meshComm.allToAll<int>(&(sendEidList[0]), &(sendEidCnts[0]), &(sendEidDisps[0]), 
              &(recvEidList[0]), &(recvEidCnts[0]), &(recvEidDisps[0]), true);
          sendEidDisps.clear();
          sendEidCnts.clear();
          sendEidList.clear();

          if(rank < globalNodeList.size()) {
            d_nodeList.resize(1, (globalNodeList[rank]));
            d_nodeList[0].setWeight(recvEidList.size());
          }
          globalNodeList.clear();

          swap(d_elemIdList, recvEidList);

          d_rankList.resize(d_elemIdList.size());
          for(int i = 0; i < npes; ++i) {
            for(int j = 0; j < recvEidCnts[i]; ++j) {
              d_rankList[recvEidDisps[i] + j] = i;
            }//end j
          }//end i
          recvEidDisps.clear();
          recvEidCnts.clear();
        } else {
          int scanResult;
          meshComm.sumScan<int>(&numInitialLocalOcts, &scanResults, 1);
          int globalOffset = scanResult - numInitialLocalOcts;
          for(size_t i = 0; i < numInitialLocalOcts; ++i) {
            tmpNodeList[i].setWeight(globalOffset + i);
          }//end i

          //PERFORMANCE IMPROVEMENT: This parallel sort can be improved. We can
          //make use of the fact that tmpNodeList is already sorted and unique
          //on each processor. SampleSort modifies the input vector too. So, we
          //need to make a copy of tmpNodeList.
          std::vector<ot::TreeNode> sortInpVec = tmpNodeList;
          std::vector<ot::TreeNode> sortOutVec;
          par::sampleSort<ot::TreeNode>(sortInpVec, sortOutVec, meshComm.getCommunicator());
          sortInpVec.clear();

          std::vector<std::vector<int> > globalIndices;
          if(!(sortOutVec.empty())) {
            globalIndices.push_back(std::vector<int>(1, sortOutVec[0].getWeight()));
            sortOutVec[0].setWeight(0);
            sortInpVec.push_back(sortOutVec[0]);
          }
          for(size_t i = 1; i < sortOutVec.size(); ++i) {
            if(sortOutVec[i - 1] == sortOutVec[i]) {
              globalIndices[globalIndices.size() - 1].push_back(sortOutVec[i].getWeight());
            } else {
              globalIndices.push_back(std::vector<int>(1, sortOutVec[i].getWeight()));
              sortOutVec[i].setWeight(0);
              sortInpVec.push_back(sortOutVec[i]);
            }
          }//end i
          sortOutVec.clear();

          int localSortInpSz = sortInpVec.size();
          int globalSortInpSz = meshComm.sumReduce<int>(localSortInpSz);
          int newNpes = npes;
          int avgSortInpSz = globalSortInpSz/newNpes;
          if(newNpes > 2) {
            while(avgSortInpSz < 2) {
              --newNpes;
              avgSortInpSz = globalSortInpSz/newNpes;
              if(newNpes == 2) {
                break;
              }
            }
          }
          int extraSortInpSz = globalSortInpSz%newNpes;
          int newLocalSortInpSz = avgSortInpSz;
          if(rank < extraSortInpSz) {
            ++newLocalSortInpSz;
          }

        }
      }

      d_stIdxList.resize(d_nodeList.size());

      ot::TreeNode firstNode;
      if(!(d_nodeList.empty())) {
        firstNode = d_nodeList[0];
        firstNode.setWeight(rank);

        d_stIdxList[0] = 0;
        for(size_t i = 1; i < d_nodeList.size(); ++i) {
          d_stIdxList[i] = d_stIdxList[i - 1] + d_nodeList[i - 1].getWeight();
        }//end i
      }
      d_mins.resize(npes);
      meshComm.allGather(firstNode, &(d_mins[0]));

      std::vector<ot::TreeNode> tmpMins;
      for(int i = 0; i < npes; ++i) {
        if(d_mins[i].getDim() > 0) {
          tmpMins.push_back(d_mins[i]);
        }
      }//end i
      swap(d_mins, tmpMins);
      tmpMins.clear();

      if(d_verbose) {
        int numFinalLocalOcts = d_nodeList.size();
        int numFinalGlobalOcts = meshComm.sumReduce(numFinalLocalOcts);
        if(!rank) {
          d_oStream<<"Total num final octants = "<<numFinalGlobalOcts <<std::endl;
        }
      }

      if(d_verbose) {
        meshComm.barrier();
      }
      setupEndTime = MPI_Wtime();
      d_timingMeasurements[Setup] = setupEndTime - setupBeginTime;

      if(d_verbose) {
        if(!rank) {
          d_oStream<<"Finished setting up DS for search in "<<(setupEndTime - setupBeginTime)<<" seconds."<<std::endl;
        }
      }
    }

    void DendroSearch::search(AMP::AMP_MPI comm, const std::vector<double> & pts) {
      const int rank = comm.getRank();
      const int npes = comm.getSize();

      double coarseSearchBeginTime, coarseSearchEndTime, fineSearchBeginTime, fineSearchEndTime;
      coarseSearchBeginTime = MPI_Wtime();

      std::vector<int> rankMap(npes);

      int myRank = -1;
      if(d_meshAdapter != NULL) {
        AMP::AMP_MPI meshComm = d_meshAdapter->getComm();
        myRank = meshComm.getRank();
      }

      comm.allGather(myRank, &(rankMap[0]));

      std::vector<int> invRankMap(npes, -1);
      for(int i = 0; i < npes; ++i) {
        if(rankMap[i] >= 0) {
          invRankMap[rankMap[i]] = i;
        }
      }//end i

      std::vector<double> bcastBuff(7);
      if(myRank == 0) {
        bcastBuff[0] = static_cast<double>(d_mins.size());
        std::copy(d_minCoords.begin(), d_minCoords.end(), &(bcastBuff[1])); 
        std::copy(d_scalingFactor.begin(), d_scalingFactor.end(), &(bcastBuff[4])); 
      }

      comm.bcast(&(bcastBuff[0]), 7, invRankMap[0]);

      int minsSize = static_cast<int>(bcastBuff[0]);
      std::copy(&(bcastBuff[1]), &(bcastBuff[1])+3, d_minCoords.begin());
      std::copy(&(bcastBuff[4]), &(bcastBuff[4])+3, d_scalingFactor.begin());
      bcastBuff.clear();

      d_mins.resize(minsSize);
      comm.bcast(&(d_mins[0]), minsSize, invRankMap[0]);

      d_numLocalPts = (pts.size())/3;

      double searchBeginTime=0., searchStep1Time=0.,
             searchStep2Time=0., searchStep3Time=0.,
             searchStep4Time=0., searchStep5Time=0., searchStep6Time=0.;

      if(d_verbose) {
        comm.barrier();
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
        ptsWrapper[i].node.setWeight(rank);
        ptsWrapper[i].values[0] = x;
        ptsWrapper[i].values[1] = y;
        ptsWrapper[i].values[2] = z;
        ptsWrapper[i].values[3] = i;     
      }//end i

      if(d_verbose) {
        comm.barrier();
        searchStep1Time = MPI_Wtime();
        if(!rank) {
          std::cout<<"Time for step-1 of search: "<<(searchStep1Time - searchBeginTime)<<" seconds."<<std::endl;
        }
      }

      //PERFORMANCE IMPROVEMENT: Should PtsWrapper be sorted or not?
      //If PtsWrapper is sorted (even just a local sort), we can skip the
      //binary searches and use binning instead.  Binning is amortized constant
      //time and using binary searches would be logarithmic. This is just a matter
      //of constants since sorting is also logarithmic.

      d_sendCnts.resize(npes);
      d_recvCnts.resize(npes);
      d_sendDisps.resize(npes);
      d_recvDisps.resize(npes);
      std::fill(d_sendCnts.begin(), d_sendCnts.end(), 0);

      std::vector<int> part(d_numLocalPts, -1);

      for(int i = 0; i < d_numLocalPts; ++i) {
        unsigned int retIdx;
        bool found = seq::maxLowerBound<ot::TreeNode>(d_mins, (ptsWrapper[i].node), retIdx, NULL, NULL);
        if(found) {
          part[i] = invRankMap[d_mins[retIdx].getWeight()];
          ++(d_sendCnts[part[i]]);
        }
      }//end i

      if(d_verbose) {
        comm.barrier();
        searchStep2Time = MPI_Wtime();
        if(!rank) {
          std::cout<<"Time for step-2 of search: "<<(searchStep2Time - searchStep1Time)<<" seconds."<<std::endl;
        }
      }

      comm.allToAll(1, &(d_sendCnts[0]), &(d_recvCnts[0]));

      d_sendDisps[0] = 0;
      d_recvDisps[0] = 0;
      for(int i = 1; i < npes; ++i) {
        d_sendDisps[i] = d_sendDisps[i - 1] + d_sendCnts[i - 1];
        d_recvDisps[i] = d_recvDisps[i - 1] + d_recvCnts[i - 1];
      }//end i

      std::vector<ot::NodeAndValues<double, 4> > sendList(d_sendDisps[npes - 1] + d_sendCnts[npes - 1]);

      std::fill(d_sendCnts.begin(), d_sendCnts.end(), 0);

      for(int i = 0; i < d_numLocalPts; ++i) {
        if(part[i] >= 0) {
          sendList[d_sendDisps[part[i]] + d_sendCnts[part[i]]] = ptsWrapper[i];
          ++(d_sendCnts[part[i]]);
        }
      }//end i
      ptsWrapper.clear();

      std::vector<ot::NodeAndValues<double, 4> > recvList(d_recvDisps[npes - 1] + d_recvCnts[npes - 1]);
      comm.allToAll((!(sendList.empty()) ? &(sendList[0]) : NULL), &(d_sendCnts[0]), &(d_sendDisps[0]),
          (!(recvList.empty()) ? &(recvList[0]) : NULL), &(d_recvCnts[0]), &(d_recvDisps[0]), true);
      sendList.clear();

      if(d_verbose) {
        comm.barrier();
        searchStep3Time = MPI_Wtime();
        if(!rank) {
          std::cout<<"Time for step-3 of search: "<<(searchStep3Time - searchStep2Time)<<" seconds."<<std::endl;
        }
      }

      std::fill(d_sendCnts.begin(), d_sendCnts.end(), 0);

      std::vector<int> ptToOctMap((recvList.size()), -1);
      for(size_t i = 0; i < recvList.size(); ++i) {
        unsigned int retIdx;
        seq::maxLowerBound<ot::TreeNode>(d_nodeList, (recvList[i].node), retIdx, NULL, NULL);
        if( d_nodeList[retIdx].isAncestor(recvList[i].node) ) {
          ptToOctMap[i] = retIdx;
          int stIdx = d_stIdxList[retIdx];
          for(size_t j = 0; j < d_nodeList[retIdx].getWeight(); ++j) {
            d_sendCnts[d_rankList[stIdx + j]]++;
          }//end j
        }
      }//end i

      if(d_verbose) {
        comm.barrier();
        searchStep4Time = MPI_Wtime();
        if(!rank) {
          std::cout<<"Time for step-4 of search: "<<(searchStep4Time - searchStep3Time)<<" seconds."<<std::endl;
        }
      }

      comm.allToAll(1, &(d_sendCnts[0]), &(d_recvCnts[0]));

      d_sendDisps[0] = 0;
      d_recvDisps[0] = 0;
      for(int i = 1; i < npes; ++i) {
        d_sendDisps[i] = d_sendDisps[i - 1] + d_sendCnts[i - 1];
        d_recvDisps[i] = d_recvDisps[i - 1] + d_recvCnts[i - 1];
      }//end i

      std::vector<double> sendPtsList(6*(d_sendDisps[npes - 1] + d_sendCnts[npes - 1]));

      std::fill(d_sendCnts.begin(), d_sendCnts.end(), 0);

      for(size_t i = 0; i < ptToOctMap.size(); ++i) {
        if(ptToOctMap[i] >= 0) {
          int stIdx = d_stIdxList[ptToOctMap[i]];
          for(size_t j = 0; j < d_nodeList[ptToOctMap[i]].getWeight(); ++j) {
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

      for(int i = 0; i < npes; ++i) {
        d_sendCnts[i] *= 6;
        d_sendDisps[i] *= 6;
        d_recvCnts[i] *= 6;
        d_recvDisps[i] *= 6;
      }//end i

      std::vector<double> recvPtsList(d_recvDisps[npes - 1] + d_recvCnts[npes - 1]);
      comm.allToAll((!(sendPtsList.empty()) ? &(sendPtsList[0]) : NULL), &(d_sendCnts[0]), &(d_sendDisps[0]), 
          (!(recvPtsList.empty()) ? &(recvPtsList[0]) : NULL), &(d_recvCnts[0]), &(d_recvDisps[0]), true);
      sendPtsList.clear();

      if(d_verbose) {
        comm.barrier();
        searchStep5Time = MPI_Wtime();
        if(!rank) {
          std::cout<<"Time for step-5 of search: "<<(searchStep5Time - searchStep4Time)<<" seconds."<<std::endl;
        }
      }

      coarseSearchEndTime = MPI_Wtime();
      d_timingMeasurements[CoarseSearch] = coarseSearchEndTime - coarseSearchBeginTime;
      fineSearchBeginTime = MPI_Wtime();

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
        if (d_volume_elements[eId]->within_bounding_box(tmpPtGlobalCoordPtr, d_tolerance)) {
          if (d_volume_elements[eId]->within_bounding_polyhedron(tmpPtGlobalCoordPtr, d_tolerance)) {
            d_volume_elements[eId]->map_global_to_local(tmpPtGlobalCoordPtr, &(tmpPtLocalCoord[0]));
            if (d_volume_elements[eId]->contains_point(&(tmpPtLocalCoord[0]), coordinates_are_local, d_tolerance)) {
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

      comm.allToAll(1, &(d_sendCnts[0]), &(d_recvCnts[0]));

      d_sendDisps[0] = 0;
      d_recvDisps[0] = 0;
      for(int i = 1; i < npes; ++i) {
        d_sendDisps[i] = d_sendDisps[i - 1] + d_sendCnts[i - 1];
        d_recvDisps[i] = d_recvDisps[i - 1] + d_recvCnts[i - 1];
      }//end i

      if(d_verbose) {
        comm.barrier();
        searchStep6Time = MPI_Wtime();
        if(!rank) {
          std::cout<<"Time for step-6 of search: "<<(searchStep6Time - searchStep5Time)<<" seconds."<<std::endl;
        }
      }

      fineSearchEndTime = MPI_Wtime();
      d_timingMeasurements[FineSearch] = fineSearchEndTime - fineSearchBeginTime;
    }

    void DendroSearch::interpolate(AMP::AMP_MPI comm, AMP::LinearAlgebra::Vector::const_shared_ptr vectorField, const unsigned int dofsPerNode,
        std::vector<double> & results, std::vector<bool> & foundPt) {
      const int rank = comm.getRank();
      const int npes = comm.getSize();

      double interpolateBeginTime, interpolateStep1Time=0., interpolateStep2Time;
      if(d_verbose) {
        comm.barrier();
      }
      interpolateBeginTime = MPI_Wtime();

      AMP_CHECK_ASSERT( vectorField->getUpdateStatus()==AMP::LinearAlgebra::Vector::UNCHANGED );
      AMP::Discretization::DOFManager::shared_ptr dofManager = vectorField->getDOFManager();

      for(int i = 0; i < npes; ++i) {
        d_sendCnts[i] *= (dofsPerNode + 1);
        d_recvCnts[i] *= (dofsPerNode + 1);
        d_sendDisps[i] *= (dofsPerNode + 1);
        d_recvDisps[i] *= (dofsPerNode + 1);
      }// end for i

      std::vector<double> sendResults(d_sendDisps[npes - 1] + d_sendCnts[npes - 1]);

      std::vector<int> tmpSendCnts(npes, 0);

      std::vector<double> basis_functions_values(8);
      for(size_t i = 0; i < d_foundPts.size(); i += 6) {
        unsigned int elemLocalId = static_cast<unsigned int>(d_foundPts[i]);
        std::vector<AMP::Mesh::MeshElement> amp_vector_support_points =
          d_localElems[elemLocalId].getElements(AMP::Mesh::Vertex);
        hex8_element_t::get_basis_functions_values(&(d_foundPts[i + 1]), &(basis_functions_values[0]));

        std::vector<double> value(dofsPerNode, 0.0);
        for (unsigned int j = 0; j < 8; ++j) {
          std::vector<size_t> globalID;
          dofManager->getDOFs(amp_vector_support_points[j].globalID(), globalID);
          AMP_CHECK_ASSERT(globalID.size() == dofsPerNode);
          for(size_t d = 0; d < dofsPerNode; ++d) {
            double vecVal = vectorField->getValueByGlobalID(globalID[d]);
            value[d] += (vecVal * basis_functions_values[j]);
          }//end d
        } // end j
        unsigned int ptProcId = static_cast<unsigned int>(d_foundPts[i + 5]);
        sendResults[d_sendDisps[ptProcId] + tmpSendCnts[ptProcId]] = d_foundPts[i + 4];
        ++(tmpSendCnts[ptProcId]);
        for(size_t d = 0; d < dofsPerNode; ++d) {
          sendResults[d_sendDisps[ptProcId] + tmpSendCnts[ptProcId]] = value[d];
          ++(tmpSendCnts[ptProcId]);
        }//end d   
      }//end i
      tmpSendCnts.clear();

      if(d_verbose) {
        comm.barrier();
        interpolateStep1Time = MPI_Wtime();
        if(!rank) {
          std::cout<<"Time for step-1 of interpolate: "<<(interpolateStep1Time - interpolateBeginTime)<<" seconds."<<std::endl;
        }
      }

      std::vector<double> recvResults(d_recvDisps[npes - 1] + d_recvCnts[npes - 1]);

      comm.allToAll((!(sendResults.empty()) ? &(sendResults[0]) : NULL), &(d_sendCnts[0]), &(d_sendDisps[0]),
          (!(recvResults.empty()) ? &(recvResults[0]) : NULL), &(d_recvCnts[0]), &(d_recvDisps[0]), true);
      sendResults.clear();

      //Points that are not found will have a result = 0. 
      results.resize(dofsPerNode*d_numLocalPts);
      for (unsigned int i = 0; i < results.size(); ++i) {
        results[i] = 0.0;
      } // end for i

      foundPt.resize(d_numLocalPts);
      for (unsigned int i = 0; i < foundPt.size(); ++i) {
        foundPt[i] = static_cast<bool>(NotFound);
      } // end for i

      for(size_t i = 0; i < recvResults.size(); i += (dofsPerNode + 1)) {
        unsigned int locId = static_cast<unsigned int>(recvResults[i]);
        foundPt[locId] = static_cast<bool>(Found);
        for(size_t d = 0; d < dofsPerNode; ++d) {
          results[(locId*dofsPerNode) + d] = recvResults[i + d + 1];
        }//end d
      }//end i

      for(int i = 0; i < npes; ++i) {
        d_sendCnts[i] /= (dofsPerNode + 1);
        d_recvCnts[i] /= (dofsPerNode + 1);
        d_sendDisps[i] /= (dofsPerNode + 1);
        d_recvDisps[i] /= (dofsPerNode + 1);
      } // end for i

      if(d_verbose) {
        comm.barrier();
      }
      interpolateStep2Time = MPI_Wtime();
      d_timingMeasurements[Interpolation] = interpolateStep2Time - interpolateBeginTime;
      if(d_verbose) {
        if(!rank) {
          std::cout<<"Time for step-2 of interpolate: "<<(interpolateStep2Time - interpolateStep1Time)<<" seconds."<<std::endl;
        }
      }
    }

    void DendroSearch::reportTiming(size_t n, TimingType const * timingTypes, double * timingMeasurements) {
      AMP_INSIST(!d_verbose, "verbose mode in DendroSearch implies calls to MPI_Barrier so timing measurements are bad!");
      for (size_t i = 0; i < n; ++i) {
        timingMeasurements[i] = d_timingMeasurements[timingTypes[i]];
      } // end for i
    }

  }
}



