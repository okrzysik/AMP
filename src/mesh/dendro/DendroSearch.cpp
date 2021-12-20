
#include "AMP/mesh/dendro/DendroSearch.h"

#include "AMP/mesh/latex_visualization_tools.h"
#include "AMP/utils/Utilities.h"
#include "externVars.h"
#include <fstream>
#include <iomanip>
#include <numeric>
#include <set>

namespace AMP {
namespace Mesh {


DendroSearch::DendroSearch( AMP::Mesh::Mesh::shared_ptr mesh, bool verbose, std::ostream &oStream )
    : d_meshAdapter( mesh ),
      d_verbose( verbose ),
      d_oStream( oStream ),
      d_timingMeasurements( std::vector<double>( numTimingTypes, -1.0 ) ),
      d_tolerance( 1.0e-12 )
{
    d_minCoords.resize( 3 );
    d_scalingFactor.resize( 3 );
    setupDSforSearch();
}


void DendroSearch::searchAndInterpolate( AMP::AMP_MPI comm,
                                         AMP::LinearAlgebra::Vector::const_shared_ptr vectorField,
                                         const unsigned int dofsPerNode,
                                         const std::vector<double> &pts,
                                         std::vector<double> &results,
                                         std::vector<bool> &foundPt )
{
    search( comm, pts );
    interpolate( comm, vectorField, dofsPerNode, results, foundPt );
}


void DendroSearch::setTolerance( double tolerance ) { d_tolerance = tolerance; }


void DendroSearch::projectOnBoundaryID(
    AMP::AMP_MPI comm,
    const int boundaryID,
    std::vector<AMP::Mesh::MeshElementID> &faceVerticesGlobalIDs,
    std::vector<double> &shiftGlobalCoords,
    std::vector<double> &projectionLocalCoordsOnGeomType::Face,
    std::vector<int> &flags )
{
    std::vector<size_t> faceLocalIndices;
    std::vector<AMP::Mesh::MeshElementID> volumeGlobalIDs;
    std::cerr << "***WARNING*** the version of DendroSearch::projectOnBoundaryID() you are using "
                 "is deprecated!"
              << std::endl;
    projectOnBoundaryID( comm,
                         boundaryID,
                         faceVerticesGlobalIDs,
                         shiftGlobalCoords,
                         projectionLocalCoordsOnGeomType::Face,
                         flags,
                         volumeGlobalIDs,
                         faceLocalIndices );
}


void DendroSearch::projectOnBoundaryID(
    AMP::AMP_MPI comm,
    const int boundaryID,
    std::vector<AMP::Mesh::MeshElementID> &faceVerticesGlobalIDs,
    std::vector<double> &shiftGlobalCoords,
    std::vector<double> &projectionLocalCoordsOnGeomType::Face,
    std::vector<int> &flags,
    std::vector<AMP::Mesh::MeshElementID> &volumeGlobalIDs,
    std::vector<size_t> &faceLocalIndices )
{
    const int rank = comm.getRank();
    const int npes = comm.getSize();

    double projectBeginTime, projectStep1Time = 0., projectStep2Time;
    if ( d_verbose ) {
        comm.barrier();
    }
    projectBeginTime = MPI_Wtime();

    std::vector<ProjectOnBoundaryData> sendData( d_sendDisps[npes - 1] + d_sendCnts[npes - 1] );

    std::vector<int> tmpSendCnts( npes, 0 );

    unsigned int const *faceOrdering = hex8_element_t::get_faces();
    std::vector<size_t> mapGeomType::Faces( 6, 6 );
    if ( !d_foundPts.empty() ) {
        AMP::Mesh::MeshElement controlGeomType::VolumeElement = d_localElems[0];
        std::vector<AMP::Mesh::MeshElement> controlGeomType::VolumeElementVertices =
            controlGeomType::VolumeElement.getElements( AMP::Mesh::GeomType::Vertex );
        AMP_CHECK_ASSERT( controlGeomType::VolumeElementVertices.size() == 8 );
        AMP::Mesh::MeshElement hex8ElementGeomType::FaceVertices[24];
        for ( size_t f = 0; f < 6; ++f ) {
            for ( size_t v = 0; v < 4; ++v ) {
                hex8ElementGeomType::FaceVertices[4 * f + v] =
                    controlGeomType::VolumeElementVertices[faceOrdering[4 * f + v]];
            } // end for v
            std::sort( &( hex8ElementGeomType::FaceVertices[4 * f] ),
                       &( hex8ElementGeomType::FaceVertices[4 * f] ) + 4 );
        } // end for f
        std::vector<AMP::Mesh::MeshElement> controlGeomType::VolumeElementGeomType::Faces =
            controlGeomType::VolumeElement.getElements( AMP::Mesh::GeomType::Face );
        AMP_CHECK_ASSERT( controlGeomType::VolumeElementGeomType::Faces.size() == 6 );
        //        std::vector<size_t> mapGeomType::Faces(6, 6);
        for ( size_t f = 0; f < 6; ++f ) {
            std::vector<AMP::Mesh::MeshElement> faceVertices =
                controlGeomType::VolumeElementGeomType::Faces[f].getElements(
                    AMP::Mesh::GeomType::Vertex );
            AMP_CHECK_ASSERT( faceVertices.size() == 4 );
            std::sort( faceVertices.begin(), faceVertices.end() );
            for ( size_t g = 0; g < 6; ++g ) {
                if ( std::equal( faceVertices.begin(),
                                 faceVertices.end(),
                                 &( hex8ElementGeomType::FaceVertices[4 * g] ) ) ) {
                    mapGeomType::Faces[f] = g;
                    //              std::cout<<f<<"  ->  "<<g<<"\n";
                    break;
                } // end if
            }     // end for g
        }         // end for f
    }             // end if

    for ( unsigned int i = 0; i < d_foundPts.size(); i += 6 ) {
        ProjectOnBoundaryData tmpData;
        const double *pointLocalCoords_ptr = &( d_foundPts[i + 1] );
        const size_t pointLocalID          = static_cast<size_t>( d_foundPts[i + 4] );
        const size_t pointOwnerRank        = static_cast<size_t>( d_foundPts[i + 5] );
        const size_t elementLocalID        = static_cast<size_t>( d_foundPts[i] );
        tmpData.d_PointLocalID             = pointLocalID;
        if ( d_localElems[elementLocalID].isOnBoundary(
                 boundaryID ) ) { // point was found and element is on boundary
            std::vector<AMP::Mesh::MeshElement> meshElementGeomType::Faces =
                d_localElems[elementLocalID].getElements( AMP::Mesh::GeomType::Face );
            AMP_CHECK_ASSERT( meshElementGeomType::Faces.size() == 6 );
            for ( size_t f = 0; f < 6; ++f ) {
                if ( meshElementGeomType::Faces[f].isOnBoundary( boundaryID ) ) {
                    tmpData.d_SearchStatus = FoundOnBoundary;
                    //              tmpData.d_GeomType::FaceLocalIndex = f;
                    tmpData.d_GeomType::FaceLocalIndex = mapGeomType::Faces[f];
                    tmpData.d_GeomType::VolumeID       = d_localElems[elementLocalID].globalID();
                    //              std::vector<AMP::Mesh::MeshElement> faceVertices =
                    //              meshElementGeomType::Faces[f].getElements(AMP::Mesh::GeomType::Vertex);
                    //              AMP_CHECK_ASSERT( faceVertices.size() == 4 );
                    std::vector<AMP::Mesh::MeshElement> meshElementVertices =
                        d_localElems[elementLocalID].getElements( AMP::Mesh::GeomType::Vertex );
                    AMP_CHECK_ASSERT( meshElementVertices.size() == 8 );
                    for ( size_t v = 0; v < 4; ++v ) {
                        //                tmpData.d_GeomType::FaceVerticesIDs[v] =
                        //                faceVertices[v].globalID();
                        tmpData.d_GeomType::FaceVerticesIDs[v] =
                            meshElementVertices[faceOrdering[4 * mapGeomType::Faces[f] + v]]
                                .globalID();
                    } // end for v
                      //              d_volume_elements[elementLocalID]->project_on_face(f,
                      //              pointLocalCoords_ptr,
                    d_volume_elements[elementLocalID]->project_on_face(
                        mapGeomType::Faces[f],
                        pointLocalCoords_ptr,
                        &( tmpData.d_ProjectionLocalCoordsOnGeomType::Face[0] ),
                        &( tmpData.d_ShiftGlobalCoords[0] ) );
                    break; // we assume only one face will be on the boundary
                }          // end if
            }              // end for f
        } else {           // point was found but element is not on boundary
            tmpData.d_SearchStatus       = FoundNotOnBoundary;
            tmpData.d_GeomType::VolumeID = d_localElems[elementLocalID].globalID();
        } // end if
        sendData[d_sendDisps[pointOwnerRank] + tmpSendCnts[pointOwnerRank]] = tmpData;
        ++tmpSendCnts[pointOwnerRank];
    } // end i
    AMP_CHECK_ASSERT( std::equal( tmpSendCnts.begin(), tmpSendCnts.end(), d_sendCnts.begin() ) );
    tmpSendCnts.clear();

    if ( d_verbose ) {
        comm.barrier();
        projectStep1Time = MPI_Wtime();
        if ( !rank ) {
            std::cout << "Time for step-1 of project on boundary: "
                      << ( projectStep1Time - projectBeginTime ) << " seconds." << std::endl;
        }
    }

    std::vector<ProjectOnBoundaryData> recvData( d_recvDisps[npes - 1] + d_recvCnts[npes - 1] );

    comm.allToAll<ProjectOnBoundaryData>( ( !( sendData.empty() ) ? &( sendData[0] ) : NULL ),
                                          &( d_sendCnts[0] ),
                                          &( d_sendDisps[0] ),
                                          ( !( recvData.empty() ) ? &( recvData[0] ) : NULL ),
                                          &( d_recvCnts[0] ),
                                          &( d_recvDisps[0] ),
                                          true );
    sendData.clear();

    faceVerticesGlobalIDs.resize( 4 * d_numLocalPts );
    std::fill(
        faceVerticesGlobalIDs.begin(), faceVerticesGlobalIDs.end(), AMP::Mesh::MeshElementID() );

    projectionLocalCoordsOnGeomType::Face.resize( 2 * d_numLocalPts );
    std::fill( projectionLocalCoordsOnGeomType::Face.begin(),
               projectionLocalCoordsOnGeomType::Face.end(),
               0.0 );

    shiftGlobalCoords.resize( 3 * d_numLocalPts );
    std::fill( shiftGlobalCoords.begin(), shiftGlobalCoords.end(), 0.0 );

    flags.resize( d_numLocalPts );
    std::fill( flags.begin(), flags.end(), NotFound );

    faceLocalIndices.resize( d_numLocalPts );
    std::fill( faceLocalIndices.begin(), faceLocalIndices.end(), 7 );

    volumeGlobalIDs.resize( d_numLocalPts );
    std::fill( volumeGlobalIDs.begin(), volumeGlobalIDs.end(), AMP::Mesh::MeshElementID() );

    for ( int i = 0; i < npes; ++i ) {
        for ( int j = 0; j < d_recvCnts[i]; ++j ) {
            const ProjectOnBoundaryData tmpData = recvData[d_recvDisps[i] + j];
            const size_t pointLocalID           = tmpData.d_PointLocalID;
            if ( tmpData.d_SearchStatus >
                 flags[pointLocalID] ) { // FoundOnBoundary overwrites FoundNotOnBoundary
                flags[pointLocalID]           = tmpData.d_SearchStatus;
                volumeGlobalIDs[pointLocalID] = tmpData.d_GeomType::VolumeID;
                if ( flags[pointLocalID] == FoundOnBoundary ) {
                    for ( size_t d = 0; d < 2; ++d ) {
                        projectionLocalCoordsOnGeomType::Face[2 * pointLocalID + d] =
                            tmpData.d_ProjectionLocalCoordsOnGeomType::Face[d];
                    } // end for d
                    for ( size_t d = 0; d < 3; ++d ) {
                        shiftGlobalCoords[3 * pointLocalID + d] = tmpData.d_ShiftGlobalCoords[d];
                    } // end for d
                    for ( size_t v = 0; v < 4; ++v ) {
                        faceVerticesGlobalIDs[4 * pointLocalID + v] =
                            tmpData.d_GeomType::FaceVerticesIDs[v];
                    } // end for v
                    faceLocalIndices[pointLocalID] = tmpData.d_GeomType::FaceLocalIndex;
                } // end if
            }     // end if
        }         // end for j
    }             // end for i

    recvData.clear();

    if ( d_verbose ) {
        comm.barrier();
    }
    projectStep2Time                             = MPI_Wtime();
    d_timingMeasurements[ProjectionOnBoundaryID] = projectStep2Time - projectBeginTime;
    if ( d_verbose ) {
        if ( !rank ) {
            std::cout << "Time for step-2 of project on boundary: "
                      << ( projectStep2Time - projectStep1Time ) << " seconds." << std::endl;
        }
    }
}


void DendroSearch::setupDSforSearch()
{
    if ( d_meshAdapter == NULL ) {
        return;
    }

    const unsigned int MaxDepth = 30;

    AMP::AMP_MPI meshComm = d_meshAdapter->getComm();
    const int rank        = meshComm.getRank();
    const int npes        = meshComm.getSize();

    double setupBeginTime, setupEndTime;
    if ( d_verbose ) {
        meshComm.barrier();
    }
    setupBeginTime = MPI_Wtime();

    std::vector<double> box = d_meshAdapter->getBoundingBox();
    for ( int i = 0; i < d_meshAdapter->getDim(); ++i ) {
        d_minCoords[i]     = box[( 2 * i ) + 0];
        double maxCoord    = box[( 2 * i ) + 1];
        d_scalingFactor[i] = 1.0 / ( 1.0e-10 + maxCoord - d_minCoords[i] );
    } // end i

    size_t globalNumElems = d_meshAdapter->numGlobalElements( AMP::Mesh::GeomType::Volume );
    if ( d_verbose ) {
        meshComm.barrier();
        if ( !rank ) {
            d_oStream << "Total number of mesh elements = " << globalNumElems << std::endl;
        }
    }

    double avgHboxInv = std::pow( globalNumElems, ( 1.0 / 3.0 ) );
    AMP_CHECK_ASSERT( avgHboxInv > 1.0 );
    d_boxLevel = binOp::fastLog2( static_cast<unsigned int>( std::ceil( avgHboxInv ) ) );
    if ( d_boxLevel > 5 ) {
        d_boxLevel = 5;
    }
    AMP_CHECK_ASSERT( d_boxLevel < MaxDepth );
    const double hBox = 1.0 / ( static_cast<double>( 1u << d_boxLevel ) );

    if ( d_verbose ) {
        meshComm.barrier();
        if ( !rank ) {
            d_oStream << "BoxLevel = " << d_boxLevel << std::endl;
        }
    }

    unsigned int twoPowFactor = ( 1u << ( MaxDepth - d_boxLevel ) );

    size_t localNumElems = d_meshAdapter->numLocalElements( AMP::Mesh::GeomType::Volume );
    AMP_CHECK_ASSERT( localNumElems > 0 );

    std::vector<ot::TreeNode> tmpNodeList;
    std::vector<std::vector<int>> tmpElemIdList;

    d_volume_elements.clear();
    d_localElems.clear();
    d_volume_elements.reserve( localNumElems );
    d_localElems.reserve( localNumElems );
    AMP::Mesh::MeshIterator el = d_meshAdapter->getIterator( AMP::Mesh::GeomType::Volume, 0 );
    for ( size_t eId = 0; eId < localNumElems; ++eId, ++el ) {
        std::vector<int> eIdSingleton( 1, eId );
        d_localElems.push_back( *el );
        std::vector<AMP::Mesh::MeshElement> vertices =
            el->getElements( AMP::Mesh::GeomType::Vertex );
        double support_points[24];
        std::vector<int> minId( 3, 0 );
        std::vector<int> maxId( 3, 0 );
        for ( size_t j = 0; j < vertices.size(); ++j ) {
            std::vector<double> pt = vertices[j].coord();
            double scaledPt[3];
            for ( int k = 0; k < 3; ++k ) {
                support_points[( 3 * j ) + k] = pt[k];
                scaledPt[k]                   = ( ( pt[k] - d_minCoords[k] ) * d_scalingFactor[k] );
                int id                        = static_cast<int>( scaledPt[k] / hBox );
                if ( j == 0 ) {
                    minId[k] = id;
                    maxId[k] = id;
                } else {
                    if ( minId[k] > id ) {
                        minId[k] = id;
                    }
                    if ( maxId[k] < id ) {
                        maxId[k] = id;
                    }
                }
            } // end k
        }     // end j
        d_volume_elements.push_back( new hex8_element_t( support_points ) );
        // PERFORMANCE IMPROVEMENT: We can skip the boxes that lie
        // completely outside the element.
        for ( int k = minId[2]; k <= maxId[2]; ++k ) {
            for ( int j = minId[1]; j <= maxId[1]; ++j ) {
                for ( int i = minId[0]; i <= maxId[0]; ++i ) {
                    unsigned int bX = i * twoPowFactor;
                    unsigned int bY = j * twoPowFactor;
                    unsigned int bZ = k * twoPowFactor;
                    ot::TreeNode box( bX, bY, bZ, d_boxLevel, 3, MaxDepth );
                    unsigned int retIdx = 0;
                    bool found =
                        seq::maxLowerBound<ot::TreeNode>( tmpNodeList, box, retIdx, NULL, NULL );
                    if ( found ) {
                        if ( tmpNodeList[retIdx] == box ) {
                            tmpElemIdList[retIdx].push_back( eId );
                        } else {
                            tmpNodeList.insert( tmpNodeList.begin() + retIdx + 1, box );
                            tmpElemIdList.insert( tmpElemIdList.begin() + retIdx + 1,
                                                  eIdSingleton );
                        }
                    } else {
                        tmpNodeList.insert( tmpNodeList.begin(), box );
                        tmpElemIdList.insert( tmpElemIdList.begin(), eIdSingleton );
                    }
                } // end i
            }     // end j
        }         // end k
    }             // end eId

    d_nodeList.clear();
    d_rankList.clear();
    d_elemIdList.clear();
    d_mins.clear();

    if ( npes == 1 ) {
        if ( d_verbose ) {
            meshComm.barrier();
            if ( !rank ) {
                d_oStream << "DendroSearch - Setup: Case A." << std::endl;
            }
        }

        swap( d_nodeList, tmpNodeList );

        for ( size_t i = 0; i < d_nodeList.size(); ++i ) {
            d_nodeList[i].setWeight( tmpElemIdList[i].size() );
            d_elemIdList.insert(
                d_elemIdList.end(), tmpElemIdList[i].begin(), tmpElemIdList[i].end() );
            tmpElemIdList[i].clear();
        } // end i
        tmpElemIdList.clear();

        d_rankList.resize( d_elemIdList.size(), 0 );

        AMP_CHECK_ASSERT( !( d_nodeList.empty() ) );

        ot::TreeNode firstNode = d_nodeList[0];
        firstNode.setWeight( rank );
        d_mins.resize( 1, firstNode );
    } else {
        int numInitialLocalOcts  = tmpNodeList.size();
        int numInitialGlobalOcts = meshComm.sumReduce<int>( numInitialLocalOcts );
        AMP_CHECK_ASSERT( numInitialGlobalOcts > 0 );
        if ( numInitialGlobalOcts <= npes ) {
            if ( d_verbose ) {
                meshComm.barrier();
                if ( !rank ) {
                    d_oStream << "DendroSearch - Setup: Case B. Total Number of Initial Octants = "
                              << numInitialGlobalOcts << std::endl;
                }
            }

            std::vector<ot::TreeNode> globalNodeList( numInitialGlobalOcts );

            ot::TreeNode *sendOctPtr = NULL;
            if ( numInitialLocalOcts > 0 ) {
                sendOctPtr = &( tmpNodeList[0] );
            }
            meshComm.allGather<ot::TreeNode>(
                sendOctPtr, numInitialLocalOcts, &( globalNodeList[0] ), NULL, NULL, false );

            seq::makeVectorUnique( globalNodeList, false );

            std::vector<int> sendEidList;
            for ( int i = 0; i < numInitialLocalOcts; ++i ) {
                sendEidList.insert(
                    sendEidList.end(), tmpElemIdList[i].begin(), tmpElemIdList[i].end() );
            } // end i

            std::vector<int> sendEidCnts( npes, 0 );
            for ( int i = 0; i < numInitialLocalOcts; ++i ) {
                unsigned int retIdx = 0;
                bool found          = seq::maxLowerBound<ot::TreeNode>(
                    globalNodeList, tmpNodeList[i], retIdx, NULL, NULL );
                AMP_CHECK_ASSERT( found );
                AMP_CHECK_ASSERT( globalNodeList[retIdx] == tmpNodeList[i] );
                sendEidCnts[retIdx] = tmpElemIdList[i].size();
            } // end i
            tmpElemIdList.clear();
            tmpNodeList.clear();

            std::vector<int> recvEidCnts( npes );
            meshComm.allToAll<int>( 1, &( sendEidCnts[0] ), &( recvEidCnts[0] ) );

            std::vector<int> sendEidDisps( npes );
            std::vector<int> recvEidDisps( npes );
            sendEidDisps[0] = 0;
            recvEidDisps[0] = 0;
            for ( int i = 1; i < npes; ++i ) {
                sendEidDisps[i] = sendEidDisps[i - 1] + sendEidCnts[i - 1];
                recvEidDisps[i] = recvEidDisps[i - 1] + recvEidCnts[i - 1];
            } // end i

            std::vector<int> recvEidList( recvEidDisps[npes - 1] + recvEidCnts[npes - 1] );
            int *sendEidListPtr = NULL;
            if ( !( sendEidList.empty() ) ) {
                sendEidListPtr = &( sendEidList[0] );
            }
            int *recvEidListPtr = NULL;
            if ( !( recvEidList.empty() ) ) {
                recvEidListPtr = &( recvEidList[0] );
            }
            meshComm.allToAll<int>( sendEidListPtr,
                                    &( sendEidCnts[0] ),
                                    &( sendEidDisps[0] ),
                                    recvEidListPtr,
                                    &( recvEidCnts[0] ),
                                    &( recvEidDisps[0] ),
                                    true );
            sendEidDisps.clear();
            sendEidCnts.clear();
            sendEidList.clear();

            if ( rank < static_cast<int>( globalNodeList.size() ) ) {
                d_nodeList.resize( 1, ( globalNodeList[rank] ) );
                d_nodeList[0].setWeight( recvEidList.size() );
            }

            swap( d_mins, globalNodeList );
            for ( size_t i = 0; i < d_mins.size(); ++i ) {
                d_mins[i].setWeight( i );
            } // end i

            swap( d_elemIdList, recvEidList );

            d_rankList.resize( d_elemIdList.size() );
            for ( int i = 0; i < npes; ++i ) {
                for ( int j = 0; j < recvEidCnts[i]; ++j ) {
                    d_rankList[recvEidDisps[i] + j] = i;
                } // end j
            }     // end i
            recvEidDisps.clear();
            recvEidCnts.clear();
        } else {
            if ( d_verbose ) {
                meshComm.barrier();
                if ( !rank ) {
                    d_oStream << "DendroSearch - Setup: Case C. Total Number of Initial Octants = "
                              << numInitialGlobalOcts << std::endl;
                }
            }

            // PERFORMANCE IMPROVEMENT: This parallel sort + unique step can be improved. We can
            // make use of the fact that tmpNodeList is already sorted and unique on each processor.
            d_nodeList = tmpNodeList;
            par::removeDuplicates<ot::TreeNode>( d_nodeList, false, meshComm.getCommunicator() );

            ot::TreeNode firstNode;
            if ( !( d_nodeList.empty() ) ) {
                firstNode = d_nodeList[0];
                firstNode.setWeight( rank );
            }
            d_mins.resize( npes );
            meshComm.allGather( firstNode, &( d_mins[0] ) );

            std::vector<ot::TreeNode> tmpMins;
            for ( int i = 0; i < npes; ++i ) {
                if ( d_mins[i].getDim() > 0 ) {
                    tmpMins.push_back( d_mins[i] );
                }
            } // end i
            swap( d_mins, tmpMins );
            tmpMins.clear();

            for ( int i = 0; i < numInitialLocalOcts; ++i ) {
                tmpNodeList[i].setWeight( tmpElemIdList[i].size() );
            } // end i

            std::vector<int> sendEidList;
            for ( int i = 0; i < numInitialLocalOcts; ++i ) {
                sendEidList.insert(
                    sendEidList.end(), tmpElemIdList[i].begin(), tmpElemIdList[i].end() );
                tmpElemIdList[i].clear();
            } // end i
            tmpElemIdList.clear();

            std::vector<int> sendOctCnts( npes, 0 );
            for ( int i = 0; i < numInitialLocalOcts; ++i ) {
                unsigned int retIdx = 0;
                bool found =
                    seq::maxLowerBound<ot::TreeNode>( d_mins, tmpNodeList[i], retIdx, NULL, NULL );
                AMP_CHECK_ASSERT( found );
                ++( sendOctCnts[d_mins[retIdx].getWeight()] );
            } // end i

            std::vector<int> recvOctCnts( npes );
            meshComm.allToAll<int>( 1, &( sendOctCnts[0] ), &( recvOctCnts[0] ) );

            std::vector<int> sendOctDisps( npes );
            std::vector<int> recvOctDisps( npes );
            sendOctDisps[0] = 0;
            recvOctDisps[0] = 0;
            for ( int i = 1; i < npes; ++i ) {
                sendOctDisps[i] = sendOctDisps[i - 1] + sendOctCnts[i - 1];
                recvOctDisps[i] = recvOctDisps[i - 1] + recvOctCnts[i - 1];
            } // end i

            std::vector<int> sendEidCnts( npes, 0 );
            for ( int i = 0; i < npes; ++i ) {
                for ( int j = 0; j < sendOctCnts[i]; ++j ) {
                    sendEidCnts[i] += ( tmpNodeList[sendOctDisps[i] + j].getWeight() );
                } // end j
            }     // end i

            std::vector<ot::TreeNode> recvOctList( recvOctDisps[npes - 1] + recvOctCnts[npes - 1] );
            ot::TreeNode *tmpNodeListPtr = NULL;
            if ( !( tmpNodeList.empty() ) ) {
                tmpNodeListPtr = &( tmpNodeList[0] );
            }
            ot::TreeNode *recvOctListPtr = NULL;
            if ( !( recvOctList.empty() ) ) {
                recvOctListPtr = &( recvOctList[0] );
            }
            meshComm.allToAll<ot::TreeNode>( tmpNodeListPtr,
                                             &( sendOctCnts[0] ),
                                             &( sendOctDisps[0] ),
                                             recvOctListPtr,
                                             &( recvOctCnts[0] ),
                                             &( recvOctDisps[0] ),
                                             true );
            tmpNodeList.clear();
            sendOctCnts.clear();
            sendOctDisps.clear();

            std::vector<int> recvEidCnts( npes, 0 );
            for ( int i = 0; i < npes; ++i ) {
                for ( int j = 0; j < recvOctCnts[i]; ++j ) {
                    recvEidCnts[i] += ( recvOctList[recvOctDisps[i] + j].getWeight() );
                } // end j
            }     // end i

            std::vector<int> sendEidDisps( npes );
            std::vector<int> recvEidDisps( npes );
            sendEidDisps[0] = 0;
            recvEidDisps[0] = 0;
            for ( int i = 1; i < npes; ++i ) {
                sendEidDisps[i] = sendEidDisps[i - 1] + sendEidCnts[i - 1];
                recvEidDisps[i] = recvEidDisps[i - 1] + recvEidCnts[i - 1];
            } // end i

            std::vector<int> recvEidList( recvEidDisps[npes - 1] + recvEidCnts[npes - 1] );
            int *sendEidListPtr = NULL;
            if ( !( sendEidList.empty() ) ) {
                sendEidListPtr = &( sendEidList[0] );
            }
            int *recvEidListPtr = NULL;
            if ( !( recvEidList.empty() ) ) {
                recvEidListPtr = &( recvEidList[0] );
            }
            meshComm.allToAll<int>( sendEidListPtr,
                                    &( sendEidCnts[0] ),
                                    &( sendEidDisps[0] ),
                                    recvEidListPtr,
                                    &( recvEidCnts[0] ),
                                    &( recvEidDisps[0] ),
                                    true );
            sendEidList.clear();
            sendEidCnts.clear();
            sendEidDisps.clear();
            recvEidCnts.clear();
            recvEidDisps.clear();

            std::vector<std::vector<int>> tmpEidList( recvOctList.size() );
            for ( size_t i = 0, j = 0; i < recvOctList.size(); ++i ) {
                tmpEidList[i].resize( recvOctList[i].getWeight() );
                for ( size_t k = 0; k < tmpEidList[i].size(); ++k, ++j ) {
                    tmpEidList[i][k] = recvEidList[j];
                } // end k
            }     // end i
            recvEidList.clear();

            std::vector<std::vector<int>> dummyElemIdList( d_nodeList.size() );
            std::vector<std::vector<int>> dummyRankList( d_nodeList.size() );

            for ( int i = 0; i < npes; ++i ) {
                for ( int j = 0; j < recvOctCnts[i]; ++j ) {
                    unsigned int retIdx = 0;
                    bool found          = seq::maxLowerBound<ot::TreeNode>(
                        d_nodeList, recvOctList[recvOctDisps[i] + j], retIdx, NULL, NULL );
                    AMP_CHECK_ASSERT( found );
                    AMP_CHECK_ASSERT( d_nodeList[retIdx] == recvOctList[recvOctDisps[i] + j] );
                    dummyElemIdList[retIdx].insert( dummyElemIdList[retIdx].end(),
                                                    tmpEidList[recvOctDisps[i] + j].begin(),
                                                    tmpEidList[recvOctDisps[i] + j].end() );
                    dummyRankList[retIdx].insert(
                        dummyRankList[retIdx].end(), tmpEidList[recvOctDisps[i] + j].size(), i );
                } // end j
            }     // end i
            recvOctCnts.clear();
            recvOctDisps.clear();
            recvOctList.clear();

            for ( size_t i = 0; i < d_nodeList.size(); ++i ) {
                d_nodeList[i].setWeight( dummyElemIdList[i].size() );
                d_elemIdList.insert(
                    d_elemIdList.end(), dummyElemIdList[i].begin(), dummyElemIdList[i].end() );
                d_rankList.insert(
                    d_rankList.end(), dummyRankList[i].begin(), dummyRankList[i].end() );
            } // end i
            dummyElemIdList.clear();
            dummyRankList.clear();
        }
    }

    d_stIdxList.resize( d_nodeList.size() );
    if ( !( d_nodeList.empty() ) ) {
        d_stIdxList[0] = 0;
        for ( size_t i = 1; i < d_nodeList.size(); ++i ) {
            d_stIdxList[i] = d_stIdxList[i - 1] + d_nodeList[i - 1].getWeight();
        } // end i
    }

    if ( d_verbose ) {
        int numFinalLocalOcts  = d_nodeList.size();
        int numFinalGlobalOcts = meshComm.sumReduce( numFinalLocalOcts );
        meshComm.barrier();
        if ( !rank ) {
            d_oStream << "Total num final octants = " << numFinalGlobalOcts << std::endl;
        }
    }

    if ( d_verbose ) {
        meshComm.barrier();
    }
    setupEndTime                = MPI_Wtime();
    d_timingMeasurements[Setup] = setupEndTime - setupBeginTime;

    if ( d_verbose ) {
        meshComm.barrier();
        if ( !rank ) {
            d_oStream << "Finished setting up DS for search in "
                      << ( setupEndTime - setupBeginTime ) << " seconds." << std::endl;
        }
    }
}


void DendroSearch::search( AMP::AMP_MPI comm, const std::vector<double> &pts )
{
    int rank = comm.getRank();
    int npes = comm.getSize();

    //    std::string fileName = "debug_dendro_" + AMP::Utilities::intToString(rank);
    //    std::fstream d_fout;
    //    d_fout.open(fileName.c_str(), std::fstream::out);
    //    d_fout<<"local elements="<<(d_meshAdapter.get() != NULL ?
    //    static_cast<int>(d_meshAdapter->numLocalElements(AMP::Mesh::GeomType::Volume)) : -1)
    //        <<"  global="<<(d_meshAdapter.get() != NULL ?
    //        static_cast<int>(d_meshAdapter->numGlobalElements(AMP::Mesh::GeomType::Volume)) :
    //        -1)<<"\n";

    double coarseSearchBeginTime = MPI_Wtime();

    std::vector<int> rankMap( npes );

    int myRank = -1;
    if ( d_meshAdapter != NULL ) {
        AMP::AMP_MPI meshComm = d_meshAdapter->getComm();
        myRank                = meshComm.getRank();
    }

    comm.allGather( myRank, &( rankMap[0] ) );

    std::vector<int> invRankMap( npes, -1 );
    for ( int i = 0; i < npes; ++i ) {
        if ( rankMap[i] >= 0 ) {
            invRankMap[rankMap[i]] = i;
        }
    } // end i
    //    d_fout<<"rankmap=";
    //    for (int i = 0; i < npes; ++i) {
    //        d_fout<<rankMap[i]<<"  ";
    //    }
    //    d_fout<<"\n";
    //    d_fout<<"invRankmap=";
    //    for (int i = 0; i < npes; ++i) {
    //        d_fout<<invRankMap[i]<<"  ";
    //    }
    //    d_fout<<"\n";

    std::vector<double> bcastBuff( 7 );
    if ( myRank == 0 ) {
        bcastBuff[0] = static_cast<double>( d_mins.size() );
        std::copy( d_minCoords.begin(), d_minCoords.end(), &( bcastBuff[1] ) );
        std::copy( d_scalingFactor.begin(), d_scalingFactor.end(), &( bcastBuff[4] ) );
    }

    comm.bcast( &( bcastBuff[0] ), 7, invRankMap[0] );

    int minsSize = static_cast<int>( bcastBuff[0] );
    std::copy( &( bcastBuff[1] ), &( bcastBuff[1] ) + 3, d_minCoords.begin() );
    std::copy( &( bcastBuff[4] ), &( bcastBuff[4] ) + 3, d_scalingFactor.begin() );
    bcastBuff.clear();

    d_mins.resize( minsSize );
    comm.bcast( &( d_mins[0] ), minsSize, invRankMap[0] );

    d_numLocalPts = ( pts.size() ) / 3;
    //    d_fout<<"d_numLocalPts="<<d_numLocalPts<<"\n";

    double searchTime[7] = { 0, 0, 0, 0, 0, 0, 0 };

    if ( d_verbose ) {
        comm.barrier();
        searchTime[0] = MPI_Wtime();
    }

    const unsigned int MaxDepth = 30;
    const unsigned int ITPMD    = ( 1u << MaxDepth );
    const double DTPMD          = static_cast<double>( ITPMD );

    std::vector<ot::NodeAndValues<double, 4>> ptsWrapper( d_numLocalPts );
    for ( int i = 0; i < d_numLocalPts; ++i ) {
        double x        = pts[3 * i];
        double y        = pts[( 3 * i ) + 1];
        double z        = pts[( 3 * i ) + 2];
        double scaledX  = ( ( x - d_minCoords[0] ) * d_scalingFactor[0] );
        double scaledY  = ( ( y - d_minCoords[1] ) * d_scalingFactor[1] );
        double scaledZ  = ( ( z - d_minCoords[2] ) * d_scalingFactor[2] );
        unsigned int pX = static_cast<unsigned int>( scaledX * DTPMD );
        unsigned int pY = static_cast<unsigned int>( scaledY * DTPMD );
        unsigned int pZ = static_cast<unsigned int>( scaledZ * DTPMD );

        ptsWrapper[i].node = ot::TreeNode( pX, pY, pZ, MaxDepth, 3, MaxDepth );
        ptsWrapper[i].node.setWeight( rank );
        ptsWrapper[i].values[0] = x;
        ptsWrapper[i].values[1] = y;
        ptsWrapper[i].values[2] = z;
        ptsWrapper[i].values[3] = i;
    } // end i

    if ( d_verbose ) {
        comm.barrier();
        searchTime[1] = MPI_Wtime();
        if ( !rank ) {
            std::cout << "Time for step-1 of search: " << searchTime[1] - searchTime[0]
                      << " seconds." << std::endl;
        }
    }

    // PERFORMANCE IMPROVEMENT: Should PtsWrapper be sorted or not?
    // If PtsWrapper is sorted (even just a local sort), we can skip the
    // binary searches and use binning instead.  Binning is amortized constant
    // time and using binary searches would be logarithmic. This is just a matter
    // of constants since sorting is also logarithmic.

    d_sendCnts.resize( npes );
    d_recvCnts.resize( npes );
    d_sendDisps.resize( npes );
    d_recvDisps.resize( npes );
    std::fill( d_sendCnts.begin(), d_sendCnts.end(), 0 );

    std::vector<int> part( d_numLocalPts, -1 );

    for ( int i = 0; i < d_numLocalPts; ++i ) {
        unsigned int retIdx = 0;
        bool found =
            seq::maxLowerBound<ot::TreeNode>( d_mins, ( ptsWrapper[i].node ), retIdx, NULL, NULL );
        if ( found ) {
            part[i] = invRankMap[d_mins[retIdx].getWeight()];
            ++( d_sendCnts[part[i]] );
        }
    } // end i

    if ( d_verbose ) {
        comm.barrier();
        searchTime[2] = MPI_Wtime();
        if ( !rank ) {
            std::cout << "Time for step-2 of search: " << searchTime[2] - searchTime[1]
                      << " seconds." << std::endl;
        }
    }

    comm.allToAll<int>( 1, &( d_sendCnts[0] ), &( d_recvCnts[0] ) );

    d_sendDisps[0] = 0;
    d_recvDisps[0] = 0;
    for ( int i = 1; i < npes; ++i ) {
        d_sendDisps[i] = d_sendDisps[i - 1] + d_sendCnts[i - 1];
        d_recvDisps[i] = d_recvDisps[i - 1] + d_recvCnts[i - 1];
    } // end i

    std::vector<ot::NodeAndValues<double, 4>> sendList( d_sendDisps[npes - 1] +
                                                        d_sendCnts[npes - 1] );

    std::fill( d_sendCnts.begin(), d_sendCnts.end(), 0 );

    for ( int i = 0; i < d_numLocalPts; ++i ) {
        if ( part[i] >= 0 ) {
            sendList[d_sendDisps[part[i]] + d_sendCnts[part[i]]] = ptsWrapper[i];
            ++( d_sendCnts[part[i]] );
        }
    } // end i
    ptsWrapper.clear();

    std::vector<ot::NodeAndValues<double, 4>> recvList( d_recvDisps[npes - 1] +
                                                        d_recvCnts[npes - 1] );
    //    d_fout<<"sendList.size()="<<sendList.size()<<"\n";
    //    d_fout<<"recvList.size()="<<recvList.size()<<"\n";
    comm.allToAll( ( !( sendList.empty() ) ? &( sendList[0] ) : NULL ),
                   &( d_sendCnts[0] ),
                   &( d_sendDisps[0] ),
                   ( !( recvList.empty() ) ? &( recvList[0] ) : NULL ),
                   &( d_recvCnts[0] ),
                   &( d_recvDisps[0] ),
                   true );
    sendList.clear();

    if ( d_verbose ) {
        comm.barrier();
        searchTime[3] = MPI_Wtime();
        if ( !rank ) {
            std::cout << "Time for step-3 of search: " << searchTime[3] - searchTime[2]
                      << " seconds." << std::endl;
        }
    }

    std::fill( d_sendCnts.begin(), d_sendCnts.end(), 0 );

    rank                       = -1;
    double fineSearchBeginTime = MPI_Wtime();
    if ( d_meshAdapter.get() != NULL ) {
        AMP::AMP_MPI meshComm = d_meshAdapter->getComm();
        rank                  = meshComm.getRank();
        npes                  = meshComm.getSize();

        std::vector<int> ptToOctMap( ( recvList.size() ), -1 );
        for ( size_t i = 0; i < recvList.size(); ++i ) {
            unsigned int retIdx = 0;
            seq::maxLowerBound<ot::TreeNode>(
                d_nodeList, ( recvList[i].node ), retIdx, NULL, NULL );
            if ( d_nodeList[retIdx].isAncestor( recvList[i].node ) ) {
                ptToOctMap[i] = retIdx;
                int stIdx     = d_stIdxList[retIdx];
                for ( size_t j = 0; j < d_nodeList[retIdx].getWeight(); ++j ) {
                    d_sendCnts[d_rankList[stIdx + j]]++;
                } // end j
            }
        } // end i

        if ( d_verbose ) {
            meshComm.barrier();
            searchTime[4] = MPI_Wtime();
            if ( !rank ) {
                std::cout << "Time for step-4 of search: " << searchTime[4] - searchTime[3]
                          << " seconds." << std::endl;
            }
        }

        meshComm.allToAll( 1, &( d_sendCnts[0] ), &( d_recvCnts[0] ) );

        d_sendDisps[0] = 0;
        d_recvDisps[0] = 0;
        for ( int i = 1; i < npes; ++i ) {
            d_sendDisps[i] = d_sendDisps[i - 1] + d_sendCnts[i - 1];
            d_recvDisps[i] = d_recvDisps[i - 1] + d_recvCnts[i - 1];
        } // end i

        std::vector<double> sendPtsList( 6 * ( d_sendDisps[npes - 1] + d_sendCnts[npes - 1] ) );

        std::fill( d_sendCnts.begin(), d_sendCnts.end(), 0 );

        for ( size_t i = 0; i < ptToOctMap.size(); ++i ) {
            if ( ptToOctMap[i] >= 0 ) {
                int stIdx = d_stIdxList[ptToOctMap[i]];
                for ( size_t j = 0; j < d_nodeList[ptToOctMap[i]].getWeight(); ++j ) {
                    int recvRank = d_rankList[stIdx + j];
                    int currIdx  = 6 * ( d_sendDisps[recvRank] + d_sendCnts[recvRank] );
                    // Local Id of this element on the processor that owns this element
                    sendPtsList[currIdx] = d_elemIdList[stIdx + j];
                    // Pt's x coordinate
                    sendPtsList[currIdx + 1] = recvList[i].values[0];
                    // Pt's y coordinate
                    sendPtsList[currIdx + 2] = recvList[i].values[1];
                    // Pt's z coordinate
                    sendPtsList[currIdx + 3] = recvList[i].values[2];
                    // Local Id of Pt on the processor that owns this Pt
                    sendPtsList[currIdx + 4] = recvList[i].values[3];
                    // rank of processor that owns Pt
                    sendPtsList[currIdx + 5] = recvList[i].node.getWeight();
                    ++( d_sendCnts[recvRank] );
                } // end j
            }
        } // end i
        recvList.clear();

        for ( int i = 0; i < npes; ++i ) {
            d_sendCnts[i] *= 6;
            d_sendDisps[i] *= 6;
            d_recvCnts[i] *= 6;
            d_recvDisps[i] *= 6;
        } // end i

        std::vector<double> recvPtsList( d_recvDisps[npes - 1] + d_recvCnts[npes - 1] );
        //        d_fout<<"sendPtsList.size()="<<sendPtsList.size()<<"\n";
        //        d_fout<<"recvPtsList.size()="<<recvPtsList.size()<<"\n";
        meshComm.allToAll( ( !( sendPtsList.empty() ) ? &( sendPtsList[0] ) : NULL ),
                           &( d_sendCnts[0] ),
                           &( d_sendDisps[0] ),
                           ( !( recvPtsList.empty() ) ? &( recvPtsList[0] ) : NULL ),
                           &( d_recvCnts[0] ),
                           &( d_recvDisps[0] ),
                           true );
        sendPtsList.clear();

        if ( d_verbose ) {
            meshComm.barrier();
            searchTime[5] = MPI_Wtime();
            if ( !rank ) {
                std::cout << "Time for step-5 of search: " << searchTime[5] - searchTime[4]
                          << " seconds." << std::endl;
            }
        }

        double coarseSearchEndTime         = MPI_Wtime();
        d_timingMeasurements[CoarseSearch] = coarseSearchEndTime - coarseSearchBeginTime;
        fineSearchBeginTime                = MPI_Wtime();

        std::fill( d_sendCnts.begin(), d_sendCnts.end(), 0 );

        int numRecvPts = recvPtsList.size() / 6;
        std::vector<double> tmpPtLocalCoord( 3 );

        d_foundPts.reserve( 6 * numRecvPts );
        unsigned int numFoundPts   = 0;
        bool coordinates_are_local = true;
        //        d_fout<<"numRecvPts="<<numRecvPts<<"\n";
        //        d_fout<<"d_volume_elements.size()="<<d_volume_elements.size()<<std::endl;
        for ( int i = 0; i < numRecvPts; ++i ) {
            double const *tmpPtGlobalCoordPtr = &( recvPtsList[6 * i] ) + 1;
            unsigned int eId                  = static_cast<unsigned int>( recvPtsList[6 * i] );
            unsigned int procId               = static_cast<unsigned int>( recvPtsList[6 * i + 5] );
            //            d_fout<<"i="<<i<<"  eid="<<eId<<"  procId="<<procId<<"
            //            loaclId="<<recvPtsList[6*i+4];
            if ( d_volume_elements[eId]->within_bounding_box( tmpPtGlobalCoordPtr, d_tolerance ) ) {
                //              d_fout<<"  bbox";
                if ( d_volume_elements[eId]->within_bounding_polyhedron( tmpPtGlobalCoordPtr,
                                                                         d_tolerance ) ) {
                    //                  d_fout<<"  bhedron";
                    d_volume_elements[eId]->map_global_to_local( tmpPtGlobalCoordPtr,
                                                                 &( tmpPtLocalCoord[0] ) );
                    if ( d_volume_elements[eId]->contains_point(
                             &( tmpPtLocalCoord[0] ), coordinates_are_local, d_tolerance ) ) {
                        //                      d_fout<<"  ###";
                        d_foundPts.push_back( recvPtsList[6 * i] );
                        for ( unsigned int d = 0; d < 3; ++d ) {
                            d_foundPts.push_back( tmpPtLocalCoord[d] );
                        }
                        d_foundPts.push_back( recvPtsList[6 * i + 4] );
                        d_foundPts.push_back( recvPtsList[6 * i + 5] );
                        ++numFoundPts;
                        ++( d_sendCnts[procId] );
                    } // end if
                }     // end if
            }         // end if
                      //            d_fout<<"\n";
                      //            if( (static_cast<unsigned int>(recvPtsList[6*i+4]) == 3)
                      //                || (static_cast<unsigned int>(recvPtsList[6*i+4]) == 4)
                      //                || (static_cast<unsigned int>(recvPtsList[6*i+4]) == 2)) {
                      //              double point_of_view[3] = { 1.0, 1.0, 1.0 };
            //              draw_point(tmpPtGlobalCoordPtr, "red", std::cout, "$\\diamond$");
            //              draw_hex8_element(d_volume_elements[eId], point_of_view, std::cout);
            //            } // end if
        } // end i
        recvPtsList.clear();
        //        d_fout<<"myRank="<<rank<<"\n";
        npes = comm.getSize();
    } // end if

    comm.allToAll( 1, &( d_sendCnts[0] ), &( d_recvCnts[0] ) );

    d_sendDisps[0] = 0;
    d_recvDisps[0] = 0;
    for ( int i = 1; i < npes; ++i ) {
        d_sendDisps[i] = d_sendDisps[i - 1] + d_sendCnts[i - 1];
        d_recvDisps[i] = d_recvDisps[i - 1] + d_recvCnts[i - 1];
    } // end i
    //    d_fout<<"d_sendCnts=";
    //    for (int i = 0; i < npes; ++i) {
    //        d_fout<<d_sendCnts[i]<<"  ";
    //    } // end for i
    //    d_fout<<"\n";
    //    d_fout<<"d_recvCnts=";
    //    for (int i = 0; i < npes; ++i) {
    //        d_fout<<d_recvCnts[i]<<"  ";
    //    } // end for i
    //    d_fout<<"\n";

    if ( d_verbose ) {
        comm.barrier();
        searchTime[6] = MPI_Wtime();
        if ( !rank ) {
            std::cout << "Time for step-6 of search: " << searchTime[6] - searchTime[5]
                      << " seconds." << std::endl;
        }
    }

    double fineSearchEndTime         = MPI_Wtime();
    d_timingMeasurements[FineSearch] = fineSearchEndTime - fineSearchBeginTime;
    comm.bcast( &( d_timingMeasurements[0] ), numTimingTypes, invRankMap[0] );
}


void DendroSearch::interpolate( AMP::AMP_MPI comm,
                                AMP::LinearAlgebra::Vector::const_shared_ptr vectorField,
                                const unsigned int dofsPerNode,
                                std::vector<double> &results,
                                std::vector<bool> &foundPt )
{
    const int rank = comm.getRank();
    const int npes = comm.getSize();

    double interpolateBeginTime, interpolateStep1Time = 0., interpolateStep2Time;
    if ( d_verbose ) {
        comm.barrier();
    }
    interpolateBeginTime = MPI_Wtime();

    AMP_CHECK_ASSERT( vectorField->getUpdateStatus() ==
                      AMP::LinearAlgebra::VectorData::UpdateState::UNCHANGED );
    std::shared_ptr<AMP::Discretization::DOFManager> dofManager = vectorField->getDOFManager();

    for ( int i = 0; i < npes; ++i ) {
        d_sendCnts[i] *= ( dofsPerNode + 1 );
        d_recvCnts[i] *= ( dofsPerNode + 1 );
        d_sendDisps[i] *= ( dofsPerNode + 1 );
        d_recvDisps[i] *= ( dofsPerNode + 1 );
    } // end for i

    std::vector<double> sendResults( d_sendDisps[npes - 1] + d_sendCnts[npes - 1] );

    std::vector<int> tmpSendCnts( npes, 0 );

    std::vector<double> basis_functions_values( 8 );
    for ( size_t i = 0; i < d_foundPts.size(); i += 6 ) {
        unsigned int elemLocalId = static_cast<unsigned int>( d_foundPts[i] );
        std::vector<AMP::Mesh::MeshElement> amp_vector_support_points =
            d_localElems[elemLocalId].getElements( AMP::Mesh::GeomType::Vertex );
        hex8_element_t::get_basis_functions_values( &( d_foundPts[i + 1] ),
                                                    &( basis_functions_values[0] ) );

        std::vector<double> value( dofsPerNode, 0.0 );
        for ( unsigned int j = 0; j < 8; ++j ) {
            std::vector<size_t> globalID;
            dofManager->getDOFs( amp_vector_support_points[j].globalID(), globalID );
            AMP_CHECK_ASSERT( globalID.size() == dofsPerNode );
            for ( size_t d = 0; d < dofsPerNode; ++d ) {
                double vecVal = vectorField->getValueByGlobalID( globalID[d] );
                value[d] += ( vecVal * basis_functions_values[j] );
            } // end d
        }     // end j
        unsigned int ptProcId = static_cast<unsigned int>( d_foundPts[i + 5] );
        sendResults[d_sendDisps[ptProcId] + tmpSendCnts[ptProcId]] = d_foundPts[i + 4];
        ++( tmpSendCnts[ptProcId] );
        for ( size_t d = 0; d < dofsPerNode; ++d ) {
            sendResults[d_sendDisps[ptProcId] + tmpSendCnts[ptProcId]] = value[d];
            ++( tmpSendCnts[ptProcId] );
        } // end d
    }     // end i
    tmpSendCnts.clear();

    if ( d_verbose ) {
        comm.barrier();
        interpolateStep1Time = MPI_Wtime();
        if ( !rank ) {
            std::cout << "Time for step-1 of interpolate: "
                      << ( interpolateStep1Time - interpolateBeginTime ) << " seconds."
                      << std::endl;
        }
    }

    std::vector<double> recvResults( d_recvDisps[npes - 1] + d_recvCnts[npes - 1] );

    comm.allToAll( ( !( sendResults.empty() ) ? &( sendResults[0] ) : NULL ),
                   &( d_sendCnts[0] ),
                   &( d_sendDisps[0] ),
                   ( !( recvResults.empty() ) ? &( recvResults[0] ) : NULL ),
                   &( d_recvCnts[0] ),
                   &( d_recvDisps[0] ),
                   true );
    sendResults.clear();

    // Points that are not found will have a result = 0.
    results.resize( dofsPerNode * d_numLocalPts );
    for ( unsigned int i = 0; i < results.size(); ++i ) {
        results[i] = 0.0;
    } // end for i

    foundPt.resize( d_numLocalPts );
    for ( unsigned int i = 0; i < foundPt.size(); ++i ) {
        foundPt[i] = static_cast<bool>( NotFound );
    } // end for i

    for ( size_t i = 0; i < recvResults.size(); i += ( dofsPerNode + 1 ) ) {
        unsigned int locId = static_cast<unsigned int>( recvResults[i] );
        foundPt[locId]     = static_cast<bool>( Found );
        for ( size_t d = 0; d < dofsPerNode; ++d ) {
            results[( locId * dofsPerNode ) + d] = recvResults[i + d + 1];
        } // end d
    }     // end i

    for ( int i = 0; i < npes; ++i ) {
        d_sendCnts[i] /= ( dofsPerNode + 1 );
        d_recvCnts[i] /= ( dofsPerNode + 1 );
        d_sendDisps[i] /= ( dofsPerNode + 1 );
        d_recvDisps[i] /= ( dofsPerNode + 1 );
    } // end for i

    if ( d_verbose ) {
        comm.barrier();
    }
    interpolateStep2Time                = MPI_Wtime();
    d_timingMeasurements[Interpolation] = interpolateStep2Time - interpolateBeginTime;
    if ( d_verbose ) {
        if ( !rank ) {
            std::cout << "Time for step-2 of interpolate: "
                      << ( interpolateStep2Time - interpolateStep1Time ) << " seconds."
                      << std::endl;
        }
    }
}


void DendroSearch::reportTiming( size_t n,
                                 TimingType const *timingTypes,
                                 double *timingMeasurements )
{
    AMP_INSIST( !d_verbose,
                "Using verbose mode in DendroSearch. This involves calls to "
                "MPI_Barrier. So timing measurements are bad!" );
    for ( size_t i = 0; i < n; ++i ) {
        timingMeasurements[i] = d_timingMeasurements[timingTypes[i]];
    } // end for i
}


} // namespace Mesh
} // namespace AMP
