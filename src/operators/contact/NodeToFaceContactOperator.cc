
#include "operators/mechanics/IsotropicElasticModel.h"
#include "operators/contact/NodeToFaceContactOperator.h"
#include "ampmesh/dendro/DendroSearch.h"
#include "ampmesh/MeshID.h"
#include "ampmesh/euclidean_geometry_tools.h"

#include <iterator>
#include <algorithm>
#include <cmath>
#include <iomanip>


namespace AMP {
namespace Operator {


void NodeToFaceContactOperator::initialize() 
{
    d_ActiveSet.clear();
    d_InactiveSet.clear();
    /** get all slave boundary vertices and tag them as inactive */
    AMP::Mesh::Mesh::shared_ptr slaveMesh = d_Mesh->Subset(d_SlaveMeshID);
    if (slaveMesh.get() != NULL) {
        AMP::Mesh::MeshIterator slaveMeshIterator = slaveMesh->getBoundaryIDIterator(AMP::Mesh::Vertex, d_SlaveBoundaryID);
        AMP::Mesh::MeshIterator slaveMeshIterator_begin = slaveMeshIterator.begin(), 
          slaveMeshIterator_end = slaveMeshIterator.end();
        size_t const nSlaveVertices = slaveMeshIterator.size();
        d_ActiveSet.reserve(nSlaveVertices);
        d_InactiveSet.reserve(nSlaveVertices);
        for (slaveMeshIterator = slaveMeshIterator_begin; slaveMeshIterator != slaveMeshIterator_end; ++slaveMeshIterator) {
          d_InactiveSet.push_back(slaveMeshIterator->globalID());
        } // end loop over the slave vertices on boundary
    } //end if
}


size_t NodeToFaceContactOperator::updateActiveSet(AMP::LinearAlgebra::Vector::shared_ptr displacementFieldVector, bool skipDisplaceMesh) 
{
    size_t const npes = d_GlobalComm.getSize();
    size_t const rank = d_GlobalComm.getRank();

    std::cout<<"### "<<rank<<" / "<<npes<<" ###"<<std::endl;

    size_t nActiveSlaveVerticesDeactivated = 0;
    size_t nInactiveSlaveVerticesActivated = 0;

    size_t const activeSetSizeInitial = d_GlobalComm.sumReduce(d_ActiveSet.size());
    size_t const inactiveSetSizeInitial = d_GlobalComm.sumReduce(d_InactiveSet.size());
    if (!rank) {
        std::cout<<"d_ActiveSet.size()="<<activeSetSizeInitial<<"  "
            <<"d_InactiveSet.size()="<<inactiveSetSizeInitial<<"\n";
    } // end if

    /** displace the mesh */
    if (!skipDisplaceMesh) {
        d_Mesh->displaceMesh(displacementFieldVector);
    } // end if

    /** get the inactive slave vertices coordinates */
    size_t const nInactiveSlaveVertices = d_InactiveSet.size();
    std::vector<double> inactiveSlaveVerticesCoord(3*nInactiveSlaveVertices); 
    AMP::Mesh::MeshElement tmpInactiveSlaveVertex;
    std::vector<double> tmpInactiveSlaveVertexCoord(3);
    for (size_t i = 0; i < nInactiveSlaveVertices; ++i) {
        tmpInactiveSlaveVertex = d_Mesh->getElement(d_InactiveSet[i]);
        tmpInactiveSlaveVertexCoord = tmpInactiveSlaveVertex.coord();
        std::copy(tmpInactiveSlaveVertexCoord.begin(), tmpInactiveSlaveVertexCoord.end(), &(inactiveSlaveVerticesCoord[3*i]));
    } // end for i
    /** perform a dendro search for all inactive slave vertices over the master mesh */
    AMP::Mesh::Mesh::shared_ptr masterMesh = d_Mesh->Subset(d_MasterMeshID);
    // TODO: read dummyVerboseFlag from input file
    bool const dummyVerboseFlag = false;
    AMP::Mesh::DendroSearch dendroSearchOnMaster(masterMesh, dummyVerboseFlag);
    dendroSearchOnMaster.setTolerance(1.0e-10);
    dendroSearchOnMaster.search(d_GlobalComm, inactiveSlaveVerticesCoord);

    std::vector<AMP::Mesh::MeshElementID> tmpMasterVerticesGlobalIDs;
    std::vector<double> tmpSlaveVerticesShift, tmpSlaveVerticesLocalCoordOnFace;
    std::vector<int> flags;
    std::vector<AMP::Mesh::MeshElementID> tmpMasterVolumesGlobalIDs;
    std::vector<size_t> tmpMasterFacesLocalIndices;
    dendroSearchOnMaster.projectOnBoundaryID(d_GlobalComm, d_MasterBoundaryID,
        tmpMasterVerticesGlobalIDs, tmpSlaveVerticesShift, tmpSlaveVerticesLocalCoordOnFace, flags, 
        tmpMasterVolumesGlobalIDs, tmpMasterFacesLocalIndices);


    size_t const nActiveSlaveVertices = d_ActiveSet.size();

    if (d_GlobalComm.anyReduce(!d_ActiveSet.empty())) {
        /** compute the normal vectors at slave nodes */

        // to be moved in the class private data
        std::vector<AMP::Mesh::MeshElementID> d_RecvMasterVolumesGlobalIDs;
        std::vector<size_t> d_RecvMasterFacesLocalIndices;

        std::vector<int> sendCnts(npes, 0);
        for (size_t i = 0; i < nActiveSlaveVertices; ++i) {
          ++sendCnts[d_MasterVolumesGlobalIDs[i].owner_rank()];
        } // end for i
        std::vector<int> sendDisps(npes, 0);
        sendDisps[0] = 0;
        for (size_t i = 1; i < npes; ++i) {
          sendDisps[i] = sendDisps[i-1] + sendCnts[i-1];
        } // end for i
        size_t const nSendData = static_cast<size_t>(sendDisps[npes-1] + sendCnts[npes-1]);
        AMP_ASSERT( nSendData == nActiveSlaveVertices );
        std::vector<ProjectionData> sendProjectionDataBuffer(nSendData);
        std::fill(sendCnts.begin(), sendCnts.end(), 0);
        std::vector<size_t> sendMap(nSendData, nSendData);
        std::vector<size_t> recvMap(nSendData, nSendData);
        for (size_t i = 0; i < nSendData; ++i) {
            size_t sendToRank = d_MasterVolumesGlobalIDs[i].owner_rank();
            sendMap[i] = static_cast<size_t>(sendDisps[sendToRank]+sendCnts[sendToRank]);
            recvMap[sendMap[i]] = i;
            sendProjectionDataBuffer[sendMap[i]].d_MasterVolumeGlobalID = d_MasterVolumesGlobalIDs[i];
            sendProjectionDataBuffer[sendMap[i]].d_MasterFaceLocalIndex = d_MasterFacesLocalIndices[i];
            hex8_element_t::get_local_coordinates_on_face(&(d_MasterShapeFunctionsValues[4*i]), sendProjectionDataBuffer[sendMap[i]].d_SlaveVertexLocalCoordOnMasterFace);
            ++sendCnts[sendToRank];
        } // end for i
        std::vector<int> recvCnts(npes, 0);
        d_GlobalComm.allToAll(1, &(sendCnts[0]), &(recvCnts[0]));
        std::vector<int> recvDisps(npes, 0);
        recvDisps[0] = 0;
        for (size_t i = 1; i < npes; ++i) {
            recvDisps[i] = recvDisps[i-1] + recvCnts[i-1];
        } // end for i
        size_t const nRecvData = static_cast<size_t>(recvDisps[npes-1] + recvCnts[npes-1]);
        std::vector<ProjectionData> recvProjectionDataBuffer(nRecvData);
        d_GlobalComm.allToAll((!(sendProjectionDataBuffer.empty()) ? &(sendProjectionDataBuffer[0]) : NULL), &(sendCnts[0]), &(sendDisps[0]),
            (!(recvProjectionDataBuffer.empty()) ? &(recvProjectionDataBuffer[0]) : NULL), &(recvCnts[0]), &(recvDisps[0]), true);
        sendProjectionDataBuffer.clear(); 

        std::vector<StressStateData> sendStressStateDataBuffer(nRecvData);
        for (size_t i = 0; i < nRecvData; ++i) {
//          std::cout<<"# "<<i<<" ";
            // retrieve master volume element then compute normal vector and strain tensor at slave vertex
            std::vector<AMP::Mesh::MeshElement> masterVolumeVertices = d_Mesh->getElement(recvProjectionDataBuffer[i].d_MasterVolumeGlobalID).getElements(AMP::Mesh::Vertex);
            std::vector<AMP::Mesh::MeshElementID> masterVolumeVerticesGlobalIDs(8);
            double masterVolumeVerticesCoordinates[24];
            AMP_ASSERT(masterVolumeVertices.size() == 8);
            for (size_t j = 0; j < masterVolumeVertices.size(); ++j) {
                std::vector<double> vertexCoord = masterVolumeVertices[j].coord();
                AMP_ASSERT(vertexCoord.size() == 3);
                std::copy(vertexCoord.begin(), vertexCoord.end(), &(masterVolumeVerticesCoordinates[3*j]));
                masterVolumeVerticesGlobalIDs[j] = masterVolumeVertices[j].globalID();
//            std::cout<<"v="<<j<<"  "<<masterVolumeVertices[j].globalID().local_id()<<"\n";
            } // end for j
            std::vector<AMP::Mesh::MeshElement> masterVolumeFaces = d_Mesh->getElement(recvProjectionDataBuffer[i].d_MasterVolumeGlobalID).getElements(AMP::Mesh::Face);
            AMP_ASSERT( masterVolumeFaces[recvProjectionDataBuffer[i].d_MasterFaceLocalIndex].isOnBoundary(d_MasterBoundaryID) );
            bool stopHere = false;
            unsigned int const facesNew[24] = {
                0, 1, 5, 4,
                1, 2, 6, 5,
                3, 2, 6, 7,
                0, 3, 7, 4,
                0, 1, 2, 3,
                4, 5, 6, 7
            };
            for (size_t f = 0; f < 6; ++f) {
                if (f != recvProjectionDataBuffer[i].d_MasterFaceLocalIndex) {
                    AMP_ASSERT( !masterVolumeFaces[f].isOnBoundary(d_MasterBoundaryID) );
                } else {
//              std::cout<<"###";
                AMP_ASSERT( f == 5 );
            } // end if
            std::vector<AMP::Mesh::MeshElement> masterFaceVertices = masterVolumeFaces[f].getElements(AMP::Mesh::Vertex);
            unsigned int const * faces = hex8_element_t::get_faces();
            for (size_t v = 0; v < 4; ++v) {
                bool justChecking = masterFaceVertices[v].globalID() == masterVolumeVerticesGlobalIDs[faces[4*f+v]];
//              std::cout<<"f="<<f<<"  v="<<v<<"  index="<<faces[4*f+v]<<"  -> "<<justChecking<<"  "<<masterFaceVertices[v].globalID().local_id()<<"  "<<masterVolumeVerticesGlobalIDs[faces[4*f+v]].local_id()<<"\n";
                AMP_ASSERT( masterFaceVertices[v].globalID().local_id() == masterVolumeVerticesGlobalIDs[facesNew[4*f+v]].local_id() );
//              AMP_ASSERT( masterFaceVertices[v].globalID().local_id() == masterVolumeVerticesGlobalIDs[faces[4*f+v]].local_id() );
                if (!justChecking) { stopHere = true; }
            } // end for v
        } // end for f
        if (stopHere) { AMP_ASSERT(false); }

        hex8_element_t masterVolumeElement(masterVolumeVerticesCoordinates);
        masterVolumeElement.compute_normal_to_face( recvProjectionDataBuffer[i].d_MasterFaceLocalIndex, 
            recvProjectionDataBuffer[i].d_SlaveVertexLocalCoordOnMasterFace, sendStressStateDataBuffer[i].d_SlaveVertexNormalVector);
        double slaveVertexLocalCoordinates[3];
        hex8_element_t::map_face_to_local( recvProjectionDataBuffer[i].d_MasterFaceLocalIndex, 
            recvProjectionDataBuffer[i].d_SlaveVertexLocalCoordOnMasterFace, slaveVertexLocalCoordinates);
        std::vector<size_t> displacementIndices;
        getVectorIndicesFromGlobalIDs(masterVolumeVerticesGlobalIDs, displacementIndices);
        AMP_ASSERT(displacementIndices.size() == 24);
        double displacementValues[24];
        displacementFieldVector->getLocalValuesByGlobalID(24, &(displacementIndices[0]), &(displacementValues[0]));
        double strainTensor[6];
        hex8_element_t::compute_strain_tensor(slaveVertexLocalCoordinates, displacementValues, strainTensor);
        // compute normal of the surface traction at slave vertex
        double constitutiveMatrix[6][6];
        double E = 1.0e6;//boost::dynamic_pointer_cast<AMP::Operator::IsotropicElasticModel>(d_MasterMechanicsMaterialModel)->getYoungsModulus();
        double Nu = 0.3;//boost::dynamic_pointer_cast<AMP::Operator::IsotropicElasticModel>(d_MasterMechanicsMaterialModel)->getPoissonsRatio();
        double G = E/(2.0*(1.0 + Nu));

        double c = G;
        double a = 2.0*c*(1.0 - Nu)/(1.0 - (2.0*Nu));
        double b = 2.0*c*Nu/(1.0 - (2.0*Nu));

        for(int ii = 0; ii < 6; ii++) {
          for(int jj = 0; jj < 6; jj++) {
            constitutiveMatrix[ii][jj] = 0.0;
          }//end for j
        }//end for i

        constitutiveMatrix[0][0] = a;
        constitutiveMatrix[1][1] = a;
        constitutiveMatrix[2][2] = a;

        constitutiveMatrix[0][1] = b;
        constitutiveMatrix[0][2] = b;
        constitutiveMatrix[1][0] = b;
        constitutiveMatrix[1][2] = b;
        constitutiveMatrix[2][0] = b;
        constitutiveMatrix[2][1] = b;

        constitutiveMatrix[3][3] = c;
        constitutiveMatrix[4][4] = c;
        constitutiveMatrix[5][5] = c;

        //  d_MasterMechanicsMaterialModel->getConstitutiveMatrix(constitutiveMatrix);
        double stressTensor[6];
        for (size_t ii = 0; ii < 6; ++ii) {
            stressTensor[ii] = 0.0;
            for (size_t jj = 0; jj < 6; ++jj) {
                stressTensor[ii] += constitutiveMatrix[ii][jj] * strainTensor[jj];
            } // end for jj
        } // end for ii
/*          std::cout<<"epsilon="
              <<strainTensor[0]<<"  "
              <<strainTensor[1]<<"  "
              <<strainTensor[2]<<"  "
              <<strainTensor[3]<<"  "
              <<strainTensor[4]<<"  "
              <<strainTensor[5]<<"  ";
          std::cout<<"sigma="
              <<stressTensor[0]<<"  "
              <<stressTensor[1]<<"  "
              <<stressTensor[2]<<"  "
              <<stressTensor[3]<<"  "
              <<stressTensor[4]<<"  "
              <<stressTensor[5]<<"  ";
*/

//          compute_stress_tensor(constitutiveMatrix, strainTensor, stressTensor);
        compute_traction(stressTensor, sendStressStateDataBuffer[i].d_SlaveVertexNormalVector, sendStressStateDataBuffer[i].d_SlaveVertexSurfaceTraction);
/*
          std::cout<<"t="
              <<sendStressStateDataBuffer[i].d_SlaveVertexSurfaceTraction[0]<<"  "
              <<sendStressStateDataBuffer[i].d_SlaveVertexSurfaceTraction[1]<<"  "
              <<sendStressStateDataBuffer[i].d_SlaveVertexSurfaceTraction[2]<<"  ";
          for (size_t i = 0; i < 8; ++i) {
            std::cout<<"u["<<i<<"]="
                <<displacementValues[3*i+0]<<"  "
                <<displacementValues[3*i+1]<<"  "
                <<displacementValues[3*i+2]<<"  ";
          }
          std::cout<<"x="
              <<slaveVertexLocalCoordinates[0]<<"  "
              <<slaveVertexLocalCoordinates[1]<<"  "
              <<slaveVertexLocalCoordinates[2]<<"  ";
          std::cout<<" #"<<std::endl;
*/
        } // end for i
        recvProjectionDataBuffer.clear();

        sendCnts.swap(recvCnts);
        sendDisps.swap(recvDisps);
        std::vector<StressStateData> recvStressStateDataBuffer(nSendData);
        d_GlobalComm.allToAll((!(sendStressStateDataBuffer.empty()) ? &(sendStressStateDataBuffer[0]) : NULL), &(sendCnts[0]), &(sendDisps[0]),
            (!(recvStressStateDataBuffer.empty()) ? &(recvStressStateDataBuffer[0]) : NULL), &(recvCnts[0]), &(recvDisps[0]), true);
        sendStressStateDataBuffer.clear(); 

        d_SlaveVerticesNormalVectors.resize(3*nActiveSlaveVertices);
        d_SlaveVerticesSurfaceTraction.resize(3*nActiveSlaveVertices);
        for (size_t i = 0; i < nSendData; ++i) {
            std::copy( recvStressStateDataBuffer[i].d_SlaveVertexNormalVector, 
                recvStressStateDataBuffer[i].d_SlaveVertexNormalVector+3, 
                &(d_SlaveVerticesNormalVectors[3*recvMap[i]]) );
            std::copy( recvStressStateDataBuffer[i].d_SlaveVertexSurfaceTraction, 
                recvStressStateDataBuffer[i].d_SlaveVertexSurfaceTraction+3, 
                &(d_SlaveVerticesSurfaceTraction[3*recvMap[i]]) );
        } // end for i

        std::vector<double>::iterator slaveVerticesShiftIterator = d_SlaveShift.begin();
        std::vector<AMP::Mesh::MeshElementID>::iterator masterVerticesGlobalIDsIterator = d_MasterVerticesGlobalIDs.begin();
        std::vector<double>::iterator masterShapeFunctionsValuesIterator = d_MasterShapeFunctionsValues.begin();
        std::vector<AMP::Mesh::MeshElementID>::iterator masterVolumesGlobalIDsIterator = d_MasterVolumesGlobalIDs.begin();
        std::vector<size_t>::iterator masterFacesLocalIndicesIterator = d_MasterFacesLocalIndices.begin();

        std::vector<AMP::Mesh::MeshElementID>::iterator activeSetIterator = d_ActiveSet.begin();
        for (size_t i = 0; i < nActiveSlaveVertices; ++i) {
            double nDotT = compute_scalar_product(&(d_SlaveVerticesNormalVectors[3*i]), &(d_SlaveVerticesSurfaceTraction[3*i]));
            std::cout<<i<<"  "<<std::setprecision(15)<<"n.t="<<nDotT<<"\n";
            std::cout<<"n="
                <<d_SlaveVerticesNormalVectors[3*i+0]<<"  "
                <<d_SlaveVerticesNormalVectors[3*i+1]<<"  "
                <<d_SlaveVerticesNormalVectors[3*i+2]<<"  "
                <<"t="
                <<d_SlaveVerticesSurfaceTraction[3*i+0]<<"  "
                <<d_SlaveVerticesSurfaceTraction[3*i+1]<<"  "
                <<d_SlaveVerticesSurfaceTraction[3*i+2]<<"\n";
            std::vector<double> xv2france = d_Mesh->getElement(*activeSetIterator).coord();
            std::cout<<"  ( "<<xv2france[0]<<" , "<<xv2france[1]<<" , "<<xv2france[2]<<" )\n";
            if (nDotT > 1.0e-14) { 
                std::cout<<"@@@@@@@@@@@@@@@@\n";
                std::cout<<"@@@@@@@@@@@@@@@@\n";
                std::cout<<"@@@@@@@@@@@@@@@@\n";
                std::cout<<"@@@@@@@@@@@@@@@@\n";
                ++nActiveSlaveVerticesDeactivated;
                d_InactiveSet.push_back(*activeSetIterator);
                activeSetIterator = d_ActiveSet.erase(activeSetIterator);
                slaveVerticesShiftIterator = d_SlaveShift.erase(slaveVerticesShiftIterator, slaveVerticesShiftIterator+3);
                masterVerticesGlobalIDsIterator = d_MasterVerticesGlobalIDs.erase(masterVerticesGlobalIDsIterator, masterVerticesGlobalIDsIterator+4);
                masterShapeFunctionsValuesIterator = d_MasterShapeFunctionsValues.erase(masterShapeFunctionsValuesIterator, masterShapeFunctionsValuesIterator+4);
                masterVolumesGlobalIDsIterator = d_MasterVolumesGlobalIDs.erase(masterVolumesGlobalIDsIterator);
                masterFacesLocalIndicesIterator = d_MasterFacesLocalIndices.erase(masterFacesLocalIndicesIterator);
            } else {
                ++activeSetIterator;
                std::advance(slaveVerticesShiftIterator, 3);
                std::advance(masterVerticesGlobalIDsIterator, 4);
                std::advance(masterShapeFunctionsValuesIterator, 4);
                ++masterVolumesGlobalIDsIterator;
                ++masterFacesLocalIndicesIterator;
            } // end if
        } // end for i
        AMP_ASSERT( activeSetIterator == d_ActiveSet.end() );
        AMP_ASSERT( slaveVerticesShiftIterator == d_SlaveShift.end() );
        AMP_ASSERT( masterVerticesGlobalIDsIterator == d_MasterVerticesGlobalIDs.end() );
        AMP_ASSERT( masterShapeFunctionsValuesIterator == d_MasterShapeFunctionsValues.end() );
        AMP_ASSERT( masterVolumesGlobalIDsIterator == d_MasterVolumesGlobalIDs.end() );
        AMP_ASSERT( masterFacesLocalIndicesIterator == d_MasterFacesLocalIndices.end() );
    } // end if
    // if (!rank) { std::cout<<"nActiveSlaveVerticesDeactivated="<<d_GlobalComm.sumReduce(nActiveSlaveVerticesDeactivated)<<"\n"; }
    std::cout<<"nActiveSlaveVerticesDeactivated="<<d_GlobalComm.sumReduce(nActiveSlaveVerticesDeactivated)<<"\n";

    size_t const nAddedConstraints = static_cast<size_t>(std::count(flags.begin(), flags.end(), AMP::Mesh::DendroSearch::FoundOnBoundary)); // std::count returns ptrdiff_t 
    size_t const nOlderConstraints = d_ActiveSet.size();

    ptrdiff_t localPtsNotFound = std::count(flags.begin(), flags.end(), AMP::Mesh::DendroSearch::NotFound);
    ptrdiff_t localPtsFoundNotOnBoundary = std::count(flags.begin(), flags.end(), AMP::Mesh::DendroSearch::FoundNotOnBoundary);
    ptrdiff_t localPtsFoundOnBoundary = std::count(flags.begin(), flags.end(), AMP::Mesh::DendroSearch::FoundOnBoundary);
    ptrdiff_t globalPtsNotFound = d_GlobalComm.sumReduce(localPtsNotFound);
    ptrdiff_t globalPtsFoundNotOnBoundary = d_GlobalComm.sumReduce(localPtsFoundNotOnBoundary);
    ptrdiff_t globalPtsFoundOnBoundary = d_GlobalComm.sumReduce(localPtsFoundOnBoundary);
    d_fout<<"Global number of points not found is "<<globalPtsNotFound<<" (local was "<<localPtsNotFound<<")"<<std::endl;
    d_fout<<"Global number of points found not on boundary is "<<globalPtsFoundNotOnBoundary<<" (local was "<<localPtsFoundNotOnBoundary<<")"<<std::endl;
    d_fout<<"Global number of points found on boundary is "<<globalPtsFoundOnBoundary<<" (local was "<<localPtsFoundOnBoundary<<")"<<std::endl;
    d_fout<<"Total number of points is "<<globalPtsNotFound+globalPtsFoundNotOnBoundary+globalPtsFoundOnBoundary<<std::endl;

    AMP_ASSERT( std::count(flags.begin(), flags.end(), AMP::Mesh::DendroSearch::FoundNotOnBoundary) == 0 ); // DendroSearch::FoundNotOnBoundary is not acceptable

    // TODO: bof bof
    size_t const nConstraints = nOlderConstraints + nAddedConstraints;

    // d_ActiveSet.resize(nOlderConstraints+nAddedConstraints, AMP::Mesh::MeshElementID());
    d_SlaveShift.resize(3*nConstraints, 0.0);
    d_MasterVerticesGlobalIDs.resize(4*nConstraints, AMP::Mesh::MeshElementID());
    d_MasterShapeFunctionsValues.resize(4*nConstraints, 0.0);
    d_MasterVolumesGlobalIDs.resize(nConstraints, AMP::Mesh::MeshElementID());
    d_MasterFacesLocalIndices.resize(nConstraints, 7);
      
    std::vector<double>::iterator slaveVerticesShiftIterator = d_SlaveShift.begin();
    std::vector<AMP::Mesh::MeshElementID>::iterator masterVerticesGlobalIDsIterator = d_MasterVerticesGlobalIDs.begin();
    double * masterShapeFunctionsValuesPointer = &(d_MasterShapeFunctionsValues[0]); // needed this way to work with hex8_element_t::get_basis_functions_values_on_face(...)
    std::vector<AMP::Mesh::MeshElementID>::iterator masterVolumesGlobalIDsIterator = d_MasterVolumesGlobalIDs.begin();
    std::vector<size_t>::iterator masterFacesLocalIndicesIterator = d_MasterFacesLocalIndices.begin();
    std::advance(slaveVerticesShiftIterator, 3*nOlderConstraints);
    std::advance(masterVerticesGlobalIDsIterator, 4*nOlderConstraints);
    std::advance(masterShapeFunctionsValuesPointer, 4*nOlderConstraints);
    std::advance(masterVolumesGlobalIDsIterator, nOlderConstraints);
    std::advance(masterFacesLocalIndicesIterator, nOlderConstraints);

    std::vector<AMP::Mesh::MeshElementID>::iterator inactiveSetIterator = d_InactiveSet.begin();
    for (size_t i = 0; i < nInactiveSlaveVertices; ++i) {
        if (flags[i] == AMP::Mesh::DendroSearch::FoundOnBoundary) {
            d_ActiveSet.push_back(*inactiveSetIterator);
            inactiveSetIterator = d_InactiveSet.erase(inactiveSetIterator);
            hex8_element_t::get_basis_functions_values_on_face(&(tmpSlaveVerticesLocalCoordOnFace[2*i]), masterShapeFunctionsValuesPointer);
            std::advance(masterShapeFunctionsValuesPointer, 4);
            std::copy(&(tmpSlaveVerticesShift[3*i]), &(tmpSlaveVerticesShift[3*i])+3, slaveVerticesShiftIterator);
            std::advance(slaveVerticesShiftIterator, 3);
            std::copy(&(tmpMasterVerticesGlobalIDs[4*i]), &(tmpMasterVerticesGlobalIDs[4*i])+4, masterVerticesGlobalIDsIterator);
            std::advance(masterVerticesGlobalIDsIterator, 4);
            *masterVolumesGlobalIDsIterator = tmpMasterVolumesGlobalIDs[i];
            ++masterVolumesGlobalIDsIterator;
            *masterFacesLocalIndicesIterator = tmpMasterFacesLocalIndices[i];
            ++masterFacesLocalIndicesIterator;
            ++nInactiveSlaveVerticesActivated;
        } else {
            ++inactiveSetIterator;
        } // end if
    } // end for i
    //  std::advance(inactiveSetIterator, nRemovConstraints);
    std::advance(inactiveSetIterator, nActiveSlaveVerticesDeactivated);
    AMP_ASSERT( inactiveSetIterator == d_InactiveSet.end() );
    AMP_ASSERT( slaveVerticesShiftIterator == d_SlaveShift.end() );
    AMP_ASSERT( masterShapeFunctionsValuesPointer == &(d_MasterShapeFunctionsValues[0])+4*(nOlderConstraints+nAddedConstraints) );
    AMP_ASSERT( masterVerticesGlobalIDsIterator == d_MasterVerticesGlobalIDs.end() );
    AMP_ASSERT( masterVolumesGlobalIDsIterator == d_MasterVolumesGlobalIDs.end() );
    AMP_ASSERT( masterFacesLocalIndicesIterator == d_MasterFacesLocalIndices.end() );
    tmpMasterVerticesGlobalIDs.clear();
    tmpSlaveVerticesShift.clear();
    tmpSlaveVerticesLocalCoordOnFace.clear();
    flags.clear();
    tmpMasterVolumesGlobalIDs.clear(); 
    tmpMasterFacesLocalIndices.clear();

    // if (!rank) { std::cout<<"nInactiveSlaveVerticesActivated="<<d_GlobalComm.sumReduce(nInactiveSlaveVerticesActivated)<<"\n"; }

    /** setup for apply */
    // size_t const npes = d_GlobalComm.getSize();
    d_SendCnts.resize(npes);
    std::fill(d_SendCnts.begin(), d_SendCnts.end(), 0);
    for (size_t i = 0; i < 4*nConstraints; ++i) {
        ++d_SendCnts[d_MasterVerticesGlobalIDs[i].owner_rank()]; 
    } // end for i
    d_SendDisps.resize(npes);
    d_SendDisps[0] = 0;
    for (size_t i = 1; i < npes; ++i) {
        d_SendDisps[i] = d_SendDisps[i-1] + d_SendCnts[i-1]; 
    } // end for i
    AMP_ASSERT( d_SendDisps[npes-1] + d_SendCnts[npes-1] == 4 * static_cast<int>(nConstraints) );

    std::vector<int> tmpSendCnts(npes, 0);
    d_MasterVerticesMap.resize(4*nConstraints, nConstraints);
    std::vector<AMP::Mesh::MeshElementID> sendMasterVerticesGlobalIDs(d_SendDisps[npes-1]+d_SendCnts[npes-1], AMP::Mesh::MeshElementID());
    for (size_t i = 0; i < 4*nConstraints; ++i) {
        size_t sendToRank = d_MasterVerticesGlobalIDs[i].owner_rank();
        d_MasterVerticesMap[i] = d_SendDisps[sendToRank] + tmpSendCnts[sendToRank];
        sendMasterVerticesGlobalIDs[d_MasterVerticesMap[i]] = d_MasterVerticesGlobalIDs[i];
        ++tmpSendCnts[sendToRank];
    } // end for i 
    AMP_ASSERT( std::equal(tmpSendCnts.begin(), tmpSendCnts.end(), d_SendCnts.begin()) );
    tmpSendCnts.clear();

    d_RecvCnts.resize(npes);
    d_GlobalComm.allToAll(1, &(d_SendCnts[0]), &(d_RecvCnts[0]));
    d_RecvDisps.resize(npes);
    d_RecvDisps[0] = 0;
    for (size_t i = 1; i < npes; ++i) {
        d_RecvDisps[i] = d_RecvDisps[i-1] + d_RecvCnts[i-1];
    } // end for i
    d_RecvMasterVerticesGlobalIDs.resize(d_RecvDisps[npes-1]+d_RecvCnts[npes-1]);
    d_GlobalComm.allToAll((!(sendMasterVerticesGlobalIDs.empty()) ? &(sendMasterVerticesGlobalIDs[0]) : NULL), &(d_SendCnts[0]), &(d_SendDisps[0]),
      (!(d_RecvMasterVerticesGlobalIDs.empty()) ? &(d_RecvMasterVerticesGlobalIDs[0]) : NULL), &(d_RecvCnts[0]), &(d_RecvDisps[0]), true);

    d_TransposeSendCnts.resize(npes);
    d_TransposeSendDisps.resize(npes);
    d_TransposeRecvCnts.resize(npes);
    d_TransposeRecvDisps.resize(npes); 
    std::copy(d_SendCnts.begin(), d_SendCnts.end(), d_TransposeSendCnts.begin());
    std::copy(d_SendDisps.begin(), d_SendDisps.end(), d_TransposeSendDisps.begin());
    std::copy(d_RecvCnts.begin(), d_RecvCnts.end(), d_TransposeRecvCnts.begin());
    std::copy(d_RecvDisps.begin(), d_RecvDisps.end(), d_TransposeRecvDisps.begin());

    std::swap_ranges(d_RecvCnts.begin(), d_RecvCnts.end(), d_SendCnts.begin());
    std::swap_ranges(d_RecvDisps.begin(), d_RecvDisps.end(), d_SendDisps.begin());

    getVectorIndicesFromGlobalIDs(d_ActiveSet, d_SlaveIndices);
    getVectorIndicesFromGlobalIDs(d_RecvMasterVerticesGlobalIDs, d_RecvMasterIndices);

    for (size_t i = 0; i < npes; ++i) {
        d_SendCnts[i] *= d_DOFsPerNode; 
        d_SendDisps[i] *= d_DOFsPerNode; 
        d_RecvCnts[i] *= d_DOFsPerNode; 
        d_RecvDisps[i] *= d_DOFsPerNode; 
        d_TransposeSendCnts[i] *= d_DOFsPerNode; 
        d_TransposeSendDisps[i] *= d_DOFsPerNode; 
        d_TransposeRecvCnts[i] *= d_DOFsPerNode; 
        d_TransposeRecvDisps[i] *= d_DOFsPerNode; 
    } // end for i

    d_fout<<std::setprecision(15);
    for (size_t i = 0; i < d_ActiveSet.size(); ++i) {
        for (size_t j = 0; j < 4; ++j) {
            d_fout<<"i="<<i<<"  "
                <<"j="<<j<<"  "
                <<"4*i+j="<<4*i+j<<"  "
                <<"d_MasterShapeFunctionsValues[4*i+j]="<<d_MasterShapeFunctionsValues[4*i+j]<<"\n";
        } // end for j
        for (size_t k = 0; k < d_DOFsPerNode; ++k) {
            d_fout<<"i="<<i<<"  "
                <<"k="<<k<<"  "
                <<"d_DOFsPerNode*i+k="<<d_DOFsPerNode*i+k<<"  "
                <<"d_SlaveIndices[d_DOFsPerNode*i+k]="<<d_SlaveIndices[d_DOFsPerNode*i+k]<<"  "
                <<"d_SlaveVerticesShift[d_DOFsPerNode*i+k]="<<d_SlaveShift[d_DOFsPerNode*i+k]<<"\n";     
        } // end for k
    } // end for i

    d_SlaveVerticesGlobalIDs = d_ActiveSet;
    nInactiveSlaveVerticesActivated = nAddedConstraints;

    size_t const activeSetSizeFinal = d_GlobalComm.sumReduce(d_ActiveSet.size());
    size_t const inactiveSetSizeFinal = d_GlobalComm.sumReduce(d_InactiveSet.size());
    size_t const nInactiveSlaveVerticesActivatedReduced = d_GlobalComm.sumReduce(nInactiveSlaveVerticesActivated);
    size_t const nActiveSlaveVerticesDeactivatedReduced = d_GlobalComm.sumReduce(nActiveSlaveVerticesDeactivated);
    if (!rank) {
        std::cout<<"nInactiveSlaveVerticesActivated="<<nInactiveSlaveVerticesActivatedReduced<<"  "
            <<"nActiveSlaveVerticesDeactivated="<<nActiveSlaveVerticesDeactivatedReduced<<"\n";
        std::cout<<"d_ActiveSet.size()="<<activeSetSizeFinal<<"  "
            <<"d_InactiveSet.size()="<<inactiveSetSizeFinal<<"\n";
    } // end if

    /** move the mesh back */
    if (!skipDisplaceMesh) {
        displacementFieldVector->scale(-1.0);
        d_Mesh->displaceMesh(displacementFieldVector);
        displacementFieldVector->scale(-1.0);
    } // end if

    return d_GlobalComm.sumReduce(nInactiveSlaveVerticesActivated + nActiveSlaveVerticesDeactivated);
}


void NodeToFaceContactOperator::reset(const boost::shared_ptr<OperatorParameters> & params) 
{
    AMP_INSIST( (params != NULL), "NULL parameter" );
    AMP_INSIST( ((params->d_db) != NULL), "NULL database" );

    AMP::Mesh::Mesh::shared_ptr mesh = params->d_Mesh;
    //  AMP::AMP_MPI comm = mesh->getComm(); 
    AMP::AMP_MPI comm = d_GlobalComm;

    /** get the boundary slave vertices coordinates and global IDs */
    size_t nSlaveVertices = 0;
    std::vector<double> tmpSlaveVerticesCoord;
    std::vector<AMP::Mesh::MeshElementID> tmpSlaveVerticesGlobalIDs;
    AMP::Mesh::Mesh::shared_ptr slaveMesh = mesh->Subset(d_SlaveMeshID);
    if (slaveMesh != NULL) {
        AMP::Mesh::MeshIterator slaveMeshIterator = slaveMesh->getBoundaryIDIterator(AMP::Mesh::Vertex, d_SlaveBoundaryID);
        AMP::Mesh::MeshIterator slaveMeshIterator_begin = slaveMeshIterator.begin(), 
          slaveMeshIterator_end = slaveMeshIterator.end();
        nSlaveVertices = slaveMeshIterator.size();
        tmpSlaveVerticesCoord.resize(3*nSlaveVertices);
        tmpSlaveVerticesGlobalIDs.resize(nSlaveVertices);
        std::vector<AMP::Mesh::MeshElementID>::iterator tmpSlaveVerticesGlobalIDsIterator = tmpSlaveVerticesGlobalIDs.begin();
        std::vector<double> tmpCoord(3);
        std::vector<double>::iterator tmpSlaveVerticesCoordIterator = tmpSlaveVerticesCoord.begin();
        for (slaveMeshIterator = slaveMeshIterator_begin; slaveMeshIterator != slaveMeshIterator_end; ++slaveMeshIterator) {
            *tmpSlaveVerticesGlobalIDsIterator = slaveMeshIterator->globalID();
            ++tmpSlaveVerticesGlobalIDsIterator;
            tmpCoord = slaveMeshIterator->coord();
            AMP_ASSERT( tmpCoord.size() == 3 );
            std::copy(tmpCoord.begin(), tmpCoord.end(), tmpSlaveVerticesCoordIterator);
            for (size_t i = 0; i < 3; ++i) { ++tmpSlaveVerticesCoordIterator; }
        } // end loop over the slave vertices on boundary
        AMP_ASSERT( tmpSlaveVerticesGlobalIDsIterator == tmpSlaveVerticesGlobalIDs.end() );
        AMP_ASSERT( tmpSlaveVerticesCoordIterator == tmpSlaveVerticesCoord.end() );
    } //end if

    /** do a dendro search for the boundary slave vertices on the master mesh */
    AMP::Mesh::Mesh::shared_ptr masterMesh = mesh->Subset(d_MasterMeshID);
    AMP::Mesh::DendroSearch dendroSearchOnMaster(masterMesh, false);
    dendroSearchOnMaster.setTolerance(1.0e-10);
    dendroSearchOnMaster.search(comm, tmpSlaveVerticesCoord);

    std::vector<AMP::Mesh::MeshElementID> tmpMasterVerticesGlobalIDs;
    std::vector<double> tmpSlaveVerticesShift, tmpSlaveVerticesLocalCoordOnFace;
    std::vector<int> flags;

    dendroSearchOnMaster.projectOnBoundaryID(comm, d_MasterBoundaryID,
        tmpMasterVerticesGlobalIDs, tmpSlaveVerticesShift, tmpSlaveVerticesLocalCoordOnFace, flags);

    AMP_ASSERT( nSlaveVertices == tmpMasterVerticesGlobalIDs.size() / 4 );
    AMP_ASSERT( nSlaveVertices == tmpSlaveVerticesShift.size() / 3 );
    AMP_ASSERT( nSlaveVertices == tmpSlaveVerticesLocalCoordOnFace.size() / 2 );
    AMP_ASSERT( nSlaveVertices == flags.size() );

    /** build the constraints */
    const unsigned int nConstraints = std::count(flags.begin(), flags.end(), AMP::Mesh::DendroSearch::FoundOnBoundary);

    unsigned int localPtsNotFound = std::count(flags.begin(), flags.end(), AMP::Mesh::DendroSearch::NotFound);
    unsigned int localPtsFoundNotOnBoundary = std::count(flags.begin(), flags.end(), AMP::Mesh::DendroSearch::FoundNotOnBoundary);
    unsigned int localPtsFoundOnBoundary = std::count(flags.begin(), flags.end(), AMP::Mesh::DendroSearch::FoundOnBoundary);
    unsigned int globalPtsNotFound = comm.sumReduce(localPtsNotFound);
    unsigned int globalPtsFoundNotOnBoundary = comm.sumReduce(localPtsFoundNotOnBoundary);
    unsigned int globalPtsFoundOnBoundary = comm.sumReduce(localPtsFoundOnBoundary);
    d_fout<<"Global number of points not found is "<<globalPtsNotFound<<" (local was "<<localPtsNotFound<<")"<<std::endl;
    d_fout<<"Global number of points found not on boundary is "<<globalPtsFoundNotOnBoundary<<" (local was "<<localPtsFoundNotOnBoundary<<")"<<std::endl;
    d_fout<<"Global number of points found on boundary is "<<globalPtsFoundOnBoundary<<" (local was "<<localPtsFoundOnBoundary<<")"<<std::endl;
    d_fout<<"Total number of points is "<<globalPtsNotFound+globalPtsFoundNotOnBoundary+globalPtsFoundOnBoundary<<std::endl;

    AMP_ASSERT( std::count(flags.begin(), flags.end(), AMP::Mesh::DendroSearch::FoundNotOnBoundary) == 0 ); // DendroSearch::FoundNotOnBoundary is not acceptable

    d_SlaveVerticesGlobalIDs.resize(nConstraints);
    std::fill(d_SlaveVerticesGlobalIDs.begin(), d_SlaveVerticesGlobalIDs.end(), AMP::Mesh::MeshElementID());
    std::vector<AMP::Mesh::MeshElementID> masterVerticesGlobalIDs(4*nConstraints);
    std::fill(masterVerticesGlobalIDs.begin(), masterVerticesGlobalIDs.end(), AMP::Mesh::MeshElementID());
    d_MasterVerticesOwnerRanks.resize(4*nConstraints);
    std::fill(d_MasterVerticesOwnerRanks.begin(), d_MasterVerticesOwnerRanks.end(), comm.getSize());
    d_MasterShapeFunctionsValues.resize(4*nConstraints);
    std::fill(d_MasterShapeFunctionsValues.begin(), d_MasterShapeFunctionsValues.end(), 0.0);
    d_SlaveVerticesShift.resize(3*nConstraints);
    std::fill(d_SlaveVerticesShift.begin(), d_SlaveVerticesShift.end(), 0.0);

    std::vector<AMP::Mesh::MeshElementID>::const_iterator tmpSlaveVerticesGlobalIDsConstIterator = tmpSlaveVerticesGlobalIDs.begin();
    std::vector<AMP::Mesh::MeshElementID>::const_iterator tmpMasterVerticesGlobalIDsConstIterator = tmpMasterVerticesGlobalIDs.begin();
    double const * tmpSlaveVerticesLocalCoordOnFacePointerToConst = &(tmpSlaveVerticesLocalCoordOnFace[0]);
    std::vector<double>::const_iterator tmpSlaveVerticesShiftConstIterator = tmpSlaveVerticesShift.begin();
    double * masterShapeFunctionsValuesPointer = &(d_MasterShapeFunctionsValues[0]);
    std::vector<double>::iterator slaveVerticesShiftIterator = d_SlaveVerticesShift.begin();
    std::vector<AMP::Mesh::MeshElementID>::iterator slaveVerticesGlobalIDsIterator = d_SlaveVerticesGlobalIDs.begin();
    std::vector<AMP::Mesh::MeshElementID>::iterator masterVerticesGlobalIDsIterator = masterVerticesGlobalIDs.begin();
    std::vector<size_t>::iterator masterVerticesOwnerRanksIterator = d_MasterVerticesOwnerRanks.begin();

    //  double basis_functions_values_on_face[4];
    std::vector<int>::const_iterator flagsIterator = flags.begin(),
    flagsIterator_end = flags.end();
    for ( ; flagsIterator != flagsIterator_end; ++flagsIterator) {
        // TODO: the following if statement is debug only
        if (*flagsIterator == AMP::Mesh::DendroSearch::NotFound) {
            std::vector<double> blackSheepCoord = (mesh->getElement(*tmpSlaveVerticesGlobalIDsConstIterator)).coord();
            d_fout<<blackSheepCoord[0]<<"  "
                <<blackSheepCoord[1]<<"  "
                <<blackSheepCoord[2]<<"\n";
        } // end if
        //    AMP_ASSERT( (*flagsIterator == AMP::Mesh::DendroSearch::NotFound) || (*flagsIterator == AMP::Mesh::DendroSearch::FoundOnBoundary) );
        if (*flagsIterator == AMP::Mesh::DendroSearch::FoundOnBoundary) {
            hex8_element_t::get_basis_functions_values_on_face(tmpSlaveVerticesLocalCoordOnFacePointerToConst, masterShapeFunctionsValuesPointer);
            for (size_t d = 0; d < 2; ++d) { ++tmpSlaveVerticesLocalCoordOnFacePointerToConst; }
            for (size_t v = 0; v < 4; ++v) { ++masterShapeFunctionsValuesPointer; }
            *slaveVerticesGlobalIDsIterator = *tmpSlaveVerticesGlobalIDsConstIterator;
            ++slaveVerticesGlobalIDsIterator;
            ++tmpSlaveVerticesGlobalIDsConstIterator;
            for (size_t d = 0; d < 3; ++d) { 
                *slaveVerticesShiftIterator = *tmpSlaveVerticesShiftConstIterator;
                ++slaveVerticesShiftIterator;
                ++tmpSlaveVerticesShiftConstIterator;
            } // end for d
            for (size_t v = 0; v < 4; ++v) {
                *masterVerticesGlobalIDsIterator = *tmpMasterVerticesGlobalIDsConstIterator;
                *masterVerticesOwnerRanksIterator = masterVerticesGlobalIDsIterator->owner_rank();
                ++masterVerticesGlobalIDsIterator;
                ++tmpMasterVerticesGlobalIDsConstIterator;
                ++masterVerticesOwnerRanksIterator;
            } // end for v
        } else {
            for (size_t d = 0; d < 2; ++d) { ++tmpSlaveVerticesLocalCoordOnFacePointerToConst; }
            ++tmpSlaveVerticesGlobalIDsConstIterator;
            for (size_t d = 0; d < 3; ++d) { ++tmpSlaveVerticesShiftConstIterator; }
            for (size_t v = 0; v < 4; ++v) { ++tmpMasterVerticesGlobalIDsConstIterator; }
        } // end if
    } // end for
    AMP_ASSERT( tmpSlaveVerticesLocalCoordOnFacePointerToConst == &(tmpSlaveVerticesLocalCoordOnFace[0])+2*nSlaveVertices );
    AMP_ASSERT( tmpSlaveVerticesShiftConstIterator == tmpSlaveVerticesShift.end() );
    AMP_ASSERT( tmpSlaveVerticesGlobalIDsConstIterator == tmpSlaveVerticesGlobalIDs.end() );
    AMP_ASSERT( tmpMasterVerticesGlobalIDsConstIterator == tmpMasterVerticesGlobalIDs.end() );
    AMP_ASSERT( slaveVerticesShiftIterator == d_SlaveVerticesShift.end() );
    AMP_ASSERT( slaveVerticesGlobalIDsIterator == d_SlaveVerticesGlobalIDs.end() );
    AMP_ASSERT( masterShapeFunctionsValuesPointer == &(d_MasterShapeFunctionsValues[0])+4*nConstraints );
    AMP_ASSERT( masterVerticesGlobalIDsIterator == masterVerticesGlobalIDs.end() );
    AMP_ASSERT( masterVerticesOwnerRanksIterator == d_MasterVerticesOwnerRanks.end() );
    tmpSlaveVerticesGlobalIDs.clear();
    tmpMasterVerticesGlobalIDs.clear();
    tmpSlaveVerticesShift.clear();
    tmpSlaveVerticesLocalCoordOnFace.clear();
    flags.clear();

    /** setup for apply */
    size_t npes = comm.getSize();
    d_SendCnts.resize(npes);
    std::fill(d_SendCnts.begin(), d_SendCnts.end(), 0);
    for (size_t i = 0; i < 4*nConstraints; ++i) {
        ++d_SendCnts[d_MasterVerticesOwnerRanks[i]]; 
    } // end for i
    d_SendDisps.resize(npes);
    d_SendDisps[0] = 0;
    for (size_t i = 1; i < npes; ++i) {
        d_SendDisps[i] = d_SendDisps[i-1] + d_SendCnts[i-1]; 
    } // end for i
    AMP_ASSERT( d_SendDisps[npes-1] + d_SendCnts[npes-1] == 4* static_cast<int>(nConstraints) );

    std::vector<int> tmpSendCnts(npes, 0);
    d_MasterVerticesMap.resize(4*nConstraints, nConstraints);
    std::vector<AMP::Mesh::MeshElementID> sendMasterVerticesGlobalIDs(d_SendDisps[npes-1]+d_SendCnts[npes-1], AMP::Mesh::MeshElementID());
    for (size_t i = 0; i < 4*nConstraints; ++i) {
        size_t sendToRank = d_MasterVerticesOwnerRanks[i];
        d_MasterVerticesMap[i] = d_SendDisps[sendToRank] + tmpSendCnts[sendToRank];
        sendMasterVerticesGlobalIDs[d_MasterVerticesMap[i]] = masterVerticesGlobalIDs[i];
        ++tmpSendCnts[sendToRank];
    } // end for i 
    AMP_ASSERT( std::equal(tmpSendCnts.begin(), tmpSendCnts.end(), d_SendCnts.begin()) );
    tmpSendCnts.clear();

    d_RecvCnts.resize(npes);
    comm.allToAll(1, &(d_SendCnts[0]), &(d_RecvCnts[0]));
    d_RecvDisps.resize(npes);
    d_RecvDisps[0] = 0;
    for (size_t i = 1; i < npes; ++i) {
        d_RecvDisps[i] = d_RecvDisps[i-1] + d_RecvCnts[i-1];
    } // end for i
    d_RecvMasterVerticesGlobalIDs.resize(d_RecvDisps[npes-1]+d_RecvCnts[npes-1]);
    comm.allToAll((!(sendMasterVerticesGlobalIDs.empty()) ? &(sendMasterVerticesGlobalIDs[0]) : NULL), &(d_SendCnts[0]), &(d_SendDisps[0]),
      (!(d_RecvMasterVerticesGlobalIDs.empty()) ? &(d_RecvMasterVerticesGlobalIDs[0]) : NULL), &(d_RecvCnts[0]), &(d_RecvDisps[0]), true);

    d_TransposeSendCnts.resize(npes);
    d_TransposeSendDisps.resize(npes);
    d_TransposeRecvCnts.resize(npes);
    d_TransposeRecvDisps.resize(npes); 
    std::copy(d_SendCnts.begin(), d_SendCnts.end(), d_TransposeSendCnts.begin());
    std::copy(d_SendDisps.begin(), d_SendDisps.end(), d_TransposeSendDisps.begin());
    std::copy(d_RecvCnts.begin(), d_RecvCnts.end(), d_TransposeRecvCnts.begin());
    std::copy(d_RecvDisps.begin(), d_RecvDisps.end(), d_TransposeRecvDisps.begin());

    std::swap_ranges(d_RecvCnts.begin(), d_RecvCnts.end(), d_SendCnts.begin());
    std::swap_ranges(d_RecvDisps.begin(), d_RecvDisps.end(), d_SendDisps.begin());

    getVectorIndicesFromGlobalIDs(d_SlaveVerticesGlobalIDs, d_SlaveIndices);
    getVectorIndicesFromGlobalIDs(d_RecvMasterVerticesGlobalIDs, d_RecvMasterIndices);

    for (size_t i = 0; i < npes; ++i) {
        d_SendCnts[i] *= d_DOFsPerNode; 
        d_SendDisps[i] *= d_DOFsPerNode; 
        d_RecvCnts[i] *= d_DOFsPerNode; 
        d_RecvDisps[i] *= d_DOFsPerNode; 
        d_TransposeSendCnts[i] *= d_DOFsPerNode; 
        d_TransposeSendDisps[i] *= d_DOFsPerNode; 
        d_TransposeRecvCnts[i] *= d_DOFsPerNode; 
        d_TransposeRecvDisps[i] *= d_DOFsPerNode; 
    } // end for i

    d_fout<<std::setprecision(15);
    for (size_t i = 0; i < d_SlaveVerticesGlobalIDs.size(); ++i) {
        for (size_t j = 0; j < 4; ++j) {
            d_fout<<"i="<<i<<"  "
                <<"j="<<j<<"  "
                <<"4*i+j="<<4*i+j<<"  "
                <<"d_MasterShapeFunctionsValues[4*i+j]="<<d_MasterShapeFunctionsValues[4*i+j]<<"\n";
        } // end for j
        for (size_t k = 0; k < d_DOFsPerNode; ++k) {
            d_fout<<"i="<<i<<"  "
                <<"k="<<k<<"  "
                <<"d_DOFsPerNode*i+k="<<d_DOFsPerNode*i+k<<"  "
                <<"d_SlaveIndices[d_DOFsPerNode*i+k]="<<d_SlaveIndices[d_DOFsPerNode*i+k]<<"  "
                <<"d_SlaveVerticesShift[d_DOFsPerNode*i+k]="<<d_SlaveVerticesShift[d_DOFsPerNode*i+k]<<"\n";     
        } // end for k
    } // end for i

    // this is really ugly
    d_SlaveShift = d_SlaveVerticesShift;

}


void NodeToFaceContactOperator::getVectorIndicesFromGlobalIDs(const std::vector<AMP::Mesh::MeshElementID> & globalIDs, 
    std::vector<size_t> & vectorIndices) 
{
    std::vector<size_t> tmpIndices;
    std::vector<AMP::Mesh::MeshElementID>::const_iterator globalIDsConstIterator = globalIDs.begin(), 
        globalIDsConstIterator_end = globalIDs.end();
    vectorIndices.resize(globalIDs.size()*d_DOFsPerNode);
    std::vector<size_t>::iterator vectorIndicesIterator = vectorIndices.begin();
    for ( ; globalIDsConstIterator != globalIDsConstIterator_end; ++globalIDsConstIterator) {
        d_DOFManager->getDOFs(*globalIDsConstIterator, tmpIndices);
        AMP_ASSERT( *globalIDsConstIterator != AMP::Mesh::MeshElementID() );
        AMP_ASSERT( tmpIndices.size() == d_DOFsPerNode );
        std::copy(tmpIndices.begin(), tmpIndices.end(), vectorIndicesIterator);
        for (size_t i = 0; i < d_DOFsPerNode; ++i) { ++vectorIndicesIterator; }
    } // end for
    AMP_ASSERT( vectorIndicesIterator == vectorIndices.end() );
}


void NodeToFaceContactOperator::copyMasterToSlave(AMP::LinearAlgebra::Vector::shared_ptr u) 
{
    /** send and receive the master values */
    AMP::AMP_MPI comm = d_GlobalComm;
    //  AMP::AMP_MPI comm = u->getComm();
    size_t npes = comm.getSize();
    //  size_t rank = comm.getRank();

    std::vector<double> sendMasterValues(d_SendDisps[npes-1]+d_SendCnts[npes-1]);
    for (size_t i = 0; i < npes; ++i) {
        for (int j = 0; j < d_SendCnts[i]; j += d_DOFsPerNode) {
            size_t k = d_SendDisps[i] + j;
            u->getLocalValuesByGlobalID(d_DOFsPerNode, &(d_RecvMasterIndices[k]), &(sendMasterValues[k]));
        } // end for j
    } // end for i

    std::vector<double> recvMasterValues(d_RecvDisps[npes-1]+d_RecvCnts[npes-1]);
    comm.allToAll((!(sendMasterValues.empty()) ? &(sendMasterValues[0]) : NULL), &(d_SendCnts[0]), &(d_SendDisps[0]),
      (!(recvMasterValues.empty()) ? &(recvMasterValues[0]) : NULL), &(d_RecvCnts[0]), &(d_RecvDisps[0]), true);
    sendMasterValues.clear();

    /** compute slave values */
    std::vector<double> slaveValues(d_SlaveIndices.size(), 0.0);

    for (size_t i = 0; i < d_SlaveVerticesGlobalIDs.size(); ++i) {
        for (size_t j = 0; j < 4; ++j) {
          for (size_t k = 0; k < d_DOFsPerNode; ++k) {
            slaveValues[d_DOFsPerNode*i+k] += d_MasterShapeFunctionsValues[4*i+j] * recvMasterValues[d_DOFsPerNode*d_MasterVerticesMap[4*i+j]+k];
          } // end for k
        } // end for j
    } // end for i

    if (!d_SlaveVerticesGlobalIDs.empty()) { 
        u->setLocalValuesByGlobalID(d_SlaveIndices.size(), &(d_SlaveIndices[0]), &(slaveValues[0]));
    } // end if
}


void NodeToFaceContactOperator::addSlaveToMaster(AMP::LinearAlgebra::Vector::shared_ptr u) 
{
    /** send and receive slave value times shape functions values */
    AMP::AMP_MPI comm = d_GlobalComm;
    //  AMP::AMP_MPI comm = r->getComm();
    size_t npes = comm.getSize();

    std::vector<double> sendAddToMasterValues(d_TransposeSendDisps[npes-1]+d_TransposeSendCnts[npes-1]);
    for (size_t i = 0; i < d_SlaveVerticesGlobalIDs.size(); ++i) {
        for (size_t j = 0; j < 4; ++j) {
          u->getLocalValuesByGlobalID(d_DOFsPerNode, &(d_SlaveIndices[d_DOFsPerNode*i]), &(sendAddToMasterValues[d_DOFsPerNode*d_MasterVerticesMap[4*i+j]])); 
          for (size_t k = 0; k < d_DOFsPerNode; ++k) {
            sendAddToMasterValues[d_DOFsPerNode*d_MasterVerticesMap[4*i+j]+k] *= d_MasterShapeFunctionsValues[4*i+j];
          } // end for k
        } // end for j
    } // end for i

    std::vector<double> recvAddToMasterValues(d_TransposeRecvDisps[npes-1]+d_TransposeRecvCnts[npes-1]);
    comm.allToAll((!(sendAddToMasterValues.empty()) ? &(sendAddToMasterValues[0]) : NULL), &(d_TransposeSendCnts[0]), &(d_TransposeSendDisps[0]),
          (!(recvAddToMasterValues.empty()) ? &(recvAddToMasterValues[0]) : NULL), &(d_TransposeRecvCnts[0]), &(d_TransposeRecvDisps[0]), true);
    sendAddToMasterValues.clear();

    /** add slave value times shape functions values to master values and set slave values to zero */
    if (!d_RecvMasterVerticesGlobalIDs.empty()) {
        u->addLocalValuesByGlobalID(d_RecvMasterIndices.size(), &(d_RecvMasterIndices[0]), &(recvAddToMasterValues[0]));
    } // end if
}


} // end namespace Operator
} // end namespace AMP


