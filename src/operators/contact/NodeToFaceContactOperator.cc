
#include "operators/mechanics/IsotropicElasticModel.h"
#include "ampmesh/MeshID.h"
#include "ampmesh/dendro/DendroSearch.h"
#include "ampmesh/euclidean_geometry_tools.h"
#include "ampmesh/latex_visualization_tools.h"
#include "operators/contact/NodeToFaceContactOperator.h"

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iterator>
#include <numeric>

namespace AMP {
namespace Operator {

void NodeToFaceContactOperator::initialize()
{
    d_ActiveSet.clear();
    d_InactiveSet.clear();
    /** get all slave boundary vertices and tag them as inactive */
    AMP::Mesh::Mesh::shared_ptr slaveMesh = d_Mesh->Subset( d_SlaveMeshID );
    if ( slaveMesh.get() != NULL ) {
        AMP::Mesh::MeshIterator slaveMeshIterator =
            slaveMesh->getBoundaryIDIterator( AMP::Mesh::Vertex, d_SlaveBoundaryID );
        AMP::Mesh::MeshIterator slaveMeshIterator_begin = slaveMeshIterator.begin(),
                                slaveMeshIterator_end   = slaveMeshIterator.end();
        size_t const nSlaveVertices                     = slaveMeshIterator.size();
        d_ActiveSet.reserve( nSlaveVertices );
        d_InactiveSet.reserve( nSlaveVertices );
        for ( slaveMeshIterator = slaveMeshIterator_begin;
              slaveMeshIterator != slaveMeshIterator_end;
              ++slaveMeshIterator ) {
            d_InactiveSet.push_back( slaveMeshIterator->globalID() );
        } // end loop over the slave vertices on boundary
    }     // end if
}

size_t NodeToFaceContactOperator::updateActiveSet(
    AMP::LinearAlgebra::Vector::shared_ptr displacementFieldVector, bool skipDisplaceMesh )
{
    size_t const npes = d_GlobalComm.getSize();
    size_t const rank = d_GlobalComm.getRank();

    size_t nActiveSlaveVerticesDeactivated          = 0;
    size_t nInactiveSlaveVerticesActivated          = 0;
    size_t const nInactiveSlaveVerticesBeforeUpdate = d_InactiveSet.size();
    size_t const nActiveSlaveVerticesBeforeUpdate   = d_ActiveSet.size();

    d_SlaveVerticesNormalVectorBeforeUpdate.resize( 3 * nActiveSlaveVerticesBeforeUpdate );
    d_SlaveVerticesSurfaceTractionBeforeUpdate.resize( 3 * nActiveSlaveVerticesBeforeUpdate );

    /** displace the mesh */
    if ( !skipDisplaceMesh ) {
        d_Mesh->displaceMesh( displacementFieldVector );
    } // end if

    /** get the inactive slave vertices coordinates */
    std::vector<double> inactiveSlaveVerticesCoord( 3 * nInactiveSlaveVerticesBeforeUpdate );
    AMP::Mesh::MeshElement tmpInactiveSlaveVertex;
    std::vector<double> tmpInactiveSlaveVertexCoord( 3 );
    for ( size_t i = 0; i < nInactiveSlaveVerticesBeforeUpdate; ++i ) {
        tmpInactiveSlaveVertex      = d_Mesh->getElement( d_InactiveSet[i] );
        tmpInactiveSlaveVertexCoord = tmpInactiveSlaveVertex.coord();
        std::copy( tmpInactiveSlaveVertexCoord.begin(),
                   tmpInactiveSlaveVertexCoord.end(),
                   &( inactiveSlaveVerticesCoord[3 * i] ) );
    } // end for i

    /** perform a dendro search for all inactive slave vertices over the master mesh */
    AMP::Mesh::Mesh::shared_ptr masterMesh = d_Mesh->Subset( d_MasterMeshID );
    // TODO: read dummyVerboseFlag from input file
    bool const dummyVerboseFlag = false;
    AMP::Mesh::DendroSearch dendroSearchOnMaster( masterMesh, dummyVerboseFlag );
    dendroSearchOnMaster.setTolerance( 1.0e-10 );
    dendroSearchOnMaster.search( d_GlobalComm, inactiveSlaveVerticesCoord );

    std::vector<AMP::Mesh::MeshElementID> tmpMasterVerticesGlobalIDs;
    std::vector<double> tmpSlaveVerticesShift, tmpSlaveVerticesLocalCoordOnFace;
    std::vector<int> flags;
    std::vector<AMP::Mesh::MeshElementID> tmpMasterVolumesGlobalIDs;
    std::vector<size_t> tmpMasterFacesLocalIndices;
    dendroSearchOnMaster.projectOnBoundaryID( d_GlobalComm,
                                              d_MasterBoundaryID,
                                              tmpMasterVerticesGlobalIDs,
                                              tmpSlaveVerticesShift,
                                              tmpSlaveVerticesLocalCoordOnFace,
                                              flags,
                                              tmpMasterVolumesGlobalIDs,
                                              tmpMasterFacesLocalIndices );

    std::vector<int> rankMap( npes );
    int myRank = -1;
    if ( masterMesh.get() != NULL ) {
        AMP::AMP_MPI meshComm = masterMesh->getComm();
        myRank                = meshComm.getRank();
    }

    d_GlobalComm.allGather( myRank, &( rankMap[0] ) );

    std::vector<int> invRankMap( npes, -1 );
    for ( unsigned int i = 0; i < npes; ++i ) {
        if ( rankMap[i] >= 0 ) {
            invRankMap[rankMap[i]] = i;
        }
    } // end i


    /** move the mesh back */
    if ( !skipDisplaceMesh ) {
        displacementFieldVector->scale( -1.0 );
        d_Mesh->displaceMesh( displacementFieldVector );
        displacementFieldVector->scale( -1.0 );
    } // end if

    if ( d_GlobalComm.anyReduce( !d_ActiveSet.empty() ) ) {
        /** compute the normal vectors at slave nodes */
        std::vector<int> sendCnts( npes, 0 );
        for ( size_t i = 0; i < nActiveSlaveVerticesBeforeUpdate; ++i ) {
            //          ++sendCnts[d_MasterVolumesGlobalIDs[i].owner_rank()];
            ++sendCnts[invRankMap[d_MasterVolumesGlobalIDs[i].owner_rank()]];
        } // end for i
        std::vector<int> sendDisps( npes, 0 );
        sendDisps[0] = 0;
        for ( size_t i = 1; i < npes; ++i ) {
            sendDisps[i] = sendDisps[i - 1] + sendCnts[i - 1];
        } // end for i
        size_t const nSendData = static_cast<size_t>( sendDisps[npes - 1] + sendCnts[npes - 1] );
        AMP_ASSERT( nSendData == nActiveSlaveVerticesBeforeUpdate );
        std::vector<ProjectionData> sendProjectionDataBuffer( nSendData );
        std::fill( sendCnts.begin(), sendCnts.end(), 0 );
        std::vector<size_t> sendMap( nSendData, nSendData );
        std::vector<size_t> recvMap( nSendData, nSendData );
        for ( size_t i = 0; i < nSendData; ++i ) {
            //          size_t sendToRank = d_MasterVolumesGlobalIDs[i].owner_rank();
            size_t sendToRank = invRankMap[d_MasterVolumesGlobalIDs[i].owner_rank()];
            sendMap[i]        = static_cast<size_t>( sendDisps[sendToRank] + sendCnts[sendToRank] );
            recvMap[sendMap[i]] = i;
            // AMP_ASSERT( i  == sendMap[i] && sendMap[i] == recvMap[i] );
            sendProjectionDataBuffer[sendMap[i]].d_MasterVolumeGlobalID =
                d_MasterVolumesGlobalIDs[i];
            sendProjectionDataBuffer[sendMap[i]].d_MasterFaceLocalIndex =
                d_MasterFacesLocalIndices[i];
            hex8_element_t::get_local_coordinates_on_face(
                &( d_MasterShapeFunctionsValues[4 * i] ),
                sendProjectionDataBuffer[sendMap[i]].d_SlaveVertexLocalCoordOnMasterFace );
            ++sendCnts[sendToRank];
        } // end for i
        AMP_ASSERT( std::find( sendMap.begin(), sendMap.end(), nSendData ) == sendMap.end() );
        AMP_ASSERT( std::find( recvMap.begin(), recvMap.end(), nSendData ) == recvMap.end() );
        std::vector<int> recvCnts( npes, 0 );
        d_GlobalComm.allToAll( 1, &( sendCnts[0] ), &( recvCnts[0] ) );
        std::vector<int> recvDisps( npes, 0 );
        recvDisps[0] = 0;
        for ( size_t i = 1; i < npes; ++i ) {
            recvDisps[i] = recvDisps[i - 1] + recvCnts[i - 1];
        } // end for i
        size_t const nRecvData = static_cast<size_t>( recvDisps[npes - 1] + recvCnts[npes - 1] );
        std::vector<ProjectionData> recvProjectionDataBuffer( nRecvData );
        d_GlobalComm.allToAll(
            ( !( sendProjectionDataBuffer.empty() ) ? &( sendProjectionDataBuffer[0] ) : NULL ),
            &( sendCnts[0] ),
            &( sendDisps[0] ),
            ( !( recvProjectionDataBuffer.empty() ) ? &( recvProjectionDataBuffer[0] ) : NULL ),
            &( recvCnts[0] ),
            &( recvDisps[0] ),
            true );
        sendProjectionDataBuffer.clear();

        std::vector<StressStateData> sendStressStateDataBuffer( nRecvData );
        for ( size_t i = 0; i < nRecvData; ++i ) {
            // retrieve master volume element then compute normal vector and strain tensor at slave
            // vertex
            std::vector<AMP::Mesh::MeshElement> masterVolumeVertices =
                d_Mesh->getElement( recvProjectionDataBuffer[i].d_MasterVolumeGlobalID )
                    .getElements( AMP::Mesh::Vertex );
            AMP_ASSERT( masterVolumeVertices.size() == 8 );
            std::vector<AMP::Mesh::MeshElementID> masterVolumeVerticesGlobalIDs( 8 );
            double masterVolumeVerticesCoordinates[24];
            std::vector<double> vertexCoord( 3 );
            for ( size_t j = 0; j < masterVolumeVertices.size(); ++j ) {
                vertexCoord = masterVolumeVertices[j].coord();
                std::copy( vertexCoord.begin(),
                           vertexCoord.end(),
                           &( masterVolumeVerticesCoordinates[3 * j] ) );
                masterVolumeVerticesGlobalIDs[j] = masterVolumeVertices[j].globalID();
            } // end for j
            hex8_element_t masterVolumeElement( masterVolumeVerticesCoordinates );
            masterVolumeElement.compute_normal_to_face(
                recvProjectionDataBuffer[i].d_MasterFaceLocalIndex,
                recvProjectionDataBuffer[i].d_SlaveVertexLocalCoordOnMasterFace,
                sendStressStateDataBuffer[i].d_SlaveVertexNormalVector );
            double slaveVertexLocalCoordinates[3];
            hex8_element_t::map_face_to_local(
                recvProjectionDataBuffer[i].d_MasterFaceLocalIndex,
                recvProjectionDataBuffer[i].d_SlaveVertexLocalCoordOnMasterFace,
                slaveVertexLocalCoordinates );
            std::vector<size_t> displacementIndices;
            getVectorIndicesFromGlobalIDs( masterVolumeVerticesGlobalIDs, displacementIndices );
            AMP_ASSERT( displacementIndices.size() == 24 );
            double displacementValues[24];
            displacementFieldVector->getValuesByGlobalID(
                24, &( displacementIndices[0] ), &( displacementValues[0] ) );
            double strainTensor[6];
            masterVolumeElement.compute_strain_tensor(
                slaveVertexLocalCoordinates, displacementValues, strainTensor );
            // compute thermal strain
            if ( d_TemperatureFieldVector.get() != NULL ) {
                std::vector<size_t> temperatureIndices;
                d_TemperatureDOFManager->getDOFs( masterVolumeVerticesGlobalIDs,
                                                  temperatureIndices );
                AMP_ASSERT( temperatureIndices.size() == 8 );
                double temperatureValues[8];
                d_TemperatureFieldVector->getValuesByGlobalID(
                    8, &( temperatureIndices[0] ), &( temperatureValues[0] ) );
                double basisFunctionsValues[8];
                hex8_element_t::get_basis_functions_values( slaveVertexLocalCoordinates,
                                                            basisFunctionsValues );
                double temperature = std::inner_product(
                    temperatureValues, temperatureValues + 8, basisFunctionsValues, 0.0 );
                std::transform( strainTensor,
                                strainTensor + 3,
                                strainTensor,
                                std::bind2nd( std::minus<double>(),
                                              d_ThermalExpansionCoefficient *
                                                  ( temperature - d_ReferenceTemperature ) ) );
            } // end if

            // compute surface traction at slave vertex
            double youngsModulus = AMP::dynamic_pointer_cast<AMP::Operator::IsotropicElasticModel>(
                                       d_MasterMechanicsMaterialModel )
                                       ->getYoungsModulus();
            double poissonsRatio = AMP::dynamic_pointer_cast<AMP::Operator::IsotropicElasticModel>(
                                       d_MasterMechanicsMaterialModel )
                                       ->getPoissonsRatio();
            double constitutiveMatrix[36];
            compute_constitutive_matrix( youngsModulus, poissonsRatio, constitutiveMatrix );
            double stressTensor[6];
            compute_stress_tensor( constitutiveMatrix, strainTensor, stressTensor );
            compute_traction( stressTensor,
                              sendStressStateDataBuffer[i].d_SlaveVertexNormalVector,
                              sendStressStateDataBuffer[i].d_SlaveVertexSurfaceTraction );
        } // end for i
        recvProjectionDataBuffer.clear();

        sendCnts.swap( recvCnts );
        sendDisps.swap( recvDisps );
        std::vector<StressStateData> recvStressStateDataBuffer( nSendData );
        d_GlobalComm.allToAll(
            ( !( sendStressStateDataBuffer.empty() ) ? &( sendStressStateDataBuffer[0] ) : NULL ),
            &( sendCnts[0] ),
            &( sendDisps[0] ),
            ( !( recvStressStateDataBuffer.empty() ) ? &( recvStressStateDataBuffer[0] ) : NULL ),
            &( recvCnts[0] ),
            &( recvDisps[0] ),
            true );
        sendStressStateDataBuffer.clear();

        for ( size_t i = 0; i < nSendData; ++i ) {
            std::transform( recvStressStateDataBuffer[i].d_SlaveVertexNormalVector,
                            recvStressStateDataBuffer[i].d_SlaveVertexNormalVector + 3,
                            &( d_SlaveVerticesNormalVectorBeforeUpdate[3 * recvMap[i]] ),
                            std::bind1st( std::multiplies<double>(), -1.0 ) );
            std::transform( recvStressStateDataBuffer[i].d_SlaveVertexSurfaceTraction,
                            recvStressStateDataBuffer[i].d_SlaveVertexSurfaceTraction + 3,
                            &( d_SlaveVerticesSurfaceTractionBeforeUpdate[3 * recvMap[i]] ),
                            std::bind1st( std::multiplies<double>(), -1.0 ) );
        } // end for i

        recvStressStateDataBuffer.clear();

        std::vector<double>::iterator slaveVerticesShiftIterator = d_SlaveShift.begin();
        std::vector<AMP::Mesh::MeshElementID>::iterator masterVerticesGlobalIDsIterator =
            d_MasterVerticesGlobalIDs.begin();
        std::vector<double>::iterator masterShapeFunctionsValuesIterator =
            d_MasterShapeFunctionsValues.begin();
        std::vector<AMP::Mesh::MeshElementID>::iterator masterVolumesGlobalIDsIterator =
            d_MasterVolumesGlobalIDs.begin();
        std::vector<size_t>::iterator masterFacesLocalIndicesIterator =
            d_MasterFacesLocalIndices.begin();
        std::vector<double>::iterator slaveVerticesNormalVectorIterator =
            d_SlaveVerticesNormalVector.begin();

        std::vector<AMP::Mesh::MeshElementID>::iterator activeSetIterator = d_ActiveSet.begin();
        for ( size_t i = 0; i < nActiveSlaveVerticesBeforeUpdate; ++i ) {
            double nDotT =
                compute_scalar_product( &( d_SlaveVerticesNormalVectorBeforeUpdate[3 * i] ),
                                        &( d_SlaveVerticesSurfaceTractionBeforeUpdate[3 * i] ) );
            // std::vector<double> slaveVertexCoordinates =
            // d_Mesh->getElement(*activeSetIterator).coord();
            if ( nDotT > 1.0e-10 ) {
                ++nActiveSlaveVerticesDeactivated;
                d_InactiveSet.push_back( *activeSetIterator );
                activeSetIterator          = d_ActiveSet.erase( activeSetIterator );
                slaveVerticesShiftIterator = d_SlaveShift.erase( slaveVerticesShiftIterator,
                                                                 slaveVerticesShiftIterator + 3 );
                masterVerticesGlobalIDsIterator = d_MasterVerticesGlobalIDs.erase(
                    masterVerticesGlobalIDsIterator, masterVerticesGlobalIDsIterator + 4 );
                masterShapeFunctionsValuesIterator = d_MasterShapeFunctionsValues.erase(
                    masterShapeFunctionsValuesIterator, masterShapeFunctionsValuesIterator + 4 );
                masterVolumesGlobalIDsIterator =
                    d_MasterVolumesGlobalIDs.erase( masterVolumesGlobalIDsIterator );
                masterFacesLocalIndicesIterator =
                    d_MasterFacesLocalIndices.erase( masterFacesLocalIndicesIterator );
                slaveVerticesNormalVectorIterator = d_SlaveVerticesNormalVector.erase(
                    slaveVerticesNormalVectorIterator, slaveVerticesNormalVectorIterator + 3 );
            }
            else {
                d_fout << "\n";
                ++activeSetIterator;
                std::advance( slaveVerticesShiftIterator, 3 );
                std::advance( masterVerticesGlobalIDsIterator, 4 );
                std::advance( masterShapeFunctionsValuesIterator, 4 );
                ++masterVolumesGlobalIDsIterator;
                ++masterFacesLocalIndicesIterator;
                std::advance( slaveVerticesNormalVectorIterator, 3 );
            } // end if
        }     // end for i
        AMP_ASSERT( activeSetIterator == d_ActiveSet.end() );
        AMP_ASSERT( slaveVerticesShiftIterator == d_SlaveShift.end() );
        AMP_ASSERT( masterVerticesGlobalIDsIterator == d_MasterVerticesGlobalIDs.end() );
        AMP_ASSERT( masterShapeFunctionsValuesIterator == d_MasterShapeFunctionsValues.end() );
        AMP_ASSERT( masterVolumesGlobalIDsIterator == d_MasterVolumesGlobalIDs.end() );
        AMP_ASSERT( masterFacesLocalIndicesIterator == d_MasterFacesLocalIndices.end() );
        AMP_ASSERT( slaveVerticesNormalVectorIterator == d_SlaveVerticesNormalVector.end() );

    } // end if

    nInactiveSlaveVerticesActivated = static_cast<size_t>(
        std::count( flags.begin(),
                    flags.end(),
                    AMP::Mesh::DendroSearch::FoundOnBoundary ) ); // std::count returns ptrdiff_t

    ptrdiff_t localPtsNotFound =
        std::count( flags.begin(), flags.end(), AMP::Mesh::DendroSearch::NotFound );
    ptrdiff_t localPtsFoundNotOnBoundary =
        std::count( flags.begin(), flags.end(), AMP::Mesh::DendroSearch::FoundNotOnBoundary );
    ptrdiff_t localPtsFoundOnBoundary =
        std::count( flags.begin(), flags.end(), AMP::Mesh::DendroSearch::FoundOnBoundary );
    ptrdiff_t globalPtsNotFound           = d_GlobalComm.sumReduce( localPtsNotFound );
    ptrdiff_t globalPtsFoundNotOnBoundary = d_GlobalComm.sumReduce( localPtsFoundNotOnBoundary );
    ptrdiff_t globalPtsFoundOnBoundary    = d_GlobalComm.sumReduce( localPtsFoundOnBoundary );
    d_fout << "Global number of points not found is " << globalPtsNotFound << " (local was "
           << localPtsNotFound << ")" << std::endl;
    d_fout << "Global number of points found not on boundary is " << globalPtsFoundNotOnBoundary
           << " (local was " << localPtsFoundNotOnBoundary << ")" << std::endl;
    d_fout << "Global number of points found on boundary is " << globalPtsFoundOnBoundary
           << " (local was " << localPtsFoundOnBoundary << ")" << std::endl;
    d_fout << "Total number of points is "
           << globalPtsNotFound + globalPtsFoundNotOnBoundary + globalPtsFoundOnBoundary
           << std::endl;

    size_t const nActiveSlaveVerticesAfterUpdate = nActiveSlaveVerticesBeforeUpdate -
                                                   nActiveSlaveVerticesDeactivated +
                                                   nInactiveSlaveVerticesActivated;
    size_t const nInactiveSlaveVerticesAfterUpdate = nInactiveSlaveVerticesBeforeUpdate +
                                                     nActiveSlaveVerticesDeactivated -
                                                     nInactiveSlaveVerticesActivated;

    size_t const nActiveSlaveVerticesTmp = d_ActiveSet.size();
    d_SlaveVerticesNormalVector.resize( 3 * nActiveSlaveVerticesAfterUpdate, 0.0 );
    d_SlaveShift.resize( 3 * nActiveSlaveVerticesAfterUpdate, 0.0 );
    d_MasterVerticesGlobalIDs.resize( 4 * nActiveSlaveVerticesAfterUpdate,
                                      AMP::Mesh::MeshElementID() );
    d_MasterShapeFunctionsValues.resize( 4 * nActiveSlaveVerticesAfterUpdate, 0.0 );
    d_MasterVolumesGlobalIDs.resize( nActiveSlaveVerticesAfterUpdate, AMP::Mesh::MeshElementID() );
    d_MasterFacesLocalIndices.resize( nActiveSlaveVerticesAfterUpdate, 6 );

    std::vector<double>::iterator slaveVerticesShiftIterator = d_SlaveShift.begin();
    std::vector<AMP::Mesh::MeshElementID>::iterator masterVerticesGlobalIDsIterator =
        d_MasterVerticesGlobalIDs.begin();
    double *masterShapeFunctionsValuesPointer =
        &( d_MasterShapeFunctionsValues[0] ); // needed this way to work with
    // hex8_element_t::get_basis_functions_values_on_face(...)
    std::vector<AMP::Mesh::MeshElementID>::iterator masterVolumesGlobalIDsIterator =
        d_MasterVolumesGlobalIDs.begin();
    std::vector<size_t>::iterator masterFacesLocalIndicesIterator =
        d_MasterFacesLocalIndices.begin();

    AMP_ASSERT( nActiveSlaveVerticesTmp ==
                nActiveSlaveVerticesBeforeUpdate - nActiveSlaveVerticesDeactivated );
    std::advance( slaveVerticesShiftIterator, 3 * nActiveSlaveVerticesTmp );
    std::advance( masterVerticesGlobalIDsIterator, 4 * nActiveSlaveVerticesTmp );
    std::advance( masterShapeFunctionsValuesPointer, 4 * nActiveSlaveVerticesTmp );
    std::advance( masterVolumesGlobalIDsIterator, nActiveSlaveVerticesTmp );
    std::advance( masterFacesLocalIndicesIterator, nActiveSlaveVerticesTmp );

    std::vector<AMP::Mesh::MeshElementID>::iterator inactiveSetIterator = d_InactiveSet.begin();
    for ( size_t i = 0; i < nInactiveSlaveVerticesBeforeUpdate; ++i ) {
        if ( flags[i] == AMP::Mesh::DendroSearch::FoundOnBoundary ) {
            d_ActiveSet.push_back( *inactiveSetIterator );
            inactiveSetIterator = d_InactiveSet.erase( inactiveSetIterator );
            hex8_element_t::get_basis_functions_values_on_face(
                &( tmpSlaveVerticesLocalCoordOnFace[2 * i] ), masterShapeFunctionsValuesPointer );
            std::advance( masterShapeFunctionsValuesPointer, 4 );
            std::copy( &( tmpSlaveVerticesShift[3 * i] ),
                       &( tmpSlaveVerticesShift[3 * i] ) + 3,
                       slaveVerticesShiftIterator );
            std::advance( slaveVerticesShiftIterator, 3 );
            std::copy( &( tmpMasterVerticesGlobalIDs[4 * i] ),
                       &( tmpMasterVerticesGlobalIDs[4 * i] ) + 4,
                       masterVerticesGlobalIDsIterator );
            std::advance( masterVerticesGlobalIDsIterator, 4 );
            *masterVolumesGlobalIDsIterator = tmpMasterVolumesGlobalIDs[i];
            ++masterVolumesGlobalIDsIterator;
            *masterFacesLocalIndicesIterator = tmpMasterFacesLocalIndices[i];
            ++masterFacesLocalIndicesIterator;
        }
        else {
            if ( flags[i] == AMP::Mesh::DendroSearch::FoundNotOnBoundary ) {
                // check what master elements the slave vertices where found in before throwing an
                // error
                d_fout << std::setprecision( 6 ) << std::fixed;
                std::vector<AMP::Mesh::MeshElement> elementVertices =
                    d_Mesh->getElement( tmpMasterVolumesGlobalIDs[i] )
                        .getElements( AMP::Mesh::Vertex );
                double verticesCoordinates[24];
                for ( size_t kk = 0; kk < 8; ++kk ) {
                    std::vector<double> vertexCoordinates = elementVertices[kk].coord();
                    std::copy( vertexCoordinates.begin(),
                               vertexCoordinates.end(),
                               &( verticesCoordinates[3 * kk] ) );
                } // end for kk
                hex8_element_t volumeElement( verticesCoordinates );
                double point_of_view[3] = { 1.0, 1.0, 1.0 };
                draw_hex8_element( &volumeElement, point_of_view, d_fout );
            } // end if
            d_fout << std::flush;

            ++inactiveSetIterator;
        } // end if
    }     // end for i
    std::advance( inactiveSetIterator, nActiveSlaveVerticesDeactivated );
    AMP_ASSERT( inactiveSetIterator == d_InactiveSet.end() );
    AMP_ASSERT( slaveVerticesShiftIterator == d_SlaveShift.end() );
    AMP_ASSERT( masterShapeFunctionsValuesPointer ==
                &( d_MasterShapeFunctionsValues[0] ) + 4 * nActiveSlaveVerticesAfterUpdate );
    AMP_ASSERT( masterVerticesGlobalIDsIterator == d_MasterVerticesGlobalIDs.end() );
    AMP_ASSERT( masterVolumesGlobalIDsIterator == d_MasterVolumesGlobalIDs.end() );
    AMP_ASSERT( masterFacesLocalIndicesIterator == d_MasterFacesLocalIndices.end() );
    AMP_ASSERT(
        std::count( flags.begin(), flags.end(), AMP::Mesh::DendroSearch::FoundNotOnBoundary ) ==
        0 ); // DendroSearch::FoundNotOnBoundary is not acceptable
    tmpMasterVerticesGlobalIDs.clear();
    tmpSlaveVerticesShift.clear();
    tmpSlaveVerticesLocalCoordOnFace.clear();
    flags.clear();
    tmpMasterVolumesGlobalIDs.clear();
    tmpMasterFacesLocalIndices.clear();

    /** compute slave vertices shift correction */
    std::vector<int> sendCnts( npes, 0 );
    for ( size_t i = nActiveSlaveVerticesTmp; i < nActiveSlaveVerticesAfterUpdate; ++i ) {
        //          ++sendCnts[d_MasterVolumesGlobalIDs[i].owner_rank()];
        ++sendCnts[invRankMap[d_MasterVolumesGlobalIDs[i].owner_rank()]];
    } // end for i
    std::vector<int> sendDisps( npes, 0 );
    sendDisps[0] = 0;
    for ( size_t i = 1; i < npes; ++i ) {
        sendDisps[i] = sendDisps[i - 1] + sendCnts[i - 1];
    } // end for i
    size_t const nSendData = static_cast<size_t>( sendDisps[npes - 1] + sendCnts[npes - 1] );
    AMP_ASSERT( nSendData == nInactiveSlaveVerticesActivated );
    std::vector<ProjectionData> sendProjectionDataBuffer( nSendData );
    std::fill( sendCnts.begin(), sendCnts.end(), 0 );
    std::vector<size_t> sendMap( nSendData, nSendData );
    //        std::vector<size_t> recvMap(nSendData, nSendData);
    for ( size_t i = 0; i < nSendData; ++i ) {
        //          size_t sendToRank =
        //          d_MasterVolumesGlobalIDs[nActiveSlaveVerticesTmp+i].owner_rank();
        size_t sendToRank =
            invRankMap[d_MasterVolumesGlobalIDs[nActiveSlaveVerticesTmp + i].owner_rank()];
        sendMap[i] = static_cast<size_t>( sendDisps[sendToRank] + sendCnts[sendToRank] );
        //          recvMap[sendMap[i]] = i;
        sendProjectionDataBuffer[sendMap[i]].d_MasterVolumeGlobalID =
            d_MasterVolumesGlobalIDs[nActiveSlaveVerticesTmp + i];
        sendProjectionDataBuffer[sendMap[i]].d_MasterFaceLocalIndex =
            d_MasterFacesLocalIndices[nActiveSlaveVerticesTmp + i];
        hex8_element_t::get_local_coordinates_on_face(
            &( d_MasterShapeFunctionsValues[4 * ( nActiveSlaveVerticesTmp + i )] ),
            sendProjectionDataBuffer[sendMap[i]].d_SlaveVertexLocalCoordOnMasterFace );
        ++sendCnts[sendToRank];
    } // end for i
    std::vector<int> recvCnts( npes, 0 );
    d_GlobalComm.allToAll( 1, &( sendCnts[0] ), &( recvCnts[0] ) );
    std::vector<int> recvDisps( npes, 0 );
    recvDisps[0] = 0;
    for ( size_t i = 1; i < npes; ++i ) {
        recvDisps[i] = recvDisps[i - 1] + recvCnts[i - 1];
    } // end for i
    size_t const nRecvData = static_cast<size_t>( recvDisps[npes - 1] + recvCnts[npes - 1] );
    std::vector<ProjectionData> recvProjectionDataBuffer( nRecvData );
    d_GlobalComm.allToAll(
        ( !( sendProjectionDataBuffer.empty() ) ? &( sendProjectionDataBuffer[0] ) : NULL ),
        &( sendCnts[0] ),
        &( sendDisps[0] ),
        ( !( recvProjectionDataBuffer.empty() ) ? &( recvProjectionDataBuffer[0] ) : NULL ),
        &( recvCnts[0] ),
        &( recvDisps[0] ),
        true );
    sendProjectionDataBuffer.clear();

    std::vector<double> sendProjOnMasterFaceDisplacements( 3 * nRecvData );
    std::vector<AnotherDataWithNoName> sendProjOnMasterFaceData( nRecvData );
    for ( size_t i = 0; i < nRecvData; ++i ) {
        std::vector<AMP::Mesh::MeshElement> masterVolumeVertices =
            d_Mesh->getElement( recvProjectionDataBuffer[i].d_MasterVolumeGlobalID )
                .getElements( AMP::Mesh::Vertex );
        std::vector<AMP::Mesh::MeshElementID> masterVolumeVerticesGlobalIDs( 8 );
        AMP_ASSERT( masterVolumeVertices.size() == 8 );
        double masterVolumeVerticesCoordinates[24];
        std::vector<double> vertexCoord( 3 );
        for ( size_t j = 0; j < masterVolumeVertices.size(); ++j ) {
            vertexCoord = masterVolumeVertices[j].coord();
            std::copy( vertexCoord.begin(),
                       vertexCoord.end(),
                       &( masterVolumeVerticesCoordinates[3 * j] ) );
            masterVolumeVerticesGlobalIDs[j] = masterVolumeVertices[j].globalID();
        } // end for j
        std::vector<size_t> displacementIndices;
        getVectorIndicesFromGlobalIDs( masterVolumeVerticesGlobalIDs, displacementIndices );
        AMP_ASSERT( displacementIndices.size() == 24 );
        double displacementValues[24];
        if ( !skipDisplaceMesh ) {
            displacementFieldVector->getValuesByGlobalID(
                24, &( displacementIndices[0] ), &( displacementValues[0] ) );
        }
        else {
            std::fill( displacementValues, displacementValues + 24, 0.0 );
        } // end if
        double basis_functions_values_on_face[4];
        hex8_element_t::get_basis_functions_values_on_face(
            recvProjectionDataBuffer[i].d_SlaveVertexLocalCoordOnMasterFace,
            basis_functions_values_on_face );
        std::fill( sendProjOnMasterFaceData[i].d_Displacement,
                   sendProjOnMasterFaceData[i].d_Displacement + 3,
                   0.0 );
        unsigned int const *faceOrdering =
            hex8_element_t::get_face( recvProjectionDataBuffer[i].d_MasterFaceLocalIndex );
        for ( size_t j = 0; j < 4; ++j ) {
            for ( size_t k = 0; k < 3; ++k ) {
                sendProjOnMasterFaceData[i].d_Displacement[k] +=
                    displacementValues[3 * faceOrdering[j] + k] * basis_functions_values_on_face[j];
            } // end for k
        }     // end for j
              // nah
              //          std::transform(masterVolumeVerticesCoordinates,
              //          masterVolumeVerticesCoordinates+24,
              //          displacementValues, masterVolumeVerticesCoordinates, std::plus<double>());
        hex8_element_t masterVolumeElement( masterVolumeVerticesCoordinates );
        masterVolumeElement.compute_normal_to_face(
            recvProjectionDataBuffer[i].d_MasterFaceLocalIndex,
            recvProjectionDataBuffer[i].d_SlaveVertexLocalCoordOnMasterFace,
            sendProjOnMasterFaceData[i].d_NormalVector );
    } // end for i
    recvProjectionDataBuffer.clear();

    sendCnts.swap( recvCnts );
    sendDisps.swap( recvDisps );

    std::vector<AnotherDataWithNoName> recvProjOnMasterFaceData( nSendData );
    d_GlobalComm.allToAll(
        ( !( sendProjOnMasterFaceData.empty() ) ? &( sendProjOnMasterFaceData[0] ) : NULL ),
        &( sendCnts[0] ),
        &( sendDisps[0] ),
        ( !( recvProjOnMasterFaceData.empty() ) ? &( recvProjOnMasterFaceData[0] ) : NULL ),
        &( recvCnts[0] ),
        &( recvDisps[0] ),
        true );
    sendProjOnMasterFaceData.clear();

    std::transform( sendCnts.begin(),
                    sendCnts.end(),
                    sendCnts.begin(),
                    std::bind1st( std::multiplies<int>(), 3 ) );
    std::transform( sendDisps.begin(),
                    sendDisps.end(),
                    sendDisps.begin(),
                    std::bind1st( std::multiplies<int>(), 3 ) );
    std::transform( recvCnts.begin(),
                    recvCnts.end(),
                    recvCnts.begin(),
                    std::bind1st( std::multiplies<int>(), 3 ) );
    std::transform( recvDisps.begin(),
                    recvDisps.end(),
                    recvDisps.begin(),
                    std::bind1st( std::multiplies<int>(), 3 ) );

    std::vector<double> recvProjOnMasterFaceDisplacements( 3 * nSendData );
    d_GlobalComm.allToAll( ( !( sendProjOnMasterFaceDisplacements.empty() ) ?
                                 &( sendProjOnMasterFaceDisplacements[0] ) :
                                 NULL ),
                           &( sendCnts[0] ),
                           &( sendDisps[0] ),
                           ( !( recvProjOnMasterFaceDisplacements.empty() ) ?
                                 &( recvProjOnMasterFaceDisplacements[0] ) :
                                 NULL ),
                           &( recvCnts[0] ),
                           &( recvDisps[0] ),
                           true );
    sendProjOnMasterFaceDisplacements.clear();

    double activatedSlaveVertexShiftCorrection[3];
    double activatedSlaveVertexDisplacementValues[3];
    for ( size_t i = 0; i < nInactiveSlaveVerticesActivated; ++i ) {
        std::vector<size_t> activatedSlaveVertexDisplacementIndices;
        getVectorIndicesFromGlobalIDs(
            std::vector<AMP::Mesh::MeshElementID>( 1, d_ActiveSet[nActiveSlaveVerticesTmp + i] ),
            activatedSlaveVertexDisplacementIndices );
        AMP_ASSERT( activatedSlaveVertexDisplacementIndices.size() == 3 );
        if ( !skipDisplaceMesh ) {
            displacementFieldVector->getValuesByGlobalID(
                3,
                &( activatedSlaveVertexDisplacementIndices[0] ),
                activatedSlaveVertexDisplacementValues );
        }
        else {
            std::fill( activatedSlaveVertexDisplacementValues,
                       activatedSlaveVertexDisplacementValues + 3,
                       0.0 );
        } // end if
        make_vector_from_two_points( activatedSlaveVertexDisplacementValues,
                                     recvProjOnMasterFaceData[sendMap[i]].d_Displacement,
                                     activatedSlaveVertexShiftCorrection );
        std::transform( &( d_SlaveShift[3 * ( nActiveSlaveVerticesTmp + i )] ),
                        &( d_SlaveShift[3 * ( nActiveSlaveVerticesTmp + i )] ) + 3,
                        activatedSlaveVertexShiftCorrection,
                        &( d_SlaveShift[3 * ( nActiveSlaveVerticesTmp + i )] ),
                        std::minus<double>() );
        std::transform( recvProjOnMasterFaceData[sendMap[i]].d_NormalVector,
                        recvProjOnMasterFaceData[sendMap[i]].d_NormalVector + 3,
                        &( d_SlaveVerticesNormalVector[3 * ( nActiveSlaveVerticesTmp + i )] ),
                        std::bind1st( std::multiplies<double>(), -1.0 ) );
    } // end for i
    recvProjOnMasterFaceDisplacements.clear();
    recvProjOnMasterFaceData.clear();

    /** setup for apply */
    size_t const nConstraints = nActiveSlaveVerticesAfterUpdate;
    d_SendCnts.resize( npes );
    std::fill( d_SendCnts.begin(), d_SendCnts.end(), 0 );
    for ( size_t i = 0; i < 4 * nConstraints; ++i ) {
        //        ++d_SendCnts[d_MasterVerticesGlobalIDs[i].owner_rank()];
        ++d_SendCnts[invRankMap[d_MasterVerticesGlobalIDs[i].owner_rank()]];
    } // end for i
    d_SendDisps.resize( npes );
    d_SendDisps[0] = 0;
    for ( size_t i = 1; i < npes; ++i ) {
        d_SendDisps[i] = d_SendDisps[i - 1] + d_SendCnts[i - 1];
    } // end for i
    AMP_ASSERT( d_SendDisps[npes - 1] + d_SendCnts[npes - 1] ==
                4 * static_cast<int>( nConstraints ) );

    std::vector<int> tmpSendCnts( npes, 0 );
    d_MasterVerticesMap.resize( 4 * nConstraints, 4 * nConstraints );
    std::vector<AMP::Mesh::MeshElementID> sendMasterVerticesGlobalIDs(
        d_SendDisps[npes - 1] + d_SendCnts[npes - 1], AMP::Mesh::MeshElementID() );
    for ( size_t i = 0; i < 4 * nConstraints; ++i ) {
        //        size_t sendToRank = d_MasterVerticesGlobalIDs[i].owner_rank();
        size_t sendToRank      = invRankMap[d_MasterVerticesGlobalIDs[i].owner_rank()];
        d_MasterVerticesMap[i] = d_SendDisps[sendToRank] + tmpSendCnts[sendToRank];
        sendMasterVerticesGlobalIDs[d_MasterVerticesMap[i]] = d_MasterVerticesGlobalIDs[i];
        ++tmpSendCnts[sendToRank];
    } // end for i
    AMP_ASSERT( std::equal( tmpSendCnts.begin(), tmpSendCnts.end(), d_SendCnts.begin() ) );
    AMP_ASSERT( std::find( d_MasterVerticesMap.begin(),
                           d_MasterVerticesMap.end(),
                           4 * nConstraints ) == d_MasterVerticesMap.end() );
    tmpSendCnts.clear();

    d_RecvCnts.resize( npes );
    d_GlobalComm.allToAll( 1, &( d_SendCnts[0] ), &( d_RecvCnts[0] ) );
    d_RecvDisps.resize( npes );
    d_RecvDisps[0] = 0;
    for ( size_t i = 1; i < npes; ++i ) {
        d_RecvDisps[i] = d_RecvDisps[i - 1] + d_RecvCnts[i - 1];
    } // end for i
    d_RecvMasterVerticesGlobalIDs.resize( d_RecvDisps[npes - 1] + d_RecvCnts[npes - 1] );
    d_GlobalComm.allToAll(
        ( !( sendMasterVerticesGlobalIDs.empty() ) ? &( sendMasterVerticesGlobalIDs[0] ) : NULL ),
        &( d_SendCnts[0] ),
        &( d_SendDisps[0] ),
        ( !( d_RecvMasterVerticesGlobalIDs.empty() ) ? &( d_RecvMasterVerticesGlobalIDs[0] ) :
                                                       NULL ),
        &( d_RecvCnts[0] ),
        &( d_RecvDisps[0] ),
        true );

    d_TransposeSendCnts.resize( npes );
    d_TransposeSendDisps.resize( npes );
    d_TransposeRecvCnts.resize( npes );
    d_TransposeRecvDisps.resize( npes );
    std::copy( d_SendCnts.begin(), d_SendCnts.end(), d_TransposeSendCnts.begin() );
    std::copy( d_SendDisps.begin(), d_SendDisps.end(), d_TransposeSendDisps.begin() );
    std::copy( d_RecvCnts.begin(), d_RecvCnts.end(), d_TransposeRecvCnts.begin() );
    std::copy( d_RecvDisps.begin(), d_RecvDisps.end(), d_TransposeRecvDisps.begin() );

    std::swap_ranges( d_RecvCnts.begin(), d_RecvCnts.end(), d_SendCnts.begin() );
    std::swap_ranges( d_RecvDisps.begin(), d_RecvDisps.end(), d_SendDisps.begin() );

    getVectorIndicesFromGlobalIDs( d_ActiveSet, d_SlaveIndices );
    getVectorIndicesFromGlobalIDs( d_RecvMasterVerticesGlobalIDs, d_RecvMasterIndices );

    for ( size_t i = 0; i < npes; ++i ) {
        d_SendCnts[i] *= d_DOFsPerNode;
        d_SendDisps[i] *= d_DOFsPerNode;
        d_RecvCnts[i] *= d_DOFsPerNode;
        d_RecvDisps[i] *= d_DOFsPerNode;
        d_TransposeSendCnts[i] *= d_DOFsPerNode;
        d_TransposeSendDisps[i] *= d_DOFsPerNode;
        d_TransposeRecvCnts[i] *= d_DOFsPerNode;
        d_TransposeRecvDisps[i] *= d_DOFsPerNode;
    } // end for i

    // TODO: need to be converted at some point...
    d_SlaveVerticesGlobalIDs = d_ActiveSet;

    AMP_ASSERT( nActiveSlaveVerticesAfterUpdate == d_ActiveSet.size() );
    AMP_ASSERT( nInactiveSlaveVerticesAfterUpdate == d_InactiveSet.size() );
    size_t const nInactiveSlaveVerticesAfterUpdateReduced =
        d_GlobalComm.sumReduce( nInactiveSlaveVerticesAfterUpdate );
    size_t const nActiveSlaveVerticesAfterUpdateReduced =
        d_GlobalComm.sumReduce( nActiveSlaveVerticesAfterUpdate );
    size_t const nInactiveSlaveVerticesBeforeUpdateReduced =
        d_GlobalComm.sumReduce( nInactiveSlaveVerticesBeforeUpdate );
    size_t const nActiveSlaveVerticesBeforeUpdateReduced =
        d_GlobalComm.sumReduce( nActiveSlaveVerticesBeforeUpdate );
    size_t const nInactiveSlaveVerticesActivatedReduced =
        d_GlobalComm.sumReduce( nInactiveSlaveVerticesActivated );
    size_t const nActiveSlaveVerticesDeactivatedReduced =
        d_GlobalComm.sumReduce( nActiveSlaveVerticesDeactivated );
    if ( !rank ) {
        std::cout << "---------------------------------\n";
        std::cout << "BEFORE UPDATING ACTIVE SET || "
                  << "ACTIVE SET SIZE=" << nActiveSlaveVerticesBeforeUpdateReduced << "  "
                  << "INACTIVE SET SIZE=" << nInactiveSlaveVerticesBeforeUpdateReduced << " || \n";
        std::cout << "AFTER UPDATING ACTIVE SET  || "
                  << "ACTIVE SET SIZE=" << nActiveSlaveVerticesAfterUpdateReduced << "  "
                  << "INACTIVE SET SIZE=" << nInactiveSlaveVerticesAfterUpdateReduced << " || \n";
        std::cout << "UPDATE ACTIVE SET -> "
                  << "# INACTIVE VERTICES ACTIVATED=" << nInactiveSlaveVerticesActivatedReduced
                  << "  "
                  << "# ACTIVE VERTICES DESCTIVATED=" << nActiveSlaveVerticesDeactivatedReduced
                  << "\n";
        std::cout << "---------------------------------\n";
    } // end if

    /** move the mesh back */
    //      if (!skipDisplaceMesh) {
    //        displacementFieldVector->scale(-1.0);
    //        d_Mesh->displaceMesh(displacementFieldVector);
    //        displacementFieldVector->scale(-1.0);
    //      } // end if

    return d_GlobalComm.sumReduce( nInactiveSlaveVerticesActivated +
                                   nActiveSlaveVerticesDeactivated );
}

void NodeToFaceContactOperator::reset( const AMP::shared_ptr<OperatorParameters> &params )
{
    AMP_INSIST( ( params != NULL ), "NULL parameter" );
    AMP_INSIST( ( ( params->d_db ) != NULL ), "NULL database" );
}

void NodeToFaceContactOperator::getVectorIndicesFromGlobalIDs(
    const std::vector<AMP::Mesh::MeshElementID> &globalIDs, std::vector<size_t> &vectorIndices )
{
    std::vector<size_t> tmpIndices;
    std::vector<AMP::Mesh::MeshElementID>::const_iterator globalIDsConstIterator =
                                                              globalIDs.begin(),
                                                          globalIDsConstIterator_end =
                                                              globalIDs.end();
    vectorIndices.resize( globalIDs.size() * d_DOFsPerNode );
    std::vector<size_t>::iterator vectorIndicesIterator = vectorIndices.begin();
    for ( ; globalIDsConstIterator != globalIDsConstIterator_end; ++globalIDsConstIterator ) {
        d_DOFManager->getDOFs( *globalIDsConstIterator, tmpIndices );
        AMP_ASSERT( *globalIDsConstIterator != AMP::Mesh::MeshElementID() );
        AMP_ASSERT( tmpIndices.size() == d_DOFsPerNode );
        std::copy( tmpIndices.begin(), tmpIndices.end(), vectorIndicesIterator );
        for ( size_t i = 0; i < d_DOFsPerNode; ++i ) {
            ++vectorIndicesIterator;
        }
    } // end for
    AMP_ASSERT( vectorIndicesIterator == vectorIndices.end() );
}

void NodeToFaceContactOperator::copyMasterToSlave( AMP::LinearAlgebra::Vector::shared_ptr u )
{
    /** send and receive the master values */
    AMP::AMP_MPI comm = d_GlobalComm;
    //  AMP::AMP_MPI comm = u->getComm();
    size_t npes = comm.getSize();
    //  size_t rank = comm.getRank();

    std::vector<double> sendMasterValues( d_SendDisps[npes - 1] + d_SendCnts[npes - 1] );
    for ( size_t i = 0; i < npes; ++i ) {
        for ( int j = 0; j < d_SendCnts[i]; j += d_DOFsPerNode ) {
            size_t k = d_SendDisps[i] + j;
            u->getLocalValuesByGlobalID(
                d_DOFsPerNode, &( d_RecvMasterIndices[k] ), &( sendMasterValues[k] ) );
        } // end for j
    }     // end for i

    std::vector<double> recvMasterValues( d_RecvDisps[npes - 1] + d_RecvCnts[npes - 1] );
    comm.allToAll( ( !( sendMasterValues.empty() ) ? &( sendMasterValues[0] ) : NULL ),
                   &( d_SendCnts[0] ),
                   &( d_SendDisps[0] ),
                   ( !( recvMasterValues.empty() ) ? &( recvMasterValues[0] ) : NULL ),
                   &( d_RecvCnts[0] ),
                   &( d_RecvDisps[0] ),
                   true );
    sendMasterValues.clear();

    /** compute slave values */
    std::vector<double> slaveValues( d_SlaveIndices.size(), 0.0 );
    if ( ( d_ContactIsFrictionless ) && ( !d_SlaveVerticesGlobalIDs.empty() ) ) {
        u->getLocalValuesByGlobalID(
            d_SlaveIndices.size(), &( d_SlaveIndices[0] ), &( slaveValues[0] ) );
    } // end if
    double tmp[3];

    for ( size_t i = 0; i < d_SlaveVerticesGlobalIDs.size(); ++i ) {
        std::fill( tmp, tmp + 3, 0.0 );
        for ( size_t j = 0; j < 4; ++j ) {
            for ( size_t k = 0; k < d_DOFsPerNode; ++k ) {
                //            slaveValues[d_DOFsPerNode*i+k] += d_MasterShapeFunctionsValues[4*i+j]
                //            *
                //            recvMasterValues[d_DOFsPerNode*d_MasterVerticesMap[4*i+j]+k];
                tmp[k] += d_MasterShapeFunctionsValues[4 * i + j] *
                          recvMasterValues[d_DOFsPerNode * d_MasterVerticesMap[4 * i + j] + k];
            } // end for k
        }     // end for j
        if ( d_ContactIsFrictionless ) {
            orthogonalize_vector_against_direction( &( d_SlaveVerticesNormalVector[3 * i] ),
                                                    &( slaveValues[d_DOFsPerNode * i] ) );
            project_vector_onto_direction( &( d_SlaveVerticesNormalVector[3 * i] ), tmp );
        } // end if
        std::transform( tmp,
                        tmp + 3,
                        &( slaveValues[d_DOFsPerNode * i] ),
                        &( slaveValues[d_DOFsPerNode * i] ),
                        std::plus<double>() );
    } // end for i

    if ( !d_SlaveVerticesGlobalIDs.empty() ) {
        u->setLocalValuesByGlobalID(
            d_SlaveIndices.size(), &( d_SlaveIndices[0] ), &( slaveValues[0] ) );
    } // end if
}

void NodeToFaceContactOperator::addSlaveToMaster( AMP::LinearAlgebra::Vector::shared_ptr u )
{
    /** send and receive slave value times shape functions values */
    AMP::AMP_MPI comm = d_GlobalComm;
    //  AMP::AMP_MPI comm = r->getComm();
    size_t npes = comm.getSize();

    std::vector<double> sendAddToMasterValues( d_TransposeSendDisps[npes - 1] +
                                               d_TransposeSendCnts[npes - 1] );
    for ( size_t i = 0; i < d_SlaveVerticesGlobalIDs.size(); ++i ) {
        for ( size_t j = 0; j < 4; ++j ) {
            u->getLocalValuesByGlobalID(
                d_DOFsPerNode,
                &( d_SlaveIndices[d_DOFsPerNode * i] ),
                &( sendAddToMasterValues[d_DOFsPerNode * d_MasterVerticesMap[4 * i + j]] ) );
            if ( d_ContactIsFrictionless ) {
                project_vector_onto_direction(
                    &( d_SlaveVerticesNormalVector[3 * i] ),
                    &( sendAddToMasterValues[d_DOFsPerNode * d_MasterVerticesMap[4 * i + j]] ) );
            } // end if
            for ( size_t k = 0; k < d_DOFsPerNode; ++k ) {
                sendAddToMasterValues[d_DOFsPerNode * d_MasterVerticesMap[4 * i + j] + k] *=
                    d_MasterShapeFunctionsValues[4 * i + j];
            } // end for k
        }     // end for j
    }         // end for i

    std::vector<double> recvAddToMasterValues( d_TransposeRecvDisps[npes - 1] +
                                               d_TransposeRecvCnts[npes - 1] );
    comm.allToAll( ( !( sendAddToMasterValues.empty() ) ? &( sendAddToMasterValues[0] ) : NULL ),
                   &( d_TransposeSendCnts[0] ),
                   &( d_TransposeSendDisps[0] ),
                   ( !( recvAddToMasterValues.empty() ) ? &( recvAddToMasterValues[0] ) : NULL ),
                   &( d_TransposeRecvCnts[0] ),
                   &( d_TransposeRecvDisps[0] ),
                   true );
    sendAddToMasterValues.clear();

    /** add slave value times shape functions values to master values and set slave values to zero
     */
    if ( !d_RecvMasterVerticesGlobalIDs.empty() ) {
        u->addLocalValuesByGlobalID( d_RecvMasterIndices.size(),
                                     &( d_RecvMasterIndices[0] ),
                                     &( recvAddToMasterValues[0] ) );
    } // end if
}

void NodeToFaceContactOperator::setSlaveToZero( AMP::LinearAlgebra::Vector::shared_ptr u )
{
    if ( ( d_ContactIsFrictionless ) && ( !d_SlaveVerticesGlobalIDs.empty() ) ) {
        AMP_ASSERT( d_DOFsPerNode == 3 );
        std::vector<double> slaveValues( d_SlaveIndices.size() );
        u->getLocalValuesByGlobalID(
            d_SlaveIndices.size(), &( d_SlaveIndices[0] ), &( slaveValues[0] ) );
        for ( size_t i = 0; i < d_SlaveVerticesGlobalIDs.size(); ++i ) {
            orthogonalize_vector_against_direction( &( d_SlaveVerticesNormalVector[3 * i] ),
                                                    &( slaveValues[d_DOFsPerNode * i] ) );
        } // end for i
        u->setLocalValuesByGlobalID(
            d_SlaveIndices.size(), &( d_SlaveIndices[0] ), &( slaveValues[0] ) );
    }
    else {
        this->ConstraintsEliminationOperator::setSlaveToZero( u );
    } // end if
}

void NodeToFaceContactOperator::addShiftToSlave( AMP::LinearAlgebra::Vector::shared_ptr u )
{
    if ( ( d_ContactIsFrictionless ) && ( !d_SlaveVerticesGlobalIDs.empty() ) ) {
        AMP_ASSERT( d_DOFsPerNode == 3 );
        std::vector<double> addToSlaveValues( d_SlaveIndices.size() );
        for ( size_t i = 0; i < d_SlaveVerticesGlobalIDs.size(); ++i ) {
            std::copy( &( d_SlaveShift[3 * i] ),
                       &( d_SlaveShift[3 * i] ) + 3,
                       &( addToSlaveValues[d_DOFsPerNode * i] ) );
            project_vector_onto_direction( &( d_SlaveVerticesNormalVector[3 * i] ),
                                           &( addToSlaveValues[d_DOFsPerNode * i] ) );
        } // end for i
        u->addLocalValuesByGlobalID(
            d_SlaveIndices.size(), &( d_SlaveIndices[0] ), &( addToSlaveValues[0] ) );
    }
    else {
        this->ConstraintsEliminationOperator::addShiftToSlave( u );
    } // end if
}

} // end namespace Operator
} // end namespace AMP
