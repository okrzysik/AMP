#include "AMP/utils/Utilities.h"
#include "AMP/vectors/VectorIndexer.h"
#include "AMP/vectors/data/VectorData.h"

#include <iostream>
#include <memory>
#include <vector>


namespace AMP::LinearAlgebra {


/************************************************************************
 * Constructors                                                          *
 ************************************************************************/
CommunicationListParameters::CommunicationListParameters() : d_comm( AMP_COMM_NULL )
{
    d_localsize = (size_t) -1;
}
CommunicationListParameters::CommunicationListParameters( const CommunicationListParameters &rhs )
    : d_comm( rhs.d_comm ), d_localsize( rhs.d_localsize ), d_remote_DOFs( rhs.d_remote_DOFs )
{
}
CommunicationList::CommunicationList() : d_iBegin( 0 ), d_iNumRows( 0 ), d_iTotalRows( 0 ) {}
CommunicationList::CommunicationList( std::shared_ptr<const CommunicationListParameters> params )
    : d_comm( params->d_comm ), d_iNumRows( params->d_localsize )
{
    // Check the input parameters
    AMP_ASSERT( ( params->d_localsize >> 48 ) == 0 );
    AMP_ASSERT( params->d_comm != AMP_MPI( AMP_COMM_NULL ) );
    // Get the partition (the total number of DOFs for all ranks <= current rank)
    std::vector<size_t> partition( d_comm.getSize(), 0 );
    d_comm.allGather<size_t>( params->d_localsize, &partition[0] );
    for ( int i = 1; i < d_comm.getSize(); i++ )
        partition[i] += partition[i - 1];
    // Get the first DOF on the current rank
    d_iBegin = partition[d_comm.getRank()] - params->d_localsize;
    // Get the total number of DOFs
    d_iTotalRows = partition[d_comm.getSize() - 1];
    // Construct the communication arrays
    buildCommunicationArrays( params->d_remote_DOFs, partition, d_comm.getRank() );
}
CommunicationList::CommunicationList( size_t local, const AMP_MPI &comm ) : d_comm( comm )
{
    size_t size    = comm.getSize();
    size_t lastRow = comm.sumScan( local );
    d_ReceiveSizes.resize( size );
    d_ReceiveDisplacements.resize( size );
    d_ReceiveDOFList.resize( 0 );
    d_SendSizes.resize( size );
    d_SendDisplacements.resize( size );
    d_SendDOFList.resize( 0 );
    d_iBegin     = lastRow - local;
    d_comm       = comm;
    d_iNumRows   = local;
    d_iTotalRows = comm.bcast( lastRow, size - 1 );
}


/************************************************************************
 * Subset                                                                *
 ************************************************************************/
std::shared_ptr<CommunicationList> CommunicationList::subset( std::shared_ptr<VectorIndexer> ndx )
{
    auto retVal = std::make_shared<CommunicationList>( 0, d_comm );

    retVal->d_ReceiveSizes.resize( std::max( d_SendSizes.size(), (size_t) 1 ) );
    retVal->d_ReceiveDisplacements.resize( std::max( d_SendSizes.size(), (size_t) 1 ) );
    retVal->d_SendSizes.resize( std::max( d_SendSizes.size(), (size_t) 1 ) );
    retVal->d_SendDisplacements.resize( std::max( d_SendSizes.size(), (size_t) 1 ) );
    size_t curOff = 0;
    for ( size_t i = 0; i != d_SendSizes.size(); i++ ) {
        retVal->d_SendDisplacements[i] = curOff;
        retVal->d_SendSizes[i]         = 0;
        for ( int j = 0; j != d_SendSizes[i]; j++ ) {
            if ( ndx->isInSub( d_SendDOFList[d_SendDisplacements[i] + j] ) ) {
                retVal->d_SendSizes[i]++;
                retVal->d_SendDOFList.push_back(
                    ndx->getSubID( d_SendDOFList[d_SendDisplacements[i] + j] ) );
                curOff++;
            }
        }
    }
    d_comm.allToAll( 1, retVal->d_SendSizes.data(), retVal->d_ReceiveSizes.data() );

    retVal->d_ReceiveDisplacements[0] = 0;
    size_t ii;
    for ( ii = 1; ii < d_SendSizes.size(); ii++ ) {
        retVal->d_ReceiveDisplacements[ii] =
            retVal->d_ReceiveSizes[ii - 1] + retVal->d_ReceiveDisplacements[ii - 1];
    }

    retVal->d_ReceiveDOFList.resize( retVal->d_ReceiveDisplacements[ii - 1] +
                                     retVal->d_ReceiveSizes[ii - 1] );

    d_comm.allToAll( retVal->d_SendDOFList.data(),
                     retVal->d_SendSizes.data(),
                     retVal->d_SendDisplacements.data(),
                     retVal->d_ReceiveDOFList.data(),
                     retVal->d_ReceiveSizes.data(),
                     retVal->d_ReceiveDisplacements.data(),
                     true );


    retVal->d_iNumRows = 0;
    for ( size_t i = getStartGID(); i != getStartGID() + numLocalRows(); i++ ) {
        if ( ndx->isInSub( i ) )
            retVal->d_iNumRows++;
    }
    d_comm.sumScan( &( retVal->d_iNumRows ), &( retVal->d_iTotalRows ), 1 );
    retVal->d_iBegin     = retVal->d_iTotalRows - retVal->d_iNumRows;
    retVal->d_iTotalRows = d_comm.bcast( retVal->d_iTotalRows, retVal->d_SendSizes.size() - 1 );

    return retVal;
}


/************************************************************************
 * Get ids                                                               *
 ************************************************************************/
size_t CommunicationList::getLocalGhostID( size_t GID ) const
{
    // Search d_ReceiveDOFList for GID
    // Note: d_ReceiveDOFList must be sorted for this to work
    AMP_INSIST( !d_ReceiveDOFList.empty(),
                "Tried to access ghost entry, but vector does not contain ghosts" );
    size_t pos = AMP::Utilities::findfirst( d_ReceiveDOFList, (size_t) GID );
    bool found = pos < d_ReceiveDOFList.size();
    if ( found ) {
        if ( d_ReceiveDOFList[pos] != GID )
            found = false;
    }
    if ( !found ) {
        std::cout << "GID = " << GID << std::endl;
        AMP_ERROR( "GID was not found in the ghost list" );
    }
    return pos;
}


/************************************************************************
 * Build the arrays                                                      *
 ************************************************************************/
void CommunicationList::buildCommunicationArrays( const std::vector<size_t> &DOFs,
                                                  const std::vector<size_t> &partitionInfo,
                                                  int commRank )
{
    d_iBegin = commRank ? partitionInfo[commRank - 1] : 0;

    size_t commSize = partitionInfo.size();

    // Check if we are working in serial
    if ( commSize == 1 ) {
        AMP_INSIST( DOFs.empty(),
                    "Error in communication list, remote DOFs are present for a serial vector" );
        d_ReceiveSizes.resize( 1, 0 );
        d_ReceiveDisplacements.resize( 1, 0 );
        d_SendSizes.resize( 1, 0 );
        d_SendDisplacements.resize( 1, 0 );
        d_ReceiveDOFList.resize( 0 );
        d_SendDOFList.resize( 0 );
        return;
    }

    // Copy (and sort) the DOFs in d_ReceiveDOFList
    d_ReceiveDOFList = DOFs;
    AMP::Utilities::quicksort( d_ReceiveDOFList );

    // Determine the number of DOFs received from each processor (requires DOFs to be sorted)
    d_ReceiveSizes.resize( commSize, 0 );
    d_ReceiveDisplacements.resize( commSize, 0 );
    size_t rank  = 0;
    size_t start = 0;
    size_t index = 0;
    while ( index < d_ReceiveDOFList.size() ) {
        if ( d_ReceiveDOFList[index] < partitionInfo[rank] ) {
            // Move to the next DOF
            index++;
        } else {
            // Store the number of DOFs with the given rank, and move to the next rank
            d_ReceiveDisplacements[rank] = start;
            d_ReceiveSizes[rank]         = index - start;
            start                        = index;
            rank++;
        }
    }
    d_ReceiveDisplacements[rank] = start;
    d_ReceiveSizes[rank]         = index - start;
    AMP_ASSERT( d_ReceiveSizes[commRank] == 0 );

    // Get the send sizes and displacements
    d_SendSizes.resize( commSize, 0 );
    d_SendDisplacements.resize( commSize, 0 );
    d_comm.allToAll( 1, &( d_ReceiveSizes[0] ), &( d_SendSizes[0] ) );
    d_SendDisplacements[0] = 0;
    for ( size_t i = 1; i < commSize; i++ )
        d_SendDisplacements[i] = d_SendDisplacements[i - 1] + d_SendSizes[i - 1];
    size_t send_buf_size = d_SendDisplacements[commSize - 1] + d_SendSizes[commSize - 1];

    // Get the send DOFs (requires the recv DOFs to be sorted)
    d_SendDOFList.resize( send_buf_size );
    d_comm.allToAll( d_ReceiveDOFList.data(),
                     d_ReceiveSizes.data(),
                     d_ReceiveDisplacements.data(),
                     d_SendDOFList.data(),
                     d_SendSizes.data(),
                     d_SendDisplacements.data(),
                     true );
}


/************************************************************************
 * set/recv data                                                         *
 ************************************************************************/
void CommunicationList::scatter_set( VectorData &vec ) const
{
    if ( d_SendSizes.empty() )
        return;
    // Pack the set buffers
    std::vector<double> send( getVectorSendBufferSize() );
    std::vector<double> recv( getVectorReceiveBufferSize() );
    vec.getLocalValuesByGlobalID( send.size(), d_SendDOFList.data(), send.data() );
    // Communicate
    d_comm.allToAll<double>( send.data(),
                             d_SendSizes.data(),
                             d_SendDisplacements.data(),
                             recv.data(),
                             (int *) d_ReceiveSizes.data(),
                             (int *) d_ReceiveDisplacements.data(),
                             true );
    // Unpack the set buffers
    vec.setGhostValuesByGlobalID( recv.size(), d_ReceiveDOFList.data(), recv.data() );
}
void CommunicationList::scatter_add( VectorData &vec ) const
{
    if ( d_SendSizes.empty() )
        return;
    // Pack the add buffers
    std::vector<double> send( getVectorReceiveBufferSize() );
    std::vector<double> recv( getVectorSendBufferSize() );
    vec.getGhostAddValuesByGlobalID( send.size(), d_ReceiveDOFList.data(), send.data() );
    // Communicate
    d_comm.allToAll<double>( send.data(),
                             d_ReceiveSizes.data(),
                             d_ReceiveDisplacements.data(),
                             recv.data(),
                             (int *) d_SendSizes.data(),
                             (int *) d_SendDisplacements.data(),
                             true );
    // Unpack the add buffers
    vec.addLocalValuesByGlobalID( recv.size(), d_SendDOFList.data(), recv.data() );
}


/************************************************************************
 * Misc. functions                                                       *
 ************************************************************************/
size_t CommunicationList::getStartGID() const { return d_iBegin; }
size_t CommunicationList::numLocalRows() const { return d_iNumRows; }
size_t CommunicationList::getTotalSize() const { return d_iTotalRows; }
size_t CommunicationList::getVectorSendBufferSize() const { return d_SendDOFList.size(); }
size_t CommunicationList::getVectorReceiveBufferSize() const { return d_ReceiveDOFList.size(); }
const std::vector<size_t> &CommunicationList::getGhostIDList() const { return d_ReceiveDOFList; }
const std::vector<size_t> &CommunicationList::getReplicatedIDList() const { return d_SendDOFList; }
const AMP_MPI &CommunicationList::getComm() const { return d_comm; }


} // namespace AMP::LinearAlgebra
