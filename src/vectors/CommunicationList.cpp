#include "AMP/utils/Utilities.h"
#include "AMP/vectors/VectorIndexer.h"
#include "AMP/vectors/data/VectorData.h"
#include <memory>

#include <iostream>
#include <vector>


namespace AMP::LinearAlgebra {


/************************************************************************
 * Some simple fuctions to get a pointer to the data in a std::vector,   *
 * where the vector may be empty.                                        *
 ************************************************************************/
template<typename T>
static T *getPtr( std::vector<T> &in )
{
    T *retVal = nullptr;
    if ( !in.empty() )
        retVal = &( in[0] );
    return retVal;
}
template<typename T>
static T *getPtr( const std::vector<T> &in )
{
    T *retVal = nullptr;
    if ( !in.empty() )
        retVal = (T *) &( in[0] );
    return retVal;
}


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
CommunicationList::CommunicationList()
    : d_iBegin( 0 ), d_iNumRows( 0 ), d_iTotalRows( 0 ), d_bFinalized( false )
{
}
CommunicationList::CommunicationList( std::shared_ptr<CommunicationListParameters> params )
    : d_comm( params->d_comm ), d_iNumRows( params->d_localsize ), d_bFinalized( false )
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
std::shared_ptr<CommunicationList> CommunicationList::createEmpty( size_t local, AMP_MPI comm )
{
    int size    = comm.getSize();
    auto retVal = new CommunicationList;
    retVal->d_ReceiveSizes.resize( size );
    retVal->d_ReceiveDisplacements.resize( size );
    retVal->d_ReceiveDOFList.resize( 0 );

    retVal->d_SendSizes.resize( size );
    retVal->d_SendDisplacements.resize( size );
    retVal->d_SendDOFList.resize( 0 );

    size_t lastRow = 0;
    comm.sumScan( &local, &( lastRow ), 1 );
    retVal->d_iBegin     = lastRow - local;
    retVal->d_comm       = comm;
    retVal->d_iNumRows   = local;
    retVal->d_iTotalRows = comm.bcast( lastRow, size - 1 );
    retVal->d_bFinalized = true;

    return std::shared_ptr<CommunicationList>( retVal );
}


/************************************************************************
 * All other functions                                                   *
 ************************************************************************/
std::shared_ptr<CommunicationList> CommunicationList::subset( std::shared_ptr<VectorIndexer> ndx )
{
    auto retVal = new CommunicationList;

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
    d_comm.allToAll( 1, getPtr( retVal->d_SendSizes ), getPtr( retVal->d_ReceiveSizes ) );

    retVal->d_ReceiveDisplacements[0] = 0;
    size_t ii;
    for ( ii = 1; ii < d_SendSizes.size(); ii++ ) {
        retVal->d_ReceiveDisplacements[ii] =
            retVal->d_ReceiveSizes[ii - 1] + retVal->d_ReceiveDisplacements[ii - 1];
    }

    retVal->d_ReceiveDOFList.resize( retVal->d_ReceiveDisplacements[ii - 1] +
                                     retVal->d_ReceiveSizes[ii - 1] );

    d_comm.allToAll( getPtr( retVal->d_SendDOFList ),
                     &( retVal->d_SendSizes[0] ),
                     &( retVal->d_SendDisplacements[0] ),
                     getPtr( retVal->d_ReceiveDOFList ),
                     &( retVal->d_ReceiveSizes[0] ),
                     &( retVal->d_ReceiveDisplacements[0] ),
                     true );

    retVal->d_comm = d_comm;

    retVal->d_iNumRows = 0;
    for ( size_t i = getStartGID(); i != getStartGID() + numLocalRows(); i++ ) {
        if ( ndx->isInSub( i ) )
            retVal->d_iNumRows++;
    }
    d_comm.sumScan( &( retVal->d_iNumRows ), &( retVal->d_iTotalRows ), 1 );
    retVal->d_iBegin     = retVal->d_iTotalRows - retVal->d_iNumRows;
    retVal->d_iTotalRows = d_comm.bcast( retVal->d_iTotalRows, retVal->d_SendSizes.size() - 1 );

    retVal->d_bFinalized = true;

    return std::shared_ptr<CommunicationList>( retVal );
}

void CommunicationList::packReceiveBuffer( std::vector<double> &recv, const VectorData &vec ) const
{
    AMP_ASSERT( recv.size() == d_ReceiveDOFList.size() );
    if ( recv.empty() )
        return;
    vec.getGhostAddValuesByGlobalID(
        (int) recv.size(), getPtr( d_ReceiveDOFList ), getPtr( recv ) );
}

void CommunicationList::packSendBuffer( std::vector<double> &send, const VectorData &vec ) const
{
    AMP_ASSERT( send.size() == d_SendDOFList.size() );
    if ( send.empty() )
        return;
    vec.getLocalValuesByGlobalID( (int) send.size(), getPtr( d_SendDOFList ), getPtr( send ) );
}

void CommunicationList::unpackReceiveBufferSet( const std::vector<double> &recv,
                                                VectorData &vec ) const
{
    AMP_ASSERT( recv.size() == d_ReceiveDOFList.size() );
    vec.setGhostValuesByGlobalID( (int) recv.size(), getPtr( d_ReceiveDOFList ), getPtr( recv ) );
}

void CommunicationList::unpackSendBufferAdd( const std::vector<double> &recv,
                                             VectorData &vec ) const
{
    AMP_ASSERT( recv.size() == d_SendDOFList.size() );
    vec.addLocalValuesByGlobalID( (int) recv.size(), getPtr( d_SendDOFList ), getPtr( recv ) );
}

void CommunicationList::unpackSendBufferSet( const std::vector<double> &recv,
                                             VectorData &vec ) const
{
    AMP_ASSERT( recv.size() == d_SendDOFList.size() );
    vec.setLocalValuesByGlobalID( (int) recv.size(), getPtr( d_SendDOFList ), getPtr( recv ) );
}

size_t CommunicationList::getLocalGhostID( size_t GID ) const
{
    // Search d_ReceiveDOFList for GID
    // Note: d_ReceiveDOFList must be sorted for this to work
    AMP_INSIST( !d_ReceiveDOFList.empty(),
                "Tried to access ghost entry, but vector does not contain ghosts" );
    size_t pos = AMP::Utilities::findfirst( d_ReceiveDOFList, (size_t) GID );
    bool found = pos < d_ReceiveDOFList.size();
    if ( found ) {
        if ( d_ReceiveDOFList[pos] != GID ) {
            found = false;
        }
    }
    if ( !found ) {
        std::cout << "GID = " << GID << std::endl;
        AMP_ERROR( "GID was not found in the ghost list" );
    }
    return pos;
}


void CommunicationList::buildCommunicationArrays( std::vector<size_t> &DOFs,
                                                  std::vector<size_t> &partitionInfo,
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

    // Determine the number of DOFs received from each processor (this requires the DOFs to be
    // sorted)
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
    d_comm.allToAll( getPtr( d_ReceiveDOFList ),
                     &( d_ReceiveSizes[0] ),
                     &( d_ReceiveDisplacements[0] ),
                     getPtr( d_SendDOFList ),
                     &( d_SendSizes[0] ),
                     &( d_SendDisplacements[0] ),
                     true );
}


void CommunicationList::scatter_set( std::vector<double> &in, std::vector<double> &out ) const
{
    if ( d_SendSizes.empty() )
        return;
    double *send_buf = getPtr( in );
    auto *send_sizes = (int *) &( d_SendSizes[0] );
    auto *send_disps = (int *) &( d_SendDisplacements[0] );
    double *recv_buf = getPtr( out );
    auto *recv_sizes = (int *) &( d_ReceiveSizes[0] );
    auto *recv_disps = (int *) &( d_ReceiveDisplacements[0] );

    d_comm.allToAll( send_buf, send_sizes, send_disps, recv_buf, recv_sizes, recv_disps, true );
}

void CommunicationList::scatter_add( std::vector<double> &in, std::vector<double> &out ) const
{
    double *send_buf = getPtr( in );
    auto *send_sizes = (int *) &( d_ReceiveSizes[0] );
    auto *send_disps = (int *) &( d_ReceiveDisplacements[0] );
    double *recv_buf = getPtr( out );
    auto *recv_sizes = (int *) &( d_SendSizes[0] );
    auto *recv_disps = (int *) &( d_SendDisplacements[0] );

    d_comm.allToAll( send_buf, send_sizes, send_disps, recv_buf, recv_sizes, recv_disps, true );
}

size_t CommunicationList::numLocalRows() const { return d_iNumRows; }

void CommunicationList::finalize()
{
    if ( !d_bFinalized ) {
        d_bFinalized = true;
    }
}

size_t CommunicationList::getTotalSize() const { return d_iTotalRows; }


CommunicationList::~CommunicationList() = default;

const std::vector<size_t> &CommunicationList::getGhostIDList() const { return d_ReceiveDOFList; }

const std::vector<size_t> &CommunicationList::getReplicatedIDList() const { return d_SendDOFList; }

size_t CommunicationList::getVectorReceiveBufferSize() const { return d_ReceiveDOFList.size(); }

size_t CommunicationList::getVectorSendBufferSize() const { return d_SendDOFList.size(); }

size_t CommunicationList::getStartGID() const { return d_iBegin; }

const AMP_MPI &CommunicationList::getComm() const { return d_comm; }

} // namespace AMP::LinearAlgebra
