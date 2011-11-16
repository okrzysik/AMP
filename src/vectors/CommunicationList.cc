#include <vector>
#include "Vector.h"
#include "utils/Utilities.h"
#include "boost/shared_ptr.hpp"

//class Vector;
//typedef  boost::shared_ptr<Vector>   Vector_shared_ptr;


namespace AMP {
namespace LinearAlgebra {


/************************************************************************
* Some simple fuctions to get a pointer to the data in a std::vector,   *
* where the vector may be empty.                                        *
************************************************************************/
template <typename T>
static T* getBufferToAvoidDebugVectorCrashing( std::vector<T> &in )
{
    T *retVal = 0;
    if ( in.size() > 0 ) retVal = &(in[0]);
    return retVal;
}
template <typename T>
static const T* getBufferToAvoidDebugVectorCrashing( const std::vector<T> &in )
{
    const T *retVal = 0;
    if ( in.size() > 0 ) retVal = &(in[0]);
    return retVal;
}


/************************************************************************
* Constructors                                                          *
************************************************************************/
CommunicationList::CommunicationList ( CommunicationListParameters::shared_ptr params ):
    d_comm ( params->d_comm ),
    d_iNumRows ( params->d_localsize ),
    d_bFinalized ( false )
{
    // Get the partition (the total number of DOFs for all ranks <= current rank)
    std::vector<unsigned int> partition(d_comm.getSize(),0);
    d_comm.allGather<unsigned int>( (unsigned int) params->d_localsize, &partition[0] );
    for (int i=1; i<d_comm.getSize(); i++)
        partition[i] += partition[i-1];
    d_comm.sumScan((int*)&d_iNumRows,(int*)&d_iTotalRows,1);
    // Get the first DOF on the current rank
    d_iBegin = partition[d_comm.getRank()] - params->d_localsize;
    // Get the total number of DOFs
    d_iTotalRows = partition[d_comm.getSize()-1];
    // Construct the communication arrays
    buildCommunicationArrays( params->d_remote_DOFs, partition, d_comm.getRank() );
}


/************************************************************************
* All other functions                                                   *
************************************************************************/
  CommunicationList::shared_ptr  CommunicationList::subset ( VectorIndexer::shared_ptr ndx )
  {
    CommunicationList *retVal = new CommunicationList;

    retVal->d_ReceiveSizes.resize ( std::max ( d_SendSizes.size() , (size_t)1 ) );
    retVal->d_ReceiveDisplacements.resize ( std::max ( d_SendSizes.size() ,(size_t) 1 ) );
    retVal->d_SendSizes.resize ( std::max ( d_SendSizes.size() ,(size_t) 1 ) );
    retVal->d_SendDisplacements.resize ( std::max ( d_SendSizes.size() ,(size_t) 1 ) );
    size_t curOff = 0;
    for ( size_t i = 0 ; i != d_SendSizes.size() ; i++ )
    {
      retVal->d_SendDisplacements[i] = curOff;
      retVal->d_SendSizes[i] = 0;
      for ( size_t j = 0 ; j != d_SendSizes[i] ; j++ )
      {
        if ( ndx->isInSub ( d_SendDOFList[d_SendDisplacements[i] + j] ) )
        {
          retVal->d_SendSizes[i]++;
          retVal->d_SendDOFList.push_back ( ndx->getSubID ( d_SendDOFList[d_SendDisplacements[i] + j] ) );
          curOff++;
        }
      }
    }
    d_comm.allToAll(1, getBufferToAvoidDebugVectorCrashing(retVal->d_SendSizes), 
        getBufferToAvoidDebugVectorCrashing(retVal->d_ReceiveSizes) );

    retVal->d_ReceiveDisplacements[0] = 0;
    size_t ii;
    for ( ii = 1 ; ii < d_SendSizes.size() ; ii++ )
    {
      retVal->d_ReceiveDisplacements[ii] = retVal->d_ReceiveSizes[ii-1] + retVal->d_ReceiveDisplacements[ii-1];
    }

    retVal->d_ReceiveDOFList.resize ( retVal->d_ReceiveDisplacements[ii-1] + retVal->d_ReceiveSizes[ii-1] );

    d_comm.allToAll( getBufferToAvoidDebugVectorCrashing ( retVal->d_SendDOFList ) ,
                    (int *)&(retVal->d_SendSizes[0]) ,
                    (int *)&(retVal->d_SendDisplacements[0] ) ,
                    getBufferToAvoidDebugVectorCrashing ( retVal->d_ReceiveDOFList ) ,
                    (int *)&(retVal->d_ReceiveSizes[0]) ,
                    (int *)&(retVal->d_ReceiveDisplacements[0] ) ,
                    true );

    retVal->d_comm = d_comm;

    retVal->d_iNumRows = 0;
    for ( size_t i = getStartGID() ; i != getStartGID() + numLocalRows() ; i++ )
    {
      if ( ndx->isInSub ( i ) )
        retVal->d_iNumRows ++;
    }
    d_comm.sumScan( (int*)&(retVal->d_iNumRows), (int*)&(retVal->d_iTotalRows), 1);
    retVal->d_iBegin = retVal->d_iTotalRows - retVal->d_iNumRows;
    retVal->d_iTotalRows = d_comm.bcast(retVal->d_iTotalRows,retVal->d_SendSizes.size()-1);

    retVal->d_bFinalized = true;

    return CommunicationList::shared_ptr  ( retVal );
  }

  void  CommunicationList::packReceiveBuffer ( std::vector<double> &send , const Vector &vec ) const
  {
    AMP_ASSERT ( send.size() == d_ReceiveDOFList.size() );
    int *doflist = (int *)getBufferToAvoidDebugVectorCrashing ( d_ReceiveDOFList );
    double *dofval = getBufferToAvoidDebugVectorCrashing ( send );
    vec.getGhostAddValuesByGlobalID ( (int) send.size() , doflist , dofval );
  }

  void  CommunicationList::packSendBuffer ( std::vector<double> &send , const Vector &vec ) const
  {
    AMP_ASSERT ( send.size() == d_SendDOFList.size() );
    int *doflist = (int *)getBufferToAvoidDebugVectorCrashing ( d_SendDOFList );
    double *dofval = getBufferToAvoidDebugVectorCrashing ( send );
    vec.getLocalValuesByGlobalID ( (int) send.size() , doflist , dofval );
  }

  void  CommunicationList::unpackReceiveBufferSet ( const std::vector<double> &recv , Vector &vec ) const
  {
    AMP_ASSERT ( recv.size() == d_ReceiveDOFList.size() );
    int *doflist = (int *)getBufferToAvoidDebugVectorCrashing ( d_ReceiveDOFList );
    const double *dofval = getBufferToAvoidDebugVectorCrashing ( recv );
    vec.setValuesByGlobalID ( (int)recv.size() , doflist , dofval );
  }

  void  CommunicationList::unpackSendBufferAdd ( const std::vector<double> &recv , Vector &vec ) const
  {
    AMP_ASSERT ( recv.size() == d_SendDOFList.size() );
    int *doflist = (int *)getBufferToAvoidDebugVectorCrashing ( d_SendDOFList );
    const double *dofval = getBufferToAvoidDebugVectorCrashing ( recv );
    vec.addLocalValuesByGlobalID ( (int)recv.size() , doflist , dofval );
  }

  void  CommunicationList::unpackSendBufferSet ( const std::vector<double> &recv , Vector &vec ) const
  {
    AMP_ASSERT ( recv.size() == d_SendDOFList.size() );
    int *doflist = (int *)getBufferToAvoidDebugVectorCrashing ( d_SendDOFList );
    const double *dofval = getBufferToAvoidDebugVectorCrashing ( recv );
    vec.setLocalValuesByGlobalID ( (int)recv.size() , doflist , dofval );
  }

  unsigned int  CommunicationList::getLocalGhostID ( int GID ) const
  {
    std::vector<unsigned int>::const_iterator pos =
                        std::lower_bound ( d_ReceiveDOFList.begin() , d_ReceiveDOFList.end() , GID );
    if ( pos == d_ReceiveDOFList.end() )
    {
      std::cout << "GID = " << GID << std::endl;
    }
    AMP_ASSERT ( pos != d_ReceiveDOFList.end() );
    return pos - d_ReceiveDOFList.begin();
  }


void CommunicationList::buildCommunicationArrays ( std::vector<unsigned int>  &DOFs , std::vector<unsigned int> &partitionInfo , int commRank )
{
    d_iBegin = commRank ? partitionInfo[commRank-1] : 0;

    size_t commSize = partitionInfo.size();

    // Check if we are working in serial
    if ( commSize == 1 ) {
        AMP_INSIST(DOFs.size()==0,"Error in communication list, remote DOFs are present for a serial vector");
        d_ReceiveSizes.resize(1,0);
        d_ReceiveDisplacements.resize(1,0);
        d_SendSizes.resize(1,0);
        d_SendDisplacements.resize(1,0);
        d_ReceiveDOFList.resize(0);
        d_SendDOFList.resize(0);
        return;
    }

    // Copy (and sort) the DOFs in d_ReceiveDOFList
    d_ReceiveDOFList = DOFs;
    AMP::Utilities::quicksort(d_ReceiveDOFList);

    // Determine the number of DOFs recieved from each processor (this requires the DOFs to be sorted)
    d_ReceiveSizes.resize( commSize, 0 );
    d_ReceiveDisplacements.resize( commSize, 0 );
    size_t rank = 0;
    size_t start = 0;
    size_t index = 0;
    while ( index < d_ReceiveDOFList.size() ) {
        if (  d_ReceiveDOFList[index] < partitionInfo[rank] ) {
            // Move to the next DOF
            index++;
        } else {
            // Store the number of DOFs with the given rank, and move to the next rank
            d_ReceiveDisplacements[rank] = start;
            d_ReceiveSizes[rank] = index-start;
            start = index;
            rank++;
        }
    }
    d_ReceiveDisplacements[rank] = start;
    d_ReceiveSizes[rank] = index-start;
    AMP_ASSERT(d_ReceiveSizes[commRank]==0);

    // Get the send sizes and displacements
    d_SendSizes.resize( commSize, 0 );
    d_SendDisplacements.resize( commSize, 0 );
    d_comm.allToAll( 1, &(d_ReceiveSizes[0]), &(d_SendSizes[0]) );
    d_SendDisplacements[0] = 0;
    for (size_t i=1; i<commSize; i++) 
        d_SendDisplacements[i] = d_SendDisplacements[i-1] + d_SendSizes[i-1];
    size_t send_buf_size = d_SendDisplacements[commSize-1] + d_SendSizes[commSize-1];

    // Get the send DOFs (requires the recv DOFs to be sorted)
    d_SendDOFList.resize( send_buf_size );
    d_comm.allToAll( getBufferToAvoidDebugVectorCrashing(d_ReceiveDOFList), (int *)&(d_ReceiveSizes[0]), (int *)&(d_ReceiveDisplacements[0]), 
                     getBufferToAvoidDebugVectorCrashing(d_SendDOFList), (int *)&(d_SendSizes[0]) , (int *)&(d_SendDisplacements[0]) , true );

}


void CommunicationList::scatter_set( std::vector<double> &in , std::vector<double> &out ) const
  {
    if ( d_SendSizes.size() == 0 )
      return;
    double *send_buf = getBufferToAvoidDebugVectorCrashing ( in );
    int *send_sizes = (int *)&(d_SendSizes[0]);
    int *send_disps = (int *)&(d_SendDisplacements[0]);
    double *recv_buf = getBufferToAvoidDebugVectorCrashing ( out );
    int *recv_sizes = (int *)&(d_ReceiveSizes[0]);
    int *recv_disps = (int *)&(d_ReceiveDisplacements[0]);

    d_comm.allToAll( send_buf, send_sizes, send_disps, recv_buf, recv_sizes, recv_disps, true );
  }

  void CommunicationList::scatter_add( std::vector<double> &in , std::vector<double> &out ) const
  {
    double *send_buf = getBufferToAvoidDebugVectorCrashing ( in );
    int *send_sizes = (int *)&(d_ReceiveSizes[0]);
    int *send_disps = (int *)&(d_ReceiveDisplacements[0]);
    double *recv_buf = getBufferToAvoidDebugVectorCrashing ( out );
    int *recv_sizes = (int *)&(d_SendSizes[0]);
    int *recv_disps = (int *)&(d_SendDisplacements[0]);

    d_comm.allToAll( send_buf, send_sizes, send_disps, recv_buf, recv_sizes, recv_disps, true );
  }

}
}

