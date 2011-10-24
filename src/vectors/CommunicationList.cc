#include <vector>
#include "Vector.h"
#include "utils/Utilities.h"
#include "boost/shared_ptr.hpp"

//class Vector;
//typedef  boost::shared_ptr<Vector>   Vector_shared_ptr;

namespace AMP {
namespace LinearAlgebra {

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

    d_ReceiveSizes.resize ( partitionInfo.size() );
    d_ReceiveDisplacements.resize ( partitionInfo.size() );
    d_SendSizes.resize ( partitionInfo.size() );
    d_SendDisplacements.resize ( partitionInfo.size() );
    d_ReceiveDOFList.resize ( 0 );
    d_SendDOFList.resize ( 0 );
    if ( partitionInfo.size() == 1 ) return;  // nothing to do in serial

    unsigned int cur_proc = 0;
    std::vector<unsigned int>::iterator  curDof = DOFs.begin();
    std::vector<unsigned int>::iterator  endDof = DOFs.end();
    while ( curDof != endDof ) {
        std::vector<unsigned int>::iterator  partStart = curDof;
        int firstDof = curDof - DOFs.begin();
        while ( curDof != endDof ) {
            if ( *curDof >= partitionInfo[cur_proc] )
                break;
            curDof++;
        }
        int endDof = curDof - DOFs.begin();
        d_ReceiveDisplacements[cur_proc] = d_ReceiveDOFList.size();
        if ( (int)cur_proc != commRank ) {
            d_ReceiveSizes[cur_proc] = endDof - firstDof;
            d_ReceiveDOFList.insert ( d_ReceiveDOFList.end() , partStart , curDof );
        }
        cur_proc++;
    }

    d_comm.allToAll( 1, &(d_ReceiveSizes[0]), &(d_SendSizes[0]) );

    unsigned int send_buf_size = 0;
    for ( unsigned int i = 0 ; i != partitionInfo.size() ; i++ ) {
        if ( i > 0 )
            d_SendDisplacements[i] = d_SendDisplacements[i-1] + d_SendSizes[i-1];
        send_buf_size += d_SendSizes[i];
    }

    d_SendDOFList.resize ( send_buf_size );
    d_ReceiveSizes[commRank] = 0;
    d_SendDisplacements[0] = 0;

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

