
namespace AMP {
namespace LinearAlgebra {

  inline
  void CommunicationList::finalizeList ( Castable & )
  {
  }

  inline
  unsigned int CommunicationList::numLocalRows () const
  {
    return d_iNumRows;
  }

  inline
  void  CommunicationList::finalize ( Castable &map )
  {
    if ( !d_bFinalized )
    {
      d_bFinalized = true;
      finalizeList ( map );
    }
  }

  inline
  unsigned int  CommunicationList::getTotalSize () const
  {
    return d_iTotalRows;
  }

  inline
  CommunicationList::CommunicationList ( CommunicationListParameters::shared_ptr params )
    : d_comm ( params->d_comm )
    , d_iNumRows ( params->d_localsize )
    , d_bFinalized ( false )
  {
    d_comm.sumScan((int*)&d_iNumRows,(int*)&d_iTotalRows,1);
    d_iBegin = d_iTotalRows - params->d_localsize;
    int size = d_comm.getSize();
    d_iTotalRows = d_comm.bcast(d_iTotalRows,size-1);
  }

  inline
  CommunicationList::shared_ptr  CommunicationList::createEmpty ( unsigned int local , AMP_MPI c )
  {
    int size = c.getSize();
    CommunicationList  *retVal = new CommunicationList;
    retVal->d_ReceiveSizes.resize ( size );
    retVal->d_ReceiveDisplacements.resize ( size );
    retVal->d_ReceiveDOFList.resize ( 0 );

    retVal->d_SendSizes.resize ( size );
    retVal->d_SendDisplacements.resize ( size );
    retVal->d_SendDOFList.resize ( 0 );

    c.sumScan((int*)&local,(int*)&(retVal->d_iTotalRows),1);
    retVal->d_iBegin = retVal->d_iTotalRows - local;
    retVal->d_comm = c;
    retVal->d_iNumRows = local;
    retVal->d_iTotalRows = c.bcast(retVal->d_iTotalRows,size-1);
    retVal->d_bFinalized = true;

    return CommunicationList::shared_ptr ( retVal );
  }

  inline
  CommunicationList::~CommunicationList ()
  {
  }

  inline
  const std::vector<unsigned int> &CommunicationList::getGhostIDList () const
  {
    return d_ReceiveDOFList;
  }

  inline
  const std::vector<unsigned int> &CommunicationList::getReplicatedIDList () const
  {
    return d_SendDOFList;
  }

  inline
  unsigned int  CommunicationList::getVectorReceiveBufferSize () const
  {
    return d_ReceiveDOFList.size(); 
  }

  inline
  unsigned int  CommunicationList::getVectorSendBufferSize () const
  { 
    return d_SendDOFList.size(); 
  }

  inline
  unsigned int  CommunicationList::getStartGID () const
  { 
    return d_iBegin; 
  }

  inline
  AMP_MPI CommunicationList::getComm() const
  {
    return d_comm;
  }

}
}

