
namespace AMP {
namespace LinearAlgebra {


  inline
  unsigned int CommunicationList::numLocalRows () const
  {
    return d_iNumRows;
  }

  inline
  void  CommunicationList::finalize ( )
  {
    if ( !d_bFinalized )
    {
      d_bFinalized = true;
    }
  }

  inline
  unsigned int  CommunicationList::getTotalSize () const
  {
    return d_iTotalRows;
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

