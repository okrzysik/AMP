
namespace AMP {
namespace LinearAlgebra {


inline size_t CommunicationList::numLocalRows() const { return d_iNumRows; }

inline void CommunicationList::finalize()
{
    if ( !d_bFinalized ) {
        d_bFinalized = true;
    }
}

inline size_t CommunicationList::getTotalSize() const { return d_iTotalRows; }


inline CommunicationList::~CommunicationList() {}

inline const std::vector<size_t> &CommunicationList::getGhostIDList() const
{
    return d_ReceiveDOFList;
}

inline const std::vector<size_t> &CommunicationList::getReplicatedIDList() const
{
    return d_SendDOFList;
}

inline size_t CommunicationList::getVectorReceiveBufferSize() const
{
    return d_ReceiveDOFList.size();
}

inline size_t CommunicationList::getVectorSendBufferSize() const { return d_SendDOFList.size(); }

inline size_t CommunicationList::getStartGID() const { return d_iBegin; }

inline const AMP_MPI& CommunicationList::getComm() const { return d_comm; }
}
}
