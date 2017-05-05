#ifndef included_AMP_CommunicationList_h
#define included_AMP_CommunicationList_h

#include "VectorIndexer.h"
#include "utils/AMP_MPI.h"
#include "utils/ParameterBase.h"
#include "utils/shared_ptr.h"
#include <map>
#include <set>
#include <vector>

namespace AMP {
namespace LinearAlgebra {

class VectorData;


/**
 * \class CommunicationListParameters
 * \brief Parameters class encapsulating the data necessary to instantiate a communication list.
 */
class CommunicationListParameters : public ParameterBase
{
public:
    //! Short hand name for the shared pointer for a communication list
    typedef AMP::shared_ptr<CommunicationListParameters> shared_ptr;

    //! The communicator over which the communication list is computed
    AMP_MPI d_comm;

    //! The number of local entities in the vector
    size_t d_localsize;

    //! The remote DOFs that we need to recieve
    std::vector<size_t> d_remote_DOFs;

    //! Default constructor
    CommunicationListParameters();
    //! Copy constructor
    CommunicationListParameters( const CommunicationListParameters & );
};


/**
 * \class CommunicationList
 * \brief What to send where and what to receive from where
 *
 * \details This interface provides communication routines to compute send and receive lists
 * for blocks of data with global indices.  For instance, a vector storing degrees of freedom
 * for nodes in a finite element analysis may share data with other cores in a parallel
 * computation.  This class tracks which local data need to be communicated with other cores
 * and which data should be received from those cores.
 */
class CommunicationList
{
public:
    //! Short hand for shared point to CommunicationList
    typedef AMP::shared_ptr<CommunicationList> shared_ptr;

    /**
     * \brief Construct a communication list
     * \param[in] params  A shared pointer to parameters for constructing the list
     * \details This will set the communicator for the communication list.  It will not
     * compute the communication lists.  Derived classes are expected to call
     * buildCommunicationArrays with appropriate data to compute the communication list
     */
    explicit CommunicationList( CommunicationListParameters::shared_ptr params );

    /**
     * \brief Destroy the communication list
     */
    virtual ~CommunicationList();

    /**
     * \brief Retrieve list of global indices shared locally stored elsewhere
     * \return A vector of indices not owned by the core but are stored locally.
     */
    const std::vector<size_t> &getGhostIDList() const;

    /**
     * \brief Subset a communication list based on a VectorIndexer
     * \param[in] sub  A VectorIndexer pointer that describes a subset
     */
    CommunicationList::shared_ptr subset( VectorIndexer::shared_ptr sub );

    /**
     * \brief Retrieve list of global indices stored here and shared elsewhere
     * \return A vector of indices owned by the core and shared on other cores.
     */
    const std::vector<size_t> &getReplicatedIDList() const;

    /**
     * \brief Retrieve the size of the buffer used to receive data from other processes
     * \details This is an alias of getGhostIDList().size()
     * \return The number of unowned entries on this core
     */
    size_t getVectorReceiveBufferSize() const;

    /**
     * \brief Retrieve the size of the buffer used to send data to other processes
     * \details This is an alias of getReplicatedIDList().size()
     * \return The number of owned entries on this core shared with other cores
     */
    size_t getVectorSendBufferSize() const;

    /**
     * \brief Scatter data stored here to processors that share the data.
     * \param[in] in  A list of data to send to other processors
     * \param[out] out The data received from other processors
     * \details  The convention is if data are set on different processes, then
     * the owner of the data has the correct value.  As such, in a scatter_set,
     * the owner of data scatters the data out which overwrites the data on cores
     * that share the data
     */
    void scatter_set( std::vector<double> &in, std::vector<double> &out ) const;

    /**
     * \brief Scatter data shared here to processors that own the data.
     * \param[in] in  A list of data to send to other processors
     * \param[out] out The data received from other processors
     * \details  When adding data to a vector, any process that shares the data
     * can contribute to the value of the data.  Therefore, this will scatter data
     * that is shared to the core that owns it.  A call to scatter_add is generally
     * followed by a call to scatter_set to ensure parallel consistency.
     */
    void scatter_add( std::vector<double> &in, std::vector<double> &out ) const;

    /**
     * \brief Given an AMP vector, this will create the appropriate send buffer used
     * in an all-to-all.
     * \param[out]  buffer  The output buffer
     * \param[in]  vector  The input vector
     * \details  This will put data owned by this core and shared elsewhere where it needs
     * to be for an all-to-all, given the communication lists computed by
     * buildCommunicationArrays.
     */
    void packSendBuffer( std::vector<double> &buffer, const VectorData &vector ) const;

    /**
     * \brief Given an AMP vector, this will create the appropriate receive buffer used
     * in an all-to-all.
     * \param[out]  buffer  The output buffer
     * \param[in]  vector  The input vector
     * \details  This will put data shared by this core and owned elsewhere where it needs
     * to be for an all-to-all, given the communication lists computed by
     * buildCommunicationArrays.
     */
    void packReceiveBuffer( std::vector<double> &buffer, const VectorData &vector ) const;

    /**
     * \brief Given a buffer from all-to-all, unpack it into a vector
     * \param[in]  buffer  The input buffer
     * \param[in]  vector  The vector to be updated
     * \details  After an all-to-all, this method will unpack the result and set the data
     * in appropriate places in the vector.
     */
    void unpackReceiveBufferSet( const std::vector<double> &buffer, VectorData &vector ) const;

    /**
     * \brief Given a buffer from all-to-all, unpack it into a vector
     * \param[in]  buffer  The input buffer
     * \param[in]  vector  The vector to be updated
     * \details  After an all-to-all, this method will unpack the result and add the data
     * in appropriate places in the vector.
     */
    void unpackSendBufferAdd( const std::vector<double> &buffer, VectorData &vector ) const;

    /**
     * \brief Given a buffer from all-to-all, unpack it into a vector
     * \param[in]  buffer  The input buffer
     * \param[in]  vector  The vector to be updated
     * \details  After an all-to-all, this method will unpack the result and set the data
     * in appropriate places in the vector.
     */
    void unpackSendBufferSet( const std::vector<double> &buffer, VectorData &vector ) const;

    /**
     * \brief  Return the first d.o.f. on this core
     * \return The first d.o.f. on this core
     */
    size_t getStartGID() const;

    /**
     * \brief  Return the total d.o.f. on entire communicator
     */
    size_t getTotalSize() const;

    /**
     * \brief  Return the number of local rows for this communication list
     * \return The number of local d.o.f. for this communication list
     */
    virtual size_t numLocalRows() const;

    /**
     * \brief  Return the local index of a shared datum.
     * \param[in] dof  The global index to get a local ghost id for
     * \details  It is assumed that data are stored in two arrays: an owned array and a shared
     * array.  This function returns the local offset of a shared datum into the shared array
     * \return The index into the shared array for the global index.
     */
    size_t getLocalGhostID( size_t dof ) const;

    /**
     * \brief  Return the communicator used for this communication list
     * \return The communicator.
     */
    const AMP_MPI& getComm() const;


    /**
     * \brief  A call to ensure the communication lists are constructed
     */
    void finalize();


    /**
     * \brief  Construct a CommunicationList with no comunication
     * \param[in]  local  The number of local elements in the vector
     * \param[in]  comm   The AMP_MPI for the vector.
     * \details  Create a communication list with no communication.
     */
    static CommunicationList::shared_ptr createEmpty( size_t local, AMP_MPI comm );

protected:
    /**
     * \brief Construct the arrays that track which data are sent/received to/from where
     * \param[out] dofs List of degrees of freedom not on this processor that will need to be
     * received
     * \param[out] partition  One more than the last degree of freedom on each processor such that
     * \f$\mathit{partition}_{i-1}\le \mathit{d.o.f.} < \mathit{partition}_i\f$.
     * \param[in] commRank  My rank in the communicator
     *
     * \details  Given dofs, the list of needed data for this core, it will compute the send
     * and receive lists for this processor so that scatters can be done in minimal space.
     */
    void buildCommunicationArrays( std::vector<size_t> &dofs,
                                   std::vector<size_t> &partition,
                                   int commRank );

private:
    std::vector<int> d_ReceiveSizes;
    std::vector<int> d_ReceiveDisplacements;
    std::vector<size_t> d_ReceiveDOFList; // Sorted DOF lists

    std::vector<int> d_SendSizes;
    std::vector<int> d_SendDisplacements;
    std::vector<size_t> d_SendDOFList;
    size_t d_iBegin;

    AMP_MPI d_comm;
    size_t d_iNumRows;
    size_t d_iTotalRows;
    bool d_bFinalized;

    CommunicationList();
};
}
}

#include "CommunicationList.inline.h"

#endif
