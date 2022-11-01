// This file includes a wrapper class for MPI functions
#ifndef included_AMP_MPI
#define included_AMP_MPI


#include <any>
#include <array>
#include <atomic>
#include <map>
#include <memory>
#include <set>
#include <string>
#include <vector>


// Add the definitions for the TPLs that are used
#include "AMP/AMP_TPLs.h"

// Include MPI if we are building with MPI
#ifdef AMP_USE_MPI
    #include "mpi.h"
#endif

// Define extra comm_world, comm_self, and comm_null ids
#define AMP_COMM_NULL ( (AMP::AMP_MPI::Comm) 0xF4000000 )
#define AMP_COMM_WORLD ( (AMP::AMP_MPI::Comm) 0xF4000001 )
#define AMP_COMM_SELF ( (AMP::AMP_MPI::Comm) 0xF4000002 )


namespace AMP {


/**
 * \class AMP_MPI
 *
 * @brief Provides C++ wrapper around MPI routines.
 *
 * Class AMP_MPI groups common MPI routines into one globally-accessible location.
 * It provides small, simple routines that are common in MPI code.
 * In some cases, the calling syntax has been simplified for convenience.
 * Moreover, there is no reason to include the preprocessor ifdef/endif guards around
 * these calls, since the MPI libraries are not called in these routines if the MPI
 * libraries are not being used (e.g., when writing serial code).
 * Note: Many of the communication routines are templated on type.
 * When using unknown types the reduce calls will fail, the send and gather calls
 * should succeed provided that the size of the data type object is a fixed size on
 * all processors.  sizeof(type) must be the same for all elements and processors.
 */
class alignas( 8 ) AMP_MPI final
{
public:
    enum class ThreadSupport : int { SINGLE, FUNNELED, SERIALIZED, MULTIPLE };

#ifdef AMP_USE_MPI
    typedef MPI_Comm Comm;
    typedef MPI_Datatype Datatype;
    typedef MPI_Request Request2;
#else
    typedef uint32_t Comm;
    typedef uint32_t Datatype;
    typedef uint32_t Request2;
#endif

    class Request final
    {
    public:
        Request( Request2 request = Request2(), std::any data = std::any() );
        ~Request();
        operator Request2() const { return d_data->first; }
        Request2 *get() { return &d_data->first; }

    private:
        std::shared_ptr<std::pair<Request2, std::any>> d_data;
    };


public: // Constructors
    /**
     *\brief  Empty constructor
     *\details  This creates an empty constructor that does not contain an MPI communicator.
     */
    AMP_MPI();


    //!  Empty destructor
    ~AMP_MPI();


    /**
     * \brief Constructor from existing MPI communicator
     * \details  This constructor creates a new communicator from an existing MPI communicator.
     *    This does not create a new internal MPI_Comm, but uses the existing comm.
     *    Note that by default, this will not free the MPI_Comm object and the user is
     * responsible
     *      for free'ing the MPI_Comm when it is no longer used.  This behavior is controlled by the
     *      optional manage argument.
     * \param[in] comm      Existing MPI communicator
     * \param[in] manage    Do we want to manage the comm (free the MPI_Comm when this object leaves
     * scope)
     */
    AMP_MPI( Comm comm, bool manage = false );


    /**
     * \brief Constructor from existing communicator
     * \details  This constructor creates a new communicator from an existing communicator.
     *   This does not create a new internal Comm, but uses the existing comm.
     * \param[in] comm Existing communicator
     */
    AMP_MPI( const AMP_MPI &comm );


    /*!
     * Move constructor
     * @param[in] rhs       Communicator to copy
     */
    AMP_MPI( AMP_MPI &&rhs );


    /**
     * \brief Assignment operator
     * \details  This operator overloads the assignment to correctly copy an communicator
     * \param[in] comm      Existing MPI object
     */
    AMP_MPI &operator=( const AMP_MPI &comm );


    /*!
     * Move assignment operator
     * @param[in] rhs       Communicator to copy
     */
    AMP_MPI &operator=( AMP_MPI &&rhs );


    /**
     * \brief Reset the object
     * \details  This resets the object to the empty state without an MPI_Comm
     */
    void reset();


public: // Member functions
    /**
     * \brief Get the node name
     * \details  This function returns a unique name for each node.
     *    It is a wrapper for MPI_Get_processor_name.
     */
    static std::string getNodeName();


    //! Function to return the number of processors available
    static int getNumberOfProcessors();


    //! Function to return the affinity of the current process
    static std::vector<int> getProcessAffinity();


    //! Function to set the affinity of the current process
    static void setProcessAffinity( const std::vector<int> &procs );


    /**
     * \brief Load balance the processes within a node
     * \details  This function will redistribute the processes within a node using the
     *    process affinities to achieve the desired load balance.
     *    Note: this is a global operation on the given comm, and it is STRONGLY
     *    recommended to use COMM_WORLD.
     * \param[in] comm      The communicator to use (Default is COMM_WORLD)
     * \param[in] method    The desired load balance method to use:
     *                      1:  Adjust the affinities so all processes share the given processors.
     *                          This effectively allows the OS to handle the load balancing
     *                          by migrating the processes as necessary.  This is recommended
     *                          for most users and use cases. (default)
     *                      2:  Adjust the affinities so that the fewest number of processes
     *                          overlap. This will try to give each process a unique set of
     *                          processors while ensuring that each process has at least N_min
     *                          processes.
     * \param[in] procs     An optional list of processors to use.
     *                      By default, setting this to an empty vector will use all available
     *                      processors on the given node. \param[in] N_min
     *                      The minimum number of processors for any process
     *                      (-1 indicates all available processors).
     * \param[in] N_max     The maximum number of processors for any process
     *                      (-1 indicates all available processors).
     *
     */
    static void balanceProcesses( const AMP_MPI &comm           = AMP_MPI( AMP_COMM_WORLD ),
                                  int method                    = 1,
                                  const std::vector<int> &procs = std::vector<int>(),
                                  int N_min                     = 1,
                                  int N_max                     = -1 );


    //! Query the level of thread support
    static ThreadSupport queryThreadSupport();


    /**
     * \brief Generate a random number
     * \details  This generates a random number that is consistent across the comm
     */
    size_t rand() const;


    /**
     * \brief Split an existing communicator
     * \details  This creates a new communicator by splitting an existing communicator.
     *   See MPI_Comm_split for information on how the underlying split will occur.
     *   Note: the underlying MPI_Comm object will be free'd automatically when it is no longer
     *   used by any MPI objects.
     * \param[in] color     Control of subset assignment (nonnegative integer).
     *                      Processes with the same color are in the same new communicator.
     *                      -1: processor will not be a member of any object
     *                          (NULL object will be returned)
     * \param[in] key       Control of rank assignment (integer).
     *                      Note that, for a fixed color, the keys need not be unique.
     *                      The processes will be sorted in ascending order according to this key,
     *                      then all the processes in a given color will have the relative rank
     *                      order as they did in their parent group. (See MPI_Comm_split)
     * \param[in] manage    Do we want to manage the comm
     *                      (free the MPI_Comm when this object leaves scope)
     */
    AMP_MPI split( int color, int key = -1, bool manage = true ) const;


    /**
     * \brief Split an existing communicator by node
     * \details  This creates a new communicator by splitting an existing communicator
     *   by the node.  This will result in a separate MPI_Comm for each physical node.
     *   Internally this will use MPI_Get_processor_name to identify the nodes.
     *   Note: the underlying MPI_Comm object will be free'd automatically when it is no longer
     *   used by any MPI objects)
     * \param[in] key       Control of rank assignment (integer).
     *                      Note that, for a fixed color, the keys need not be unique.
     *                      The processes will be sorted in ascending order according to this key,
     *                      then all the processes in a given color will have the relative rank
     *                      order as they did in their parent group. (See MPI_Comm_split)
     * \param[in] manage    Do we want to manage the comm
     *                      (free the MPI_Comm when this object leaves scope)
     */
    AMP_MPI splitByNode( int key = -1, bool manage = true ) const;


    /**
     * \brief Duplicate an existing communicator
     * \details  This creates a new communicator by duplicating an existing communicator.
     *   The resulting communicator will exist over the same processes, but have a different
     * context.
     *   Note: the underlying MPI_Comm object will be free'd automatically when it is no longer
     *   used by any MPI objects.
     * \param[in] manage    Do we want to manage the comm
     *                      (free the MPI_Comm when this object leaves scope)
     */
    AMP_MPI dup( bool manage = true ) const;


    /**
     * \brief Create a communicator from the intersection of two communicators
     * \details  This creates a new communicator by intersecting two existing communicators.
     *   Any processors that do not contain the both communicators will receive a NULL communicator.
     *   There are 3 possible cases:
     *      The communicators are disjoint (a null communicator will be returned on all processors).
     *      One communicator is a sub communicator of another.  This will require communication on
     *          the smaller communicator only.
     *      The communicators partially overlap.  This will require communication on the first
     * communicator.
     * \param[in] comm1     First communicator
     * \param[in] comm2     First communicator
     */
    static AMP_MPI intersect( const AMP_MPI &comm1, const AMP_MPI &comm2 );


    /**
     * Check if the current communicator is NULL
     */
    inline bool isNull() const { return d_isNull; }


    /**
     * \brief Return the global ranks for the comm
     * \details  This returns a vector which contains the global ranks for each
     *   member of the communicator.  The global ranks are defined according to WORLD comm.
     */
    std::vector<int> globalRanks() const;


    /**
     *  Get the current MPI communicator.
     *  Note: The underlying MPI_Comm object may be free'd by the object when it is no
     *  longer used by any communicators.  If the user has made a copy using the
     *  getCommunicator routine, then it may be free'd without user knowledge.  The
     *  user is responsible for checking if the communicator is valid, or keeping a
     *  copy of the communicator that provided the MPI_Communicator.
     */
    inline const Comm &getCommunicator() const { return d_comm; }


    /**
     * \brief Overload operator ==
     * \details  Overload operator comm1 == comm2.  Two MPI objects are == if they share the same
     * communicator.
     *   Note: this is a local operation.
     */
    bool operator==( const AMP_MPI & ) const;


    /**
     * \brief Overload operator !=
     * \details  Overload operator comm1 != comm2.  Two MPI objects are != if they
     *   do not share the same communicator.
     *   Note: this is a local operation.
     */
    bool operator!=( const AMP_MPI & ) const;


    /**
     * \brief Overload operator <
     * \details  Overload operator comm1 < comm2.  One MPI object is < another iff all the
     *   processors in the first object are also in the second.  Additionally, the second
     *   object must contain at least one processor that is not in the first object.
     *   This is a collective operation, based on the first communicator.
     *   As a result all processors on the first communicator will return the same value,
     *   while any processors that are not on the first communicator will return an unknown value.
     *   Additionally, all processors on the first object MUST call this routine and will be
     *   synchronized through this call (there is an internalallReduce).
     */
    bool operator<( const AMP_MPI & ) const;


    /**
     * \brief Overload operator <=
     * \details  Overload operator comm1 <= comm2.  One MPI object is <= another iff all the
     *   processors in the first object are also in the second.  This is a collective operation,
     *   based on the first communicator.  As a result all processors on the first communicator
     *   will return the same value, while any processors that are not on the first communicator
     *   will return an unknown value.  Additionally, all processors on the first object MUST
     *   call this routine and will be synchronized through this call (there is an internal
     *   allReduce).
     */
    bool operator<=( const AMP_MPI & ) const;


    /**
     * \brief Overload operator >
     * \details  Overload operator comm1 > comm2.  One MPI object is > another iff all the
     *   processors in the second object are also in the first.  Additionally, the first object
     *   must contain at least one processor that is not in the second object.
     *   This is a collective operation, based on the first communicator.
     *   As a result all processors on the first communicator will return the same value,
     *   while any processors that are not on the first communicator will return an unknown value.
     *   Additionally, all processors on the first object MUST call this routine and will be
     *   synchronized through this call (there is an internal allReduce).
     */
    bool operator>( const AMP_MPI & ) const;


    /**
     * \brief Overload operator >=
     * \details  Overload operator comm1 >= comm2.  One MPI object is > another iff all the
     *   processors in the second object are also in the first.  Additionally, the first object
     *   must contain at least one processor that is not in the second object.
     *   This is a collective operation, based on the first communicator.
     *   As a result all processors on the first communicator will return the same value, while any
     *   processors that are not on the first communicator will return an unknown value.
     *   Additionally, all processors on the first object MUST call this routine and will be
     *   synchronized through this call (there is an internal allReduce).
     */
    bool operator>=( const AMP_MPI & ) const;


    /**
     * \brief Compare to another communicator
     * \details  This compares the current communicator to another communicator.
     *   This returns 1 if the two communicators are equal (they share the same MPI communicator),
     *   2 if the contexts and groups are the same, 3 if different contexts but identical groups,
     *   4 if different contexts but similar groups, and 0 otherwise.
     *   Note: this is a local operation.
     */
    int compare( const AMP_MPI & ) const;


    /**
     * Return the processor rank (identifier) from 0 through the number of
     * processors minus one.
     */
    inline int getRank() const { return d_rank; }


    /**
     * Return the number of processors.
     */
    inline int getSize() const { return d_size; }


    /**
     * Return the maximum tag
     */
    inline int maxTag() const { return d_maxTag; }


    /**
     * \brief   Return a new tag
     * \details This routine will return an unused tag for communication.
     *   Note that this tag may match a user tag, but this function will
     *   not return two duplicate tags.  This is a global operation.
     */
    int newTag();


    /**
     * Call MPI_Abort or exit depending on whether running with one or more
     * processes and value set by function above, if called.  The default is
     * to call exit(-1) if running with one processor and to call MPI_Abort()
     * otherwise.  This function avoids having to guard abort calls in
     * application code.
     */
    void abort() const;


    /**
     * Set boolean flag indicating whether exit or abort is called when running
     * with one processor.  Calling this function influences the behavior of
     * calls to abort().  By default, the flag is true meaning that
     * abort() will be called.  Passing false means exit(-1) will be called.
     */
    void setCallAbortInSerialInsteadOfExit( bool flag = true );


    /**
     * \brief   Boolean all reduce
     * \details This function performs a boolean all reduce across all processors.
     *   It returns true iff all processor are true;
     * \param[in] value     The input value for the all reduce
     */
    bool allReduce( const bool value ) const;


    /**
     * \brief   Boolean any reduce
     * \details This function performs a boolean any reduce across all processors.
     *   It returns true if any processor is true;
     * \param[in] value     The input value for the all reduce
     */
    bool anyReduce( const bool value ) const;


    /**
     * \brief   Boolean all reduce
     * \details This function performs a boolean all reduce across all processors.
     *   It returns true iff all processor are true;
     * \param[in] value     The input value for the all reduce
     */
    void allReduce( std::vector<bool> &value ) const;


    /**
     * \brief   Boolean any reduce
     * \details This function performs a boolean any reduce across all processors.
     *   It returns true if any processor is true;
     * \param[in] value     The input value for the all reduce
     */
    void anyReduce( std::vector<bool> &value ) const;


    /**
     * \brief   Sum Reduce
     * \details This function performs a sum all reduce across all processor.
     *   It returns the sum across all processors;
     * \param[in] value     The input value for the all reduce
     */
    template<class type>
    type sumReduce( const type value ) const;


    /**
     * \brief   Sum Reduce
     * \details Perform an array sum Reduce across all nodes.  Each
     * processor contributes an array of values, and the
     * element-wise sum is returned in the same array.
     * \param[in] x         The input/output array for the reduce
     * \param[in] n         The number of values in the array (must match on all nodes)
     */
    template<class type>
    void sumReduce( type *x, int n = 1 ) const;


    /**
     * \brief   Sum Reduce
     * \details Perform an array sum Reduce across all nodes.  Each
     * processor contributes an array of values, and the
     * element-wise sum is returned in the same array.
     * \param[in] x         The input array for the reduce
     * \param[in] y         The output array for the reduce
     * \param[in] n         The number of values in the array (must match on all nodes)
     */
    template<class type>
    void sumReduce( const type *x, type *y, int n = 1 ) const;


    /**
     * \brief   Min Reduce
     * \details This function performs a min all reduce across all processor.
     *   It returns the minimum value across all processors;
     * \param[in] value     The input value for the all reduce
     */
    template<class type>
    type minReduce( const type value ) const;


    /**
     * \brief   Min Reduce
     * \details Perform an array min Reduce across all nodes.  Each
     * processor contributes an array of values, and the
     * element-wise minimum is returned in the same array.
     * \param[in] x         The input/output array for the reduce
     * \param[in] n         The number of values in the array (must match on all nodes)
     */
    template<class type>
    void minReduce( type *x, int n ) const;


    /**
     * \brief   Min Reduce
     * \details Perform an array min Reduce across all nodes.  Each
     * processor contributes an array of values, and the
     * element-wise minimum is returned in the same array.
     * \param[in] x         The input array for the reduce
     * \param[in] y         The output array for the reduce
     * \param[in] n         The number of values in the array (must match on all nodes)
     */
    template<class type>
    void minReduce( const type *x, type *y, int n ) const;

    /**
     * \brief   Min Reduce
     * \details Perform an array min Reduce across all nodes.  Each
     * processor contributes an array of values, and the
     * element-wise minimum is returned in the same array.
     * \param[in] x         The input/output array for the reduce
     * \param[in] n         The number of values in the array (must match on all nodes)
     * \param[out] rank     Output array indicating the rank of the processor containing
     *                      the minimum value
     */
    template<class type>
    void minReduce( type *x, int n, int *rank ) const;


    /**
     * \brief   Sum Reduce
     * \details Perform an array min Reduce across all nodes.  Each
     * processor contributes an array of values, and the
     * element-wise minimum is returned in the same array.
     * \param[in] x         The input array for the reduce
     * \param[in] y         The output array for the reduce
     * \param[in] n         The number of values in the array (must match on all nodes)
     * \param[out] rank     Output array indicating the rank of the processor containing
     *                      the minimum value
     */
    template<class type>
    void minReduce( const type *x, type *y, int n, int *rank ) const;


    /**
     * \brief   Max Reduce
     * \details This function performs a max all reduce across all processor.
     *   It returns the maximum value across all processors;
     * \param[in] value     The input value for the all reduce
     */
    template<class type>
    type maxReduce( const type value ) const;


    /**
     * \brief   Max Reduce
     * \details Perform an array max Reduce across all nodes.  Each
     * processor contributes an array of values, and the
     * element-wise maximum is returned in the same array.
     *
     * If a 'rank_of_min' argument is provided, it will set the array to the
     * rank of process holding the minimum value.  Like the double argument,
     * the size of the supplied 'rank_of_min' array should be n.
     * \param[in] x         The input/output array for the reduce
     * \param[in] n         The number of values in the array (must match on all nodes)
     */
    template<class type>
    void maxReduce( type *x, int n ) const;


    /**
     * \brief   Max Reduce
     * \details Perform an array max Reduce across all nodes.  Each
     * processor contributes an array of values, and the
     * element-wise maximum is returned in the same array.
     *
     * If a 'rank_of_min' argument is provided, it will set the array to the
     * rank of process holding the minimum value.  Like the double argument,
     * the size of the supplied 'rank_of_min' array should be n.
     * \param[in] x         The input array for the reduce
     * \param[in] y         The output array for the reduce
     * \param[in] n         The number of values in the array (must match on all nodes)
     */
    template<class type>
    void maxReduce( const type *x, type *y, int n ) const;

    /**
     * \brief   Max Reduce
     * \details Perform an array max Reduce across all nodes.  Each
     * processor contributes an array of values, and the
     * element-wise maximum is returned in the same array.
     *
     * If a 'rank_of_min' argument is provided, it will set the array to the
     * rank of process holding the minimum value.  Like the double argument,
     * the size of the supplied 'rank_of_min' array should be n.
     * \param[in] x         The input/output array for the reduce
     * \param[in] n         The number of values in the array (must match on all nodes)
     * \param[out] rank     Output array indicating the rank of the processor containing
     *                      the minimum value
     */
    template<class type>
    void maxReduce( type *x, int n, int *rank ) const;


    /**
     * \brief   Max Reduce
     * \details Perform an array max Reduce across all nodes.  Each
     * processor contributes an array of values, and the
     * element-wise maximum is returned in the same array.
     *
     * If a 'rank_of_min' argument is provided, it will set the array to the
     * rank of process holding the minimum value.  Like the double argument,
     * the size of the supplied 'rank_of_min' array should be n.
     * \param[in] x         The input array for the reduce
     * \param[in] y         The output array for the reduce
     * \param[in] n         The number of values in the array (must match on all nodes)
     * \param[out] rank     Output array indicating the rank of the processor containing
     *                      the minimum value
     */
    template<class type>
    void maxReduce( const type *x, type *y, int n, int *rank ) const;


    /**
     * \brief    Scan Sum Reduce
     * \details  Computes the sum scan (partial reductions) of data on a collection of processes.
     *   See MPI_Scan for more information.
     * \param[in] x         The input value for the scan
     */
    template<class type>
    type sumScan( const type &x ) const;


    /**
     * \brief    Scan Sum Reduce
     * \details  Computes the sum scan (partial reductions) of data on a collection of processes.
     *   See MPI_Scan for more information.
     * \param[in] x         The input array for the scan
     * \param[in] y         The output array for the scan
     * \param[in] n         The number of values in the array (must match on all nodes)
     */
    template<class type>
    void sumScan( const type *x, type *y, int n ) const;


    /**
     * \brief    Scan Min Reduce
     * \details  Computes the min scan (partial reductions) of data on a collection of processes.
     *   See MPI_Scan for more information.
     * \param[in] x         The input value for the scan
     */
    template<class type>
    type minScan( const type &x ) const;


    /**
     * \brief    Scan Min Reduce
     * \details  Computes the min scan (partial reductions) of data on a collection of processes.
     *   See MPI_Scan for more information.
     * \param[in] x         The input array for the scan
     * \param[in] y         The output array for the scan
     * \param[in] n         The number of values in the array (must match on all nodes)
     */
    template<class type>
    void minScan( const type *x, type *y, int n ) const;


    /**
     * \brief    Scan Max Reduce
     * \details  Computes the max scan (partial reductions) of data on a collection of processes.
     *   See MPI_Scan for more information.
     * \param[in] x         The input value for the scan
     */
    template<class type>
    type maxScan( const type &x ) const;


    /**
     * \brief    Scan Max Reduce
     * \details  Computes the max scan (partial reductions) of data on a collection of processes.
     *   See MPI_Scan for more information.
     * \param[in] x         The input array for the scan
     * \param[in] y         The output array for the scan
     * \param[in] n     The number of values in the array (must match on all nodes)
     */
    template<class type>
    void maxScan( const type *x, type *y, int n ) const;


    /**
     * \brief   Broadcast
     * \details This function broadcasts a value from root to all processors
     * \param[in] value     The input value for the broadcast.
     * \param[in] root      The processor performing the broadcast
     */
    template<class type>
    type bcast( const type &value, int root ) const;


    /**
     * \brief   Broadcast
     * \details This function broadcasts an array from root to all processors
     * \param[in] value     The input/output array for the broadcast
     * \param[in] n         The number of values in the array (must match on all nodes)
     * \param[in] root      The processor performing the broadcast
     */
    template<class type>
    void bcast( type *value, int n, int root ) const;


    /**
     * Perform a global barrier across all processors.
     */
    void barrier() const;


    /**
     *
     */
    /**
     * \brief   Perform a global barrier putting idle processors to sleep
     * \details  This function uses an MPI_Ibarrier to start the barrier then
     *    waits for completion using sleep.  Note: this version will sleep in
     *    increments of 10 ms and should not be used in performance critical
     *    regions.  It's primary use is to allow the CPU to go idle if unused
     *    for a long time.
     */
    void sleepBarrier() const;


    /*!
     * @brief This function sends an MPI message with an array to another processor.
     *
     * If the receiving processor knows in advance the length
     * of the array, use "send_length = false;"  otherwise,
     * this processor will first send the length of the array,
     * then send the data.  This call must be paired with a
     * matching call to recv.
     *
     * @param[in] data      Data to send
     * @param[in] recv      Receiving processor number.
     * @param[in] tag       Optional integer argument specifying an integer tag
     *                      to be sent with this message.  Default tag is 0.
     *                      The matching recv must share this tag.
     */
    template<class type>
    void send( const type &data, int recv, int tag = 0 ) const;


    /*!
     * @brief This function sends an MPI message with an array to another processor.
     *
     * If the receiving processor knows in advance the length
     * of the array, use "send_length = false;"  otherwise,
     * this processor will first send the length of the array,
     * then send the data.  This call must be paired with a
     * matching call to recv.
     *
     * @param[in] buf       Pointer to array buffer with length integers.
     * @param[in] length    Number of integers in buf that we want to send.
     * @param[in] recv      Receiving processor number.
     * @param[in] tag       Optional integer argument specifying an integer tag
     *                      to be sent with this message.  Default tag is 0.
     *                      The matching recv must share this tag.
     */
    template<class type>
    void send( const type *buf, int length, int recv, int tag = 0 ) const;


    /*!
     * @brief This function sends an MPI message with an array of bytes
     * (MPI_BYTES) to receiving_proc_number.
     *
     * This call must be paired with a matching call to recvBytes.
     *
     * @param[in] buf       Void pointer to an array of number_bytes bytes to send.
     * @param[in] N_bytes   Integer number of bytes to send.
     * @param[in] recv      Receiving processor number.
     * @param[in] tag       Optional integer argument specifying an integer tag
     *                      to be sent with this message.  Default tag is 0.
     *                      The matching recv must share this tag.
     */
    void sendBytes( const void *buf, int N_bytes, int recv, int tag = 0 ) const;


    /*!
     * @brief This function sends an MPI message with an array
     *   to another processor using a non-blocking call.
     *   The receiving processor must know the length of the array.
     *   This call must be paired  with a matching call to Irecv.
     *
     * @param[in] buf       Data to send
     * @param[in] recv_proc Receiving processor number.
     * @param[in] tag       Integer argument specifying an integer tag
     *                      to be sent with this message.
     */
    template<class type>
    Request Isend( const type &data, int recv_proc, int tag ) const;


    /*!
     * @brief This function sends an MPI message with an array
     *   to another processor using a non-blocking call.
     *   The receiving processor must know the length of the array.
     *   This call must be paired  with a matching call to Irecv.
     *
     * @param[in] buf       Pointer to array buffer with length integers.
     * @param[in] length    Number of integers in buf that we want to send.
     * @param[in] recv_proc Receiving processor number.
     * @param[in] tag       Integer argument specifying an integer tag
     *                      to be sent with this message.
     */
    template<class type>
    Request Isend( const type *buf, int length, int recv_proc, int tag ) const;


    /*!
     * @brief This function sends an MPI message with an array of bytes
     *   (MPI_BYTES) to receiving_proc_number using a non-blocking call.
     *   The receiving processor must know the number of bytes to receive.
     *   This call must be paired with a matching call to IrecvBytes.
     *
     * @param[in] buf       Void pointer to an array of number_bytes bytes to send.
     * @param[in] N_bytes   Integer number of bytes to send.
     * @param[in] recv_proc Receiving processor number.
     * @param[in] tag       Integer argument specifying an integer tag
     *                  to be sent with this message.
     */
    Request IsendBytes( const void *buf, int N_bytes, int recv_proc, int tag ) const;


    /*!
     * @brief This function receives an MPI message with a data array from another processor.
     *    This call must be paired with a matching call to send.
     *
     * @param[in] send      Processor number of sender.
     * @param[in] tag       Optional integer argument specifying a tag which must be matched
     *                      by the tag of the incoming message. Default tag is 0.
     */
    template<class type>
    type recv( int send, int tag = 0 ) const;


    /*!
     * @brief This function receives an MPI message with a data array from another processor.
     *    This call must be paired with a matching call to send.
     *
     * @param[in] buf       Pointer to integer array buffer with capacity of length integers.
     * @param[in] length    The number of elements to be received.
     * @param[in] send      Processor number of sender.
     * @param[in] tag       Optional integer argument specifying a tag which must be matched
     *                      by the tag of the incoming message. Default tag is 0.
     */
    template<class type>
    void recv( type *buf, int length, int send, int tag = 0 ) const;


    /*!
     * @brief This function receives an MPI message with a data
     * array from another processor.
     *
     * If this processor knows in advance the length of the array,
     * use "get_length = false;" otherwise we will get the return size.
     * This call must be paired with a matching call to send.
     *
     * @param[in] buf       Pointer to integer array buffer with capacity of length integers.
     * @param[in] length    If get_length==true: The number of elements to be received, otherwise
     *                      the maximum number of values that can be stored in buf.
     *                      On output the number of received elements.
     * @param[in] send      Processor number of sender.
     * @param[in] get_length Optional boolean argument specifying if we first
     *                      need to check the message size to get the size of the array.
     *                      Default value is false.
     * @param[in] tag       Optional integer argument specifying a tag which must be matched
     *                      by the tag of the incoming message. Default tag is 0.
     */
    template<class type>
    void recv( type *buf, int &length, int send, bool get_length, int tag = 0 ) const;


    /*!
     * @brief This function receives an MPI message with an array of
     * max size number_bytes (MPI_BYTES) from any processor.
     *
     * This call must be paired with a matching call to sendBytes.
     *
     * @param[in] buf       Void pointer to a buffer of size number_bytes bytes.
     * @param[in] N_bytes   Integer number specifying size of buf in bytes.
     * @param[in] send      Integer number specifying size of buf in bytes.
     * @param[in] tag       Optional integer argument specifying a tag which
     *   must be matched by the tag of the incoming message. Default
     *   tag is 0.
     */
    void recvBytes( void *buf, int N_bytes, int send, int tag = 0 ) const;


    /*!
     * @brief This function receives an MPI message with a data
     * array from another processor using a non-blocking call.
     *
     * @param[in] data       Data to receive
     * @param[in] send_proc  Processor number of sender.
     * @param[in] tag        Optional integer argument specifying a tag which must
     *                      be matched by the tag of the incoming message.
     */
    template<class type>
    Request Irecv( type &data, int send_proc, int tag ) const;


    /*!
     * @brief This function receives an MPI message with a data
     * array from another processor using a non-blocking call.
     *
     * @param[in] buf        Recieve buffer
     * @param[in] length     Maximum number of values that can be stored in buf.
     * @param[in] send_proc  Processor number of sender.
     * @param[in] tag        Optional integer argument specifying a tag which must
     *                      be matched by the tag of the incoming message.
     */
    template<class type>
    Request Irecv( type *buf, int length, int send_proc, int tag ) const;


    /*!
     * @brief This function receives an MPI message with an array of
     * max size number_bytes (MPI_BYTES) from any processor.
     *
     * This call must be paired with a matching call to sendBytes.
     *
     * @param[in] buf       Void pointer to a buffer of size number_bytes bytes.
     * @param[in] N_bytes   Integer number specifying size of buf in bytes.
     * @param[in] send_proc Processor number of sender.
     * @param[in] tag       Integer argument specifying a tag which must
     *                      be matched by the tag of the incoming message.
     */
    Request IrecvBytes( void *buf, int N_bytes, int send_proc, int tag ) const;


    /*!
     * @brief This function sends and recieves data using a blocking call
     *
     * @param[in] sendbuf   Initial address of send buffer (choice).
     * @param[in] sendcount Number of elements to send (integer).
     * @param[in] dest      Rank of destination (integer).
     * @param[in] sendtag   Send tag (integer).
     * @param[out] recvbuf  Initial address of recv buffer (choice).
     * @param[in] recvcount Maximum number of elements to receive (integer).
     * @param[in] source    Rank of source (integer).
     * @param[in] recvtag   Receive tag (integer).
     */
    template<class type>
    void sendrecv( const type *sendbuf,
                   int sendcount,
                   int dest,
                   int sendtag,
                   type *recvbuf,
                   int recvcount,
                   int source,
                   int recvtag ) const;


    /*!
     * Each processor sends every other processor a single value.
     * @param[in] x      Input value for allGather
     * @return           Output array for allGather
     */
    template<class type>
    std::vector<type> allGather( const type &x ) const;


    /*!
     * Each processor sends every other processor an array
     * @param[in] x_in   Input array for allGather
     * @return           Output array for allGather
     */
    template<class type>
    std::vector<type> allGather( const std::vector<type> &x_in ) const;


    /*!
     * Each processor sends every other processor a single value.
     * The x_out array should be preallocated to a length equal
     * to the number of processors.
     * @param[in] x_in      Input value for allGather
     * @param[in] x_out     Output array for allGather (must be preallocated to the size of the
     * communicator)
     */
    template<class type>
    void allGather( const type &x_in, type *x_out ) const;


    /*!
     * Each processor sends an array of data to all other processors.
     * Each processor receives the values from all processors and gathers them
     * to a single array.  If successful, the total number of received
     * elements will be returned.
     * @param[in] send_data     Input array
     * @param[in] send_cnt      The number of values to send
     * @param[in] recv_data     Output array of received values
     * @param[in] recv_cnt      The number of values to receive from each processor (N).
     *                          If known, this should be provided as an input.  Otherwise
     *                          it is an optional output that will return the number of
     *                          received values from each processor.
     * @param[in] recv_disp     The displacement (relative to the start of the array)
     *                          from which to store the data received from processor i.
     *                          If known, this should be provided as an input.  Otherwise
     *                          it is an optional output that will return the starting location
     *                          (relative to the start of the array) for the received data from
     *                          processor i.
     * @param[in] known_recv    Are the received counts and displacements known.
     *                          If the received sizes are known, then they must be provided,
     *                          and an extra communication step is not necessary.  If the received
     *                          sizes are not known, then an extra communication step will occur
     *                          internally and the sizes and displacements will be returned (if
     * desired).
     */
    template<class type>
    int allGather( const type *send_data,
                   int send_cnt,
                   type *recv_data,
                   int *recv_cnt   = nullptr,
                   int *recv_disp  = nullptr,
                   bool known_recv = false ) const;


    /*!
     * This function combines sets from different processors to create a single master set
     * @param[in] set       Input/Output std::set for the gather.
     */
    template<class type>
    void setGather( std::set<type> &set ) const;


    /*!
     * This function combines std::maps from different processors to create a single master std::map
     * If two or more ranks share the same key, the lowest rank will be used
     * @param[in] map       Input/Output std::map for the gather.
     */
    template<class KEY, class DATA>
    void mapGather( std::map<KEY, DATA> &map ) const;


    /*!
     * Each processor sends a value to root
     * @param[in] x      Input value to send
     * @param[in] root   The processor receiving the data
     * @return           Output array for gather (empty if not root)
     */
    template<class type>
    std::vector<type> gather( const type &x, int root ) const;


    /*!
     * Each processor sends every other processor an array
     * @param[in] x      Input array to send
     * @param[in] root   The processor receiving the data
     * @return           Output array for gather (empty if not root)
     */
    template<class type>
    std::vector<type> gather( const std::vector<type> &x, int root ) const;


    /*!
     * Each processor sends multiple values to root
     * @param[in] send_data     Input array
     * @param[in] send_cnt      The number of values to send


     * @param[in] send_data     Input array
     * @param[in] send_cnt      The number of values to send
     * @param[in] recv_data     Output array of received values
     * @param[in] recv_cnt      The number of values to receive from each processor (N).
     *                          If known, this should be provided as an input.
     * @param[in] recv_disp     The displacement (relative to the start of the array)
     *                          from which to store the data received from processor i.
     *                          If known, this should be provided as an input.
     */
    template<class type>
    void gather( const type *sendbuf,
                 int sendcount,
                 type *recvbuf,
                 const int *recv_cnt,
                 const int *recv_disp,
                 int root ) const;


    /*!
     * Each processor sends an array of n values to each processor.
     * Each processor sends an array of n values to each processor.
     * The jth block of data is sent from processor i to processor j and placed
     * in the ith block on the receiving processor.  In the variable
     * description, N is the size of the communicator.  Note that this is a
     * blocking global communication.
     * @param[in] n             The number of elements in each data block to send.
     * @param[in] send_data     Input array (nxN)
     * @param[in] recv_data     Output array of received values (nxN)
     */
    template<class type>
    void allToAll( int n, const type *send_data, type *recv_data ) const;


    /*!
     * Each processor sends an array of data to the different processors.
     * Each processor may send any size array to any processor.  In the variable
     * description, N is the size of the communicator.  Note that this is a
     * blocking global communication.  If successful, the total number of received
     * elements will be returned.
     * @param[in] send_data     Input array
     * @param[in] send_cnt      The number of values to send to each processor (N)
     * @param[in] send_disp     The displacement (relative to the start of the array)
     *                          from which to send to processor i
     * @param[in] recv_data     Output array of received values
     * @param[in] recv_cnt      The number of values to receive from each processor (N).
     *                          If known, this should be provided as an input.  Otherwise
     *                          it is an optional output that will return the number of
     *                          received values from each processor.
     * @param[in] recv_disp     The displacement (relative to the start of the array)
     *                          from which to send to processor i.
     *                          If known, this should be provided as an input.  Otherwise
     *                          it is an optional output that will return the starting location
     *                          (relative to the start of the array) for the received data from
     *                          processor i.
     * @param[in] known_recv    Are the received counts and displacements known.
     *                          If the received sizes are known, then they must be provided,
     *                          and an extra communication step is not necessary.  If the received
     *                          sizes are not know, then an extra communication step will occur
     *                          internally and the sizes and displacements will be returned (if
     * desired).
     */
    template<class type>
    int allToAll( const type *send_data,
                  const int send_cnt[],
                  const int send_disp[],
                  type *recv_data,
                  int *recv_cnt   = nullptr,
                  int *recv_disp  = nullptr,
                  bool known_recv = false ) const;

    /*!
     * \brief   Send a list of proccesor ids to communicate
     * \details This function communicates a list of proccesors to communicate.
     *    Given a list of ranks that we want to send/receieve data to/from, this routine
     *    will communicate that set to the other ranks returning the list of processors
     *    that want to communication with the current rank.
     *    Note: this routine will involved global communication
     * \param[in] ranks     List of ranks that the current rank wants to communicate with
     * \return              List of ranks that want to communicate with the current processor
     */
    std::vector<int> commRanks( const std::vector<int> &ranks ) const;


    /*!
     * \brief   Wait for a communication to finish
     * \details Wait for a communication to finish.
     *    Note: this does not require a communicator.
     * \param[in] request    Communication request to wait for (returned for Isend or Irecv)
     */
    static void wait( const Request &request );


    /*!
     * \brief   Wait for a communication to finish
     * \details Wait for a communication to finish.
     *    Note: this does not require a communicator.
     * \param[in] request    Communication request to wait for (returned for Isend or Irecv)
     */
    static void wait( Request2 request );


    /*!
     * \brief   Wait for any communication to finish.
     * \details This function waits for any of the given communication requests to finish.
     *    It returns the index of the communication request that finished.
     *    Note: this does not require a communicator.
     * \param[in] count      Number of communications to check
     * \param[in] request    Array of communication requests to wait for (returned for Isend or
     * Irecv)
     */
    static int waitAny( int count, const Request *request );


    /*!
     * \brief   Wait for any communication to finish.
     * \details This function waits for any of the given communication requests to finish.
     *    It returns the index of the communication request that finished.
     *    Note: this does not require a communicator.
     * \param[in] count      Number of communications to check
     * \param[in] request    Array of communication requests to wait for (returned for Isend or
     * Irecv)
     */
    static int waitAny( int count, Request2 *request );


    /*!
     * \brief   Wait for all communications to finish.
     * \details This function waits for all of the given communication requests to finish.
     *    Note: this does not require a communicator.
     * \param[in] count      Number of communications to check
     * \param[in] request    Array of communication requests to wait for (returned for Isend or
     * Irecv)
     */
    static void waitAll( int count, const Request *request );

    /*!
     * \brief   Wait for all communications to finish.
     * \details This function waits for all of the given communication requests to finish.
     *    Note: this does not require a communicator.
     * \param[in] count      Number of communications to check
     * \param[in] request    Array of communication requests to wait for (returned for Isend or
     * Irecv)
     */
    static void waitAll( int count, Request2 *request );


    /*!
     * \brief   Wait for some communications to finish.
     * \details This function waits for one (or more) communications to finish.
     *    It returns an array of the indicies that have finished.
     *    Note: this does not require a communicator.
     * \param[in] count      Number of communications to check
     * \param[in] request    Array of communication requests to wait for (returned for Isend or
     * Irecv)
     */
    static std::vector<int> waitSome( int count, const Request *request );


    /*!
     * \brief   Wait for some communications to finish.
     * \details This function waits for one (or more) communications to finish.
     *    It returns an array of the indicies that have finished.
     *    Note: this does not require a communicator.
     * \param[in] count      Number of communications to check
     * \param[in] request    Array of communication requests to wait for (returned for Isend or
     * Irecv)
     */
    static std::vector<int> waitSome( int count, Request2 *request );


    /*!
     * \brief   Nonblocking test for a message
     * \details This function performs a non-blocking test for a message.
     *    It will return the number of bytes in the message if a message with
     *    the specified source and tag (on the current communicator) is available.
     *    Otherwise it will return -1.
     * \param[in] source      source rank (-1: any source)
     * \param[in] tag         tag (-1: any tag)
     */
    int Iprobe( int source = -1, int tag = -1 ) const;


    /*!
     * \brief   Blocking test for a message
     * \details This function performs a blocking test for a message.
     *    It will return the number of bytes in the message when a message with
     *    the specified source and tag (on the current communicator) is available
     * \param[in] source      source rank (-1: any source)
     * \param[in] tag         tag (-1: any tag)
     */
    int probe( int source = -1, int tag = -1 ) const;


    /*!
     * \brief   Start a serial region
     * \details This function will serialize MPI processes so that they run
     *    one at a time.  A call to serializeStart must be followed by a call
     *    to serializeStop after the commands to be executed.
     *    Note: the ranks will be run in order.
     */
    void serializeStart();


    /*!
     * \brief   Stop a serial region
     * \details Stop a serial region.  See serializeStart for more information.
     */
    void serializeStop();


    /*!
     * \brief   Elapsed time
     * \details This function returns the elapsed time on the calling processor
     *    since an arbitrary point in the past (seconds).  It is a wrapper to MPI_Wtime.
     *    See "tick" for the timer resolution in seconds.
     *    The time may or may not be synchronized across processors depending on the MPI
     *    implementation.  Refer to MPI documentation for the desired platform for more information.
     */
    static double time();


    /*!
     * \brief   Timer resolution
     * \details This function returns the timer resolution used by "time"
     */
    static double tick();


    /*!
     * \brief   Change the level of the internal timers
     * \details This function changes the level of the timers used to profile MPI
     * \param[in] level         New level of the timers
     */
    static void changeProfileLevel( int level ) { profile_level = level; }


    //! Return the total number of MPI_Comm objects that have been created
    static inline size_t MPI_Comm_created() { return N_MPI_Comm_created; }

    //! Return the total number of MPI_Comm objects that have been destroyed
    static inline size_t MPI_Comm_destroyed() { return N_MPI_Comm_destroyed; }

    //! Return details about MPI
    static std::string info();

    //! Return the MPI version number { major, minor }
    static std::array<int, 2> version();

    //! Check if MPI is active
    static bool MPI_Active();

    //! Start MPI
    static void start_MPI( int argc_in, char *argv_in[], int profile_level = 0 );

    //! Stop MPI
    static void stop_MPI();


private: // data members
    // The internal MPI communicator
    Comm d_comm;

    // Is the communicator NULL
    bool d_isNull;

    // Do we want to manage this communicator
    bool d_manage;

    // Do we want to call MPI_abort instead of exit
    bool d_call_abort;

    // The rank and size of the communicatorcommunicator
    int d_rank, d_size;

    // Some attributes
    int d_maxTag;
    int *volatile d_currentTag;

    // The ranks of the comm in the global comm
    mutable int *volatile d_ranks;

    /* How many objects share the same underlying MPI communicator.
     * When the count goes to 0, the MPI comm will be free'd (assuming it was created
     * by an communicator).  This may not be perfect, but is likely to be good enough.
     * Note that for thread safety, any access to this variable should be blocked for thread safety.
     * The value of count MUST be volatile to ensure the correct value is always used.
     */
    std::atomic_int *volatile d_count;


private: // static data members
    // The level for the profiles of MPI
    static short profile_level;

    /* We want to keep track of how many MPI_Comm objects we have created over time.
     * Like the count, for thread safety this should be blocked, however the most likely error
     * caused by not blocking is a slight error in the MPI count.  Since this is just for reference
     * we do not need to block (recognizing that the value may not be 100% accurate).
     */
    static volatile uint32_t N_MPI_Comm_created;
    static volatile uint32_t N_MPI_Comm_destroyed;
};


//! Return the underlying MPI class for the object
template<class TYPE>
AMP_MPI getComm( const TYPE &obj );


} // namespace AMP


#endif
