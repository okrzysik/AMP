#ifndef included_AMP_Mutex
#define included_AMP_Mutex

#include <atomic>
#include <memory>
#include <mutex>


namespace AMP {


class AMP_MPI; // Forward decleration of AMP_MPI


/** \class Mutex
 * \brief Functions for locking/unlocking a mutex
 * \details This class provides basic routines for creating,
 *    locking, and unlocking a mutex <BR>
 *    The lock may be recursive, meaning that the same thread
 *    may lock and unlock the lock multiple times before releasing it.
 *    In this case unlock must be called the same number of times before
 *    another thread may lock the mutex.
 */
class Mutex
{
public:
    //! Empty constructor (equivalent to Mutex(false) )
    Mutex();
    /** Default constructor
     * \param recursive     If set to true a thread may repeated lock a mutex.
     *                      If set to false an attept to repeatedly lock will throw an error.*/
    explicit Mutex( bool recursive );
    //! Destructor
    ~Mutex() = default;
    //! Copy constructor
    Mutex( const Mutex & ) = delete;
    //! Assignment operator
    Mutex &operator=( const Mutex & ) = delete;
    //! Lock the mutex
    void lock();
    //! Unlock the mutex
    void unlock();
    //! Try to lock the mutex and return true if successful
    bool tryLock();
    //! Return true if we already own the lock
    bool ownLock() const;

private:
    bool d_recursive;     // Is the lock recursive (this attribute cannot be changed)
    volatile int d_count; // lock_count
    volatile int d_id;    // thread_id
    std::mutex d_mutex;   // internal mutex
};


/*!
 * \brief Function to sycronize locking a mutex across a MPI communicator
 * \details This routine will sycronize locking the given mutex across the given MPI communicator.
 * This routine will return only when all threads across the given MPI communicator
 *   have acquired the mutex.  It is assumed that multiple threads with different
 *   communicators will attempt to lock the same mutex using this function
 */
void lock_MPI_Mutex( AMP::Mutex &mutex, const AMP::AMP_MPI &comm );


} // namespace AMP

#endif
