#ifndef included_Mutex
#define included_Mutex

#include <atomic>
#include <memory>
#include <mutex>


namespace AMP {


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
    //! Empty constructor (equivilent to Mutex(false) )
    Mutex();
    /** Default constructor
     * \param recursive     If set to true a thread may repeated lock a mutex.
     *                      If set to false an attept to repeatedly lock will throw an error.*/
    explicit Mutex( bool recursive );
    //! Destructor
    ~Mutex();
    //! Copy constructor
    Mutex( const Mutex & );
    //! Assignment operator
    Mutex &operator=( const Mutex & );
    //! Lock the mutex
    void lock() const;
    //! Unlock the mutex
    void unlock() const;
    //! Try to lock the mutex and return true if successful
    bool tryLock() const;
    //! Return true if we already own the lock
    bool ownLock() const;
    //! Invalidate and clear the mutex (advanced interface, use with caution)
    void invalidate() { d_data.reset(); }

private:
    struct data_struct {
        bool recursive;     // Is the lock recursive (this attribute cannot be changed)
        volatile int count; // lock_count
        volatile int id;    // thread_id
        std::mutex mutex;   // internal mutex
    };
    std::shared_ptr<data_struct> d_data;
};


} // namespace AMP

#endif
