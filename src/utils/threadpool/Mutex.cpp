#include "AMP/utils/threadpool/Mutex.h"

#include <cmath>
#include <cstdlib>
#include <thread>


// Get the id for the current thread
static inline int getThreadId()
{
    auto id = std::this_thread::get_id();
    static std::hash<std::thread::id> hasher;
    int hash = static_cast<int>( hasher( id ) );
    hash     = std::abs( hash );
    return hash;
}


/******************************************************************
 * Mutex class                                                     *
 ******************************************************************/
Mutex::Mutex()
{
    d_data.reset( new data_struct );
    d_data->recursive = false;
    d_data->count     = 0;
    d_data->id        = -1;
}
Mutex::Mutex( bool recursive )
{
    d_data.reset( new data_struct );
    d_data->recursive = recursive;
    d_data->count     = 0;
    d_data->id        = -1;
}
Mutex::Mutex( const Mutex &rhs ) : d_data( rhs.d_data ) {}
Mutex &Mutex::operator=( const Mutex &rhs )
{
    if ( this == &rhs ) // protect against invalid self-assignment
        return *this;
    this->d_data = rhs.d_data;
    return *this;
}
Mutex::~Mutex() {}
void Mutex::lock() const
{
    // Check if we already own the lock
    int id = getThreadId();
    if ( d_data->count > 0 && d_data->id == id ) {
        if ( !d_data->recursive )
            throw std::logic_error( "Lock is already locked and non-recursive" );
        // Increment the lock count
        d_data->count++;
        return;
    }
    // Acquire the lock
    d_data->mutex.lock();
    if ( d_data->count != 0 ) // If we are getting the lock, the count must be 0
        throw std::logic_error( "Internal error" );
    d_data->count = 1; // Change lock count after acquiring mutex
    d_data->id    = id;
}
bool Mutex::tryLock() const
{
    // Check if we already own the lock
    int id = getThreadId();
    if ( d_data->count > 0 && d_data->id == id ) {
        if ( !d_data->recursive )
            return false;
        // Increment the lock count and return
        d_data->count++;
        return true;
    }
    // Try and acquire the lock
    bool success = d_data->mutex.try_lock();
    if ( success ) {
        if ( d_data->count != 0 ) // If we are getting the lock, the count must be 0
            throw std::logic_error( "Internal error" );
        d_data->count = 1; // Change lock count after acquiring mutex
        d_data->id    = id;
    }
    return success;
}
void Mutex::unlock() const
{
    // Check if we already own the lock
    int id = getThreadId();
    if ( d_data->count <= 0 )
        throw std::logic_error( "Trying to release a lock that has not been locked" );
    if ( d_data->id != id )
        throw std::logic_error( "Thread that does not own lock is attempting to release" );
    // Release the lock
    d_data->count--; // Change lock count before releasing mutex
    if ( d_data->count == 0 ) {
        d_data->id = -1;
        d_data->mutex.unlock();
    }
}
bool Mutex::ownLock() const
{
    int id = getThreadId();
    if ( d_data->count > 0 && d_data->id == id )
        return true;
    return false;
}
