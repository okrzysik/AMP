#ifndef included_AMP_MemoryPool
#define included_AMP_MemoryPool


namespace AMP {


/** \class MemoryPool
 *
 * \brief Pool allocator
 * \details This class implements a basic fast pool allocator that is thread-safe.
 */
template<class TYPE, class INT_TYPE = int>
class MemoryPool final
{
public:
    //! Default constructor
    explicit MemoryPool( size_t size );

    //! destructor
    ~MemoryPool();

    /*!
     * \brief   Allocate an object
     * \details Allocates a new object from the pool
     * @return          Return the new pointer, or nullptr if there is no more room in the pool
     */
    inline TYPE *allocate();

    /*!
     * \brief   Insert an item
     * \details Insert an item into the list
     * @param ptr       The pointer to free
     */
    inline void free( TYPE *ptr );

private:
    // Data members
    volatile TYPE *d_objects;
    volatile AtomicOperations::int32_atomic d_next;

private:
    MemoryPool( const MemoryPool & );
    MemoryPool &operator=( const MemoryPool & );
};


} // namespace AMP

#include "AMP/threadpool/memory_pool.hpp"

#endif
