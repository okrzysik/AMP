#ifndef included_AMP_AtomicList
#define included_AMP_AtomicList

#include <atomic>
#include <csignal>
#include <functional>


namespace AMP {


/** \class AtomicList
 *
 * \brief Maintain a sorted list of entries
 * \details This class implements a basic sorted list that is thread-safe and lock-free.
 *    Entries are stored smallest to largest according to the compare operator
 */
template<class TYPE, class COMPARE = std::less<TYPE>>
class AtomicList final
{
public:
    //! Default constructor
    AtomicList( size_t capacity           = 1024,
                const TYPE &default_value = TYPE(),
                const COMPARE &comp       = COMPARE() );

    //! Destructor
    ~AtomicList();

    /*!
     * \brief   Remove an item from the list
     * \details Find and remove first entry that meets the given criteria
     * @return          Return the item that matches the criteria,
     *                  or the default item if no item matches
     * @param compare   Comparison function object (i.e. an object that satisfies
     *                  the requirements of Compare) which returns true if the
     *                  given value meets the selection criteria.
     *                  The signature of the comparison function should be equivalent to:
     *                      bool cmp( const TYPE& value, ... );
     * @param args      Additional arguments for the comparison
     */
    template<typename Compare, class... Args>
    inline TYPE remove( Compare compare, const Args &... args );

    //! Remove the first from the list
    inline TYPE remove_first();

    /*!
     * \brief   Insert an item
     * \details Insert an item into the list
     * @param x         Item to insert
     */
    inline void insert( const TYPE &x );

    /*!
     * \brief   Return the size of the list
     * \details Return the number of items in the list
     */
    inline size_t size() const { return d_N.load(); }

    /*!
     * \brief   Check if the list is empty
     * \details Return true if the list is empty
     */
    inline bool empty() const { return d_N.load() == 0; }

    /*!
     * \brief   Clear the list
     * \details Removes all entries from the list
     */
    inline void clear()
    {
        while ( !empty() )
            remove_first();
    }

    /*!
     * \brief   Return the capacity of the list
     * \details Return the maximum number of items the list can hold
     */
    inline constexpr size_t capacity() const { return d_capacity; }

    /*!
     * \brief   Return the availible space
     * \details Return the current availible space (capacity - size)
     */
    inline constexpr size_t availible() const { return d_capacity - size(); }


    /*!
     * \brief   Check the list
     * \details Perform a series of checks to verify the list is in a stable state.
     *    Note: This function is only partially thread-safe: it will block all other
     *    operations on the list, but check may fail if we caught a thread modifying the list.
     *    It is intended for debugging purposes only!
     * @return          This function returns true if the list is in a good working state
     */
    inline bool check();


    //! Return the total number of inserts since object creation
    inline size_t N_insert() const { return d_N_insert.load(); }


    //! Return the total number of removals since object creation
    inline size_t N_remove() const { return d_N_remove.load(); }

private:
    // Data members
    COMPARE d_compare;
    const size_t d_capacity;
    volatile TYPE d_default;
    volatile TYPE *d_objects;
    volatile std::atomic_int32_t d_N;
    volatile std::atomic_int32_t *d_next;
    volatile std::atomic_int32_t d_unused;
    volatile std::atomic_int64_t d_N_insert;
    volatile std::atomic_int64_t d_N_remove;

private:
    inline int lock( int i )
    {
        if ( i == -1 )
            return -1;
        int tmp = 0;
        do {
            tmp = d_next[i].fetch_and( 0 );
        } while ( tmp == 0 );
        return tmp;
    }
    inline void unlock( int i, int value )
    {
        if ( i != -1 )
            d_next[i].fetch_or( value );
    }
    inline int get_unused()
    {
        int i = 0;
        do {
            i = d_unused.fetch_and( 0 );
        } while ( i == 0 );
        d_unused.fetch_or( -( d_next[i] + 4 ) + 1 );
        d_next[i] = -3;
        return i;
    }
    inline void put_unused( int i )
    {
        int j = 0;
        do {
            j = d_unused.exchange( j );
        } while ( j == 0 );
        d_next[i] = -3 - j;
        d_unused.fetch_or( i );
    }


public:
    AtomicList( const AtomicList & ) = delete;
    AtomicList &operator=( const AtomicList & ) = delete;
};


} // namespace AMP

#include "AMP/utils/threadpool/AtomicList.hpp"

#endif
