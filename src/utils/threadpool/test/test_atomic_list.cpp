#include "AMP/utils/AMPManager.h"
#include "AMP/utils/UnitTest.h"
#include "AMP/utils/Utilities.h"
#include "AMP/utils/threadpool/AtomicList.h"

#include <algorithm>
#include <atomic>
#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <functional>
#include <iostream>
#include <set>
#include <string>
#include <thread>
#include <vector>

using namespace AMP;


// Testing class to wrap a std::multiset to compare performance
template<class TYPE>
class AtomicList2 final
{
public:
    AtomicList2() : d_lock( 1 ) {}
    template<typename Compare, class... Args>
    inline TYPE remove( Compare compare, const Args &... args )
    {
        lock();
        auto it = d_data.begin();
        while ( !compare( *it, args... ) ) {
            ++it;
        }
        auto ans = *it;
        d_data.erase( it );
        unlock();
        return ans;
    }
    inline TYPE remove_first()
    {
        lock();
        auto it  = d_data.begin();
        auto ans = *it;
        d_data.erase( it );
        unlock();
        return ans;
    }
    inline void insert( const TYPE &x )
    {
        lock();
        d_data.insert( x );
        unlock();
    }
    inline size_t size()
    {
        lock();
        auto ans = d_data.size();
        unlock();
        return ans;
    }
    inline bool empty()
    {
        lock();
        auto ans = d_data.empty();
        unlock();
        return ans;
    }
    inline void clear()
    {
        lock();
        d_data.clear();
        unlock();
    }
    inline bool check() { return true; }

private:
    std::multiset<TYPE> d_data;
    AtomicOperations::int32_atomic d_lock;
    inline void lock()
    {
        int tmp = 0;
        do {
            tmp = AtomicOperations::atomic_fetch_and_and( &d_lock, 0 );
        } while ( tmp == 0 );
    }
    inline void unlock() { AtomicOperations::atomic_fetch_and_or( &d_lock, 1 ); }
};


// Modify the list
template<class LIST>
static void modify_list( LIST &list, const std::vector<int> &rnd )
{
    for ( size_t i = 0; i < rnd.size(); i++ ) {
        int r   = rnd[i];
        auto v1 = list.remove_first();
        auto v2 = list.remove( []( int ) { return true; } );
        auto v3 = list.remove( [r]( int v ) { return v >= ( r / 8 ); } );
        auto v4 = list.remove( [r]( int v ) { return v >= ( r / 4 ); } );
        auto v5 = list.remove( [r]( int v ) { return v >= ( r / 2 ); } );
        if ( v1 != -1 )
            list.insert( v1 );
        if ( v2 != -1 )
            list.insert( v2 );
        if ( v3 != -1 )
            list.insert( v3 );
        if ( v4 != -1 )
            list.insert( v4 );
        if ( v5 != -1 )
            list.insert( v5 );
    }
}


// Check that the entries in the list are correct
template<class LIST>
static bool check_list( const std::vector<int> &x, LIST &list )
{
    // insert the items in the list
    if ( x.size() != list.size() )
        return false;
    if ( !list.check() )
        return false;
    // Check that they are in sorted order
    auto x2 = x;
    std::sort( x2.begin(), x2.end() );
    bool pass = true;
    for ( int i : x2 )
        pass = pass && i == list.remove_first();
    return pass;
}


// Test the cost to insert items in the list (ns)
template<class LIST>
static int timeInsert( const std::vector<int> &x, LIST &list )
{
    const int N_it = 20;
    std::chrono::duration<double> time;
    time = time.zero();
    for ( int it = 0; it < N_it; it++ ) {
        list.clear();
        auto start = std::chrono::high_resolution_clock::now();
        for ( int i : x )
            list.insert( i );
        auto stop = std::chrono::high_resolution_clock::now();
        time += ( stop - start );
    }
    return 1e9 * time.count() / ( N_it * x.size() );
}

// Test the cost to remove (first)
inline void remove_first( AtomicList<int> &list ) { list.remove_first(); }
inline void remove_first( AtomicList2<int> &list ) { list.remove_first(); }
inline void remove_first( std::multiset<int> &list ) { list.erase( list.begin() ); }
template<class LIST>
static int timeRemoveFirst( const std::vector<int> &x, LIST &list )
{
    const int N_it = 10;
    std::chrono::duration<double> time;
    time = time.zero();
    for ( int it = 0; it < N_it; it++ ) {
        list.clear();
        for ( int i : x )
            list.insert( i );
        auto start = std::chrono::high_resolution_clock::now();
        for ( size_t i = 0; i < x.size(); i++ )
            remove_first( list );
        auto stop = std::chrono::high_resolution_clock::now();
        time += ( stop - start );
        if ( !list.empty() )
            throw std::logic_error( "List not empty" );
    }
    return 1e9 * time.count() / ( N_it * x.size() );
}


// Test the cost to remove (in order)
inline void remove_true( AtomicList<int> &list )
{
    list.remove( []( int ) { return true; } );
}
inline void remove_true( AtomicList2<int> &list )
{
    list.remove( []( int ) { return true; } );
}
inline void remove_true( std::multiset<int> &list ) { list.erase( list.begin() ); }
template<class LIST>
static int timeRemoveOrdered( const std::vector<int> &x, LIST &list )
{
    const int N_it = 10;
    std::chrono::duration<double> time;
    time = time.zero();
    time = time.zero();
    for ( int it = 0; it < N_it; it++ ) {
        list.clear();
        for ( int i : x )
            list.insert( i );
        auto start = std::chrono::high_resolution_clock::now();
        for ( size_t i = 0; i < x.size(); i++ )
            remove_true( list );
        auto stop = std::chrono::high_resolution_clock::now();
        time += ( stop - start );
        if ( !list.empty() )
            throw std::logic_error( "List not empty" );
    }
    return 1e9 * time.count() / ( N_it * x.size() );
}


// Test the cost to remove (out of order)
inline void remove_v( AtomicList<int> &list, int value )
{
    list.remove( [value]( int v ) { return v == value; } );
}
inline void remove_v( AtomicList2<int> &list, int value )
{
    list.remove( [value]( int v ) { return v == value; } );
}
inline void remove_v( std::multiset<int> &list, int value ) { list.erase( list.find( value ) ); }
template<class LIST>
static int timeRemoveUnordered( const std::vector<int> &x, LIST &list )
{
    const int N_it = 10;
    std::chrono::duration<double> time;
    time = time.zero();
    for ( int it = 0; it < N_it; it++ ) {
        list.clear();
        for ( int i : x )
            list.insert( i );
        auto start = std::chrono::high_resolution_clock::now();
        for ( int tmp : x )
            remove_v( list, tmp );
        auto stop = std::chrono::high_resolution_clock::now();
        time += ( stop - start );
    }
    return 1e9 * time.count() / ( N_it * x.size() );
}


// Test the cost to run
template<class LIST>
static int runSerial( const std::vector<int> &x, LIST &list )
{
    const int N_it  = 3;
    const int count = 50000;
    std::chrono::duration<double> time;
    time = time.zero();
    std::vector<int> rnd( count );
    for ( int i = 0; i < count; i++ )
        rnd[i] = rand();
    for ( int it = 0; it < N_it; it++ ) {
        list.clear();
        for ( int i : x )
            list.insert( i );
        auto start = std::chrono::high_resolution_clock::now();
        modify_list( list, rnd );
        auto stop = std::chrono::high_resolution_clock::now();
        time += ( stop - start );
    }
    return 1e9 * time.count() / ( count * N_it );
}
template<class LIST>
static int runParallel( const std::vector<int> &x, LIST &list, int N_threads )
{
    const int N_it  = 1;
    const int count = 20000;
    std::vector<std::thread> threads( N_threads );
    std::chrono::duration<double> time;
    time = time.zero();
    std::vector<std::vector<int>> rnd( N_threads );
    for ( int j = 0; j < N_threads; j++ ) {
        rnd[j].resize( count );
        for ( int i = 0; i < count; i++ )
            rnd[j][i] = rand();
    }
    for ( int it = 0; it < N_it; it++ ) {
        list.clear();
        for ( int i : x )
            list.insert( i );
        auto start = std::chrono::high_resolution_clock::now();
        for ( int i = 0; i < N_threads; i++ )
            threads[i] = std::thread( modify_list<LIST>, std::ref( list ), rnd[i] );
        for ( int i = 0; i < N_threads; i++ )
            threads[i].join();
        auto stop = std::chrono::high_resolution_clock::now();
        time += ( stop - start );
    }
    return 1e9 * time.count() / ( N_threads * count * N_it );
}


/******************************************************************
 * The main program                                                *
 ******************************************************************/
int main( int argc, char *argv[] )
{
    AMP::AMPManager::startup( argc, argv );
    UnitTest ut;

    const int N_threads = 8;

    // Create the list
    AtomicList<int> list( 1024, -1 );
    AtomicList2<int> list2;
    std::multiset<int> multiset;
    if ( list.size() == 0 && list.empty() && list.check() )
        ut.passes( "Initialize" );
    else
        ut.failure( "Initialize" );

    // Initialize the list with some empty values
    for ( int i = 0; i < 80; i++ )
        list.insert( rand() );
    list.insert( 2 );
    list.insert( 1 );
    list.insert( rand() );

    // Try to pull off a couple of values
    int v1 = list.remove( []( int a ) { return a == 1; } ); // Find the entry with 1
    int v2 = list.remove( []( int ) { return true; } );     // Get the first entry
    int v3 = list.remove( []( int ) { return false; } );    // Fail to get an entry
    if ( v1 == 1 && v2 == 2 && v3 == -1 && list.size() == 81 && list.check() )
        ut.passes( "Basic sanity test" );
    else
        ut.failure( "Basic sanity test" );

    // Create a vector of unordered data
    std::vector<int> x( 5 * N_threads );
    for ( int &i : x )
        i = rand();

    // Run some more checks
    list.clear();
    for ( int i : x )
        list.insert( i );
    check_list( x, list );

    // Test the cost to insert
    int t1 = timeInsert( x, list );
    int t2 = timeInsert( x, list2 );
    int t3 = timeInsert( x, multiset );
    printf( "insert time/item (AtomicList)  = %i ns\n", t1 );
    printf( "insert time/item (AtomicList2) = %i ns\n", t2 );
    printf( "insert time/item (multiset)    = %i ns\n", t3 );

    // Test the cost to remove (first)
    t1 = timeRemoveFirst( x, list );
    t2 = timeRemoveFirst( x, list2 );
    t3 = timeRemoveFirst( x, multiset );
    printf( "remove (first) time/item (AtomicList)  = %i ns\n", t1 );
    printf( "remove (first) time/item (AtomicList2) = %i ns\n", t2 );
    printf( "remove (first) time/item (multiset)    = %i ns\n", t3 );

    // Test the cost to remove (in order)
    t1 = timeRemoveOrdered( x, list );
    t2 = timeRemoveOrdered( x, list2 );
    t3 = timeRemoveOrdered( x, multiset );
    printf( "remove (ordered) time/item (AtomicList)  = %i ns\n", t1 );
    printf( "remove (ordered) time/item (AtomicList2) = %i ns\n", t2 );
    printf( "remove (ordered) time/item (multiset)    = %i ns\n", t3 );

    // Test the cost to remove (out of order)
    t1 = timeRemoveUnordered( x, list );
    t2 = timeRemoveUnordered( x, list2 );
    t3 = timeRemoveUnordered( x, multiset );
    printf( "remove (unordered) time/item (AtomicList)  = %i ns\n", t1 );
    printf( "remove (unordered) time/item (AtomicList2) = %i ns\n", t2 );
    printf( "remove (unordered) time/item (multiset)    = %i ns\n", t3 );

    // Read/write to the list and check the results
    t1 = runSerial( x, list );
    t2 = runSerial( x, list2 );
    if ( check_list( x, list ) && check_list( x, list2 ) )
        ut.passes( "Serial get/insert" );
    else
        ut.failure( "Serial get/insert" );
    printf( "\n" );
    printf( "serial time/item (AtomicList)  = %i ns\n", t1 );
    printf( "serial time/item (AtomicList2) = %i ns\n", t2 );

    // Have multiple threads reading/writing to the list simultaneously
    t1 = runParallel( x, list, N_threads );
    t2 = runParallel( x, list2, N_threads );
    if ( check_list( x, list ) && check_list( x, list2 ) )
        ut.passes( "Parallel get/insert" );
    else
        ut.failure( "Parallel get/insert" );
    printf( "parallel time/item (AtomicList)  = %i ns\n", t1 );
    printf( "parallel time/item (AtomicList2) = %i ns\n", t2 );

    // Try to over-fill the list
    while ( !list.empty() )
        list.remove_first();
    for ( size_t i = 1; i <= list.capacity(); i++ )
        list.insert( i );
    try {
        list.insert( list.capacity() + 1 );
        ut.failure( "List overflow" );
    } catch ( const std::exception &e ) {
        ut.passes( "List overflow" );
    } catch ( ... ) {
        ut.failure( "List overflow (unknown exception)" );
    }

    // Finished
    ut.report();
    auto N_errors = static_cast<int>( ut.NumFailGlobal() );
    ut.reset();
    AMP::AMPManager::shutdown();
    return N_errors;
}
