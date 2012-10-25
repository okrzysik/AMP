#ifndef included_ProfilerApp
#define included_ProfilerApp

#include "utils/ProfilerAppMacros.h"

#include <stdio.h>
#include <stdlib.h>
#include <iostream>


#if defined(WIN32) || defined(_WIN32) || defined(WIN64) || defined(_WIN64)
    #define USE_WINDOWS
#else
    #define USE_LINUX
#endif


#ifdef USE_WINDOWS
    // Windows
    #define _CRT_SECURE_NO_WARNINGS		// Supress depreciated warnings for visual studio
    #define NOMINMAX                    // Supress min max from being defined
    #include <windows.h>
    #include <string>
    #define TIME_TYPE LARGE_INTEGER
#elif defined(USE_LINUX)
    // Linux
    #include <sys/time.h>
    #include <pthread.h>
    #include <string.h>
    #define TIME_TYPE timeval
#else
    #error Unknown OS
#endif


namespace AMP {


#define BIT_WORD size_t                         // A unsigned integer data type (the larger the word size, the better the performance)
#define TRACE_SIZE 100                          // The maximum number of timers that will be checked for the trace logs
                                                // The actual number of timers is TRACE_SIZE * number of bits of BIT_WORD
                                                // Note: this only affects the trace logs, the number of timers is unlimited
#define MAX_TRACE_TRACE 65536                   // The maximum number of stored start and stop times per trace
                                                // Note: this is only used if store_trace is set, and should be a power of 2
                                                // Note: the maximum ammount of memory used per trace is 16*MAX_TRACE_TRACE bytes (plus the trace itself)
#define THREAD_HASH_SIZE 32                     // The size of the hash table to store the threads
#define TIMER_HASH_SIZE 128                     // The size of the hash table to store the timers


/** \class ProfilerApp
  *
  * This class provides some basic timing and profiling capabilities.
  * It works be providing start, stop, and save functions that allow the user to record 
  * the total time required for blocks of code.  The results are written to an ASCII data 
  * file that can be viewed directly or processed.  It is compatible with MPI and SAMRAI.
  * The overhead of the function is ~1us per call for start and stop.  the time for save 
  * depends of the amount of data to be written and the speed of the writes.  Additional 
  * details can be provided with set_store_trace which records the actual time (from startup)
  * for each start and stop call.  This method significantly increases the memory requirements.
  * Several preprocessor define statement are given for easy incorporation: \verbatim
  *    PROFILE_START(NAME) - Start profiling a block of code with the given name
  *                          The name must be unique to the file and there must be a corresponding stop statement.
  *    PROFILE_STOP(NAME)  - Stop profiling a block of code with the given name
  *                          The name must match the given start block, and there must only be one PROFILE_STOP for each 
  *                          PROFILE_START.  This records the current line number as the final line number for the block of code.
  *    PROFILE_STOP2(NAME) - Stop profiling a block of code with the given name
  *                          The name must match the given start block.  This is a special stop that does not use the current 
  *                          line as the final line of the block, but is provided for cases where there are multiple exit
  *                          paths for a given function or block of code.
  *    PROFILE_SAVE(FILE)  - Save the results of profiling to a file.   
  * \endverbatim
  * Note that these commands are global and will create a global profiler.  It is possible
  * for a user to create multiple profilers and this should not create any problems, but the 
  * class interface should be used. <BR>
  * All start/stop and enable support an optional argument level that specifies the level of detail 
  * for the timers.  All timers with a number greater than the current level in the profiler will be ignored.
  * The macros PROFILE_START and PROFILE_STOP automatically check the level for performance and calling an
  * unused timer adds ~10ns per call. <BR>
  * For repeated calls the timer adds ~ 5us per call with without trace info, and ~25us per call with full trace info. 
  * Most of this overhead is not in the time returned by the timer.
  * Note that when a timer is created the cost may be significantly higher, but this only occurs once per timer.  <BR>
  * Example usage:
  *    void my_function(void *arg) {
  *       PROFILE_START("my function");
  *       int k;
  *       for (int i=0; i<10; i++) {
  *          PROFILE_START("sub call");
  *          // Do some work
  *          if ( special_case1 ) {
  *             PROFILE_STOP2("sub call");
  *             break;
  *          }
  *          if ( special_case2 ) {
  *             PROFILE_STOP2("sub call");
  *             PROFILE_STOP2("my function");
  *             return;
  *          }
  *          // Some more work
  *          PROFILE_STOP2("sub call");
  *       }
  *       PROFILE_STOP("my function");
  *    }
  *    // Some where at the end of the calculation
  *    PROFILE_SAVE("filename");
  */
class ProfilerApp {
public:
    //! Constructor
    ProfilerApp( );

    //! Destructor
    ~ProfilerApp();

    /*!
     * \brief  Function to start profiling a block of code
     * \details  This function starts profiling a block of code until a corresponding stop is called.
     *   It is recommended to use PROFILE_START(message) to call this routine.  It will 
     *   automatically fill in the file name and the line number.  
     * @param message       Message to uniquely identify the block of code being profiled.
     *                      It must be a unique message to all start called within the same file.
     * @param filename      Name of the file containing the code
     * @param line          Line number containing the start command
     * @param level         Level of detail to include this timer (default is 0)
     *                      Only timers whos level is <= the level of the specified by enable will be included.
     */
    void start( const std::string& message, const char* filename, const int line, const int level=0 );

    /*!
     * \brief  Function to stop profiling a block of code
     * \details  This function stop profiling a block of code until a corresponding stop is called.
     *   It is recommended to use PROFILE_STOP(message) to call this routine.  It will 
     *   automatically fill in the file name and the line number.  
     * @param message       Message to uniquely identify the block of code being profiled.
     *                      It must match a start call.
     * @param filename      Name of the file containing the code
     * @param line          Line number containing the stop command
     * @param level         Level of detail to include this timer (default is 0)
     *                      Only timers whos level is <= the level of the specified by enable will be included.
     *                      Note: this must match the level in start
     */
    void stop( const std::string& message, const char* filename, const int line, const int level=0 );

    /*!
     * \brief  Function to save the profiling info
     * \details  This will save the current timer info
     * Note: .x.timer will automatically be appended to the filename, where x is the rank+1 of the process.
     * Note: .x.trace will automatically be appended to the filename when detailed traces are used.
     * @param filename      File name for saving the results
     */
    void save( const std::string& filename );

    /*!
     * \brief  Function to enable the timers
     * \details  This function will enable the current timer clase.  It supports an optional level argument
     * that specifies the level of detail to use for the timers. 
     * @param level         Level of detail to include this timer (default is 0)
     *                      Only timers whos level is <= the level of the specified by enable will be included.
     */
    void enable( int level=0 );

    //! Function to enable the timers (all current timers will be deleted)
    void disable( );

    /*!
     * \brief  Function to change if we are storing detailed trace information
     * \details  This function will change if we are storing detailed trace information (must be called before any start)
     *  Note: Enabling this option will store the starting and ending time for each call.
     *  This will allow the user to look at the detailed results to get trace information.
     *  However this will significantly increase the memory requirements for any traces
     *  that get called repeatedly and may negitivly impact the performance.
     * @param profile       Do we want to store detailed profiling data
     */
    void set_store_trace(bool profile=false);

    inline int get_level( ) const { return d_level; }

private:

    // Structure to store the info for a trace log
    struct store_trace {
        int N_calls;                // Number of calls to this block
        unsigned int id;            // This is a (hopefully) unique id that we can use for comparison
        BIT_WORD trace[TRACE_SIZE]; // Store the trace
        store_trace *next;          // Pointer to the next entry in the list
        double min_time;            // Store the minimum time spent in the given block (seconds)
        double max_time;            // Store the maximum time spent in the given block (seconds)
        double total_time;          // Store the total time spent in the given block (seconds)
        double *start_time;         // Store when start was called for the given trace (seconds from constructor call)
        double *end_time;           // Store when stop was called for the given trace (seconds from constructor call)
        // Constructor
        store_trace() {
            N_calls = 0;
            min_time = 0.0;
            max_time = 0.0;
            total_time = 0.0;
            next = NULL;
            start_time = NULL;
            end_time = NULL;
        }
        // Copy constuctor
        store_trace(const store_trace& rhs);
        // De-constructor
		~store_trace() {
            if ( start_time == NULL )
                delete [] start_time;
            if ( end_time == NULL )
                delete [] end_time;
            start_time = NULL;
            end_time = NULL;
		}
    };

    // Structure to store the timing information for a single block of code
    struct store_timer {
        bool is_active;             // Are we currently running a timer
        unsigned int id;            // A unique id for each timer
        unsigned int trace_index;   // The index of the current timer in the trace
        int N_calls;                // Number of calls to this block
        BIT_WORD trace[TRACE_SIZE]; // Store the current trace
        double min_time;            // Store the minimum time spent in the given block (seconds)
        double max_time;            // Store the maximum time spent in the given block (seconds)
        double total_time;          // Store the total time spent in the given block (seconds)
        store_trace *trace_head;    // Head of the trace-log list
        store_timer *next;          // Pointer to the next entry in the list
        TIME_TYPE start_time;       // Store when start was called for the given block
        // Constructor used to initialize key values
		store_timer() {
			is_active = false;
            id = 0;
            trace_index = 0;
            N_calls = 0;
            min_time = 0.0;
            max_time = 0.0;
            total_time = 0.0;
            trace_head = NULL;
            next = NULL;
            for (int i=0; i<TRACE_SIZE; i++)
                trace[i] = 0;
		}
    };
    
    // Structure to store the timing information for a single block of code
    struct store_timer_info {
        int start_line;             // The starting line for the timer
        int stop_line;              // The ending line for the timer
        size_t id;                  // A unique id for each timer
        std::string message;        // The message to identify the block of code
        std::string filename;       // The file containing the block of code to be timed
        volatile store_timer_info *next; // Pointer to the next entry in the list
        // Constructor used to initialize key values
		store_timer_info() {
            id = 0;
            message = "";
            filename = "";
            start_line = -1;
            stop_line = -1;
		}
    };
    
    // Structure to store thread specific information
    struct thread_info {
        size_t id;                          // The id of the calling thread
        int thread_num;                     // The internal id of the thread
        unsigned int N_timers;              // The number of timers seen by the current thread
        volatile thread_info *next;         // Pointer to the next entry in the head list
        BIT_WORD active[TRACE_SIZE];        // Store the current active traces
        store_timer *head[TIMER_HASH_SIZE]; // Store the timers in a hash table
        // Constructor used to initialize key values
		thread_info() {
            #ifdef USE_WINDOWS
                id = NULL;
            #elif defined(USE_LINUX)
                id = 0;
            #endif
            N_timers = 0;
            next = NULL;
            for (int i=0; i<TRACE_SIZE; i++)
                active[i] = 0;
            for (int i=0; i<TIMER_HASH_SIZE; i++)
                head[i] = NULL;
		}
    };
    
    // Store thread specific info (use a small hash table to make searching faster)
    volatile int N_threads;
    volatile thread_info *thread_head[THREAD_HASH_SIZE];

    // Store the global timer info in a hash table
    volatile int N_timers;
    volatile store_timer_info *timer_table[TIMER_HASH_SIZE];

    // Function to return a pointer to the thread info (or create it if necessary)
    thread_info* get_thread_data( );

    // Function to return the appropriate timer block
    store_timer* get_block( const char* message, const char* filename, const int start, const int stop );

    // Function to return a hopefully unique id based on the message and filename
    static size_t get_timer_id( const char* message, const char* filename );

    // Function to return a hopefully unique id based on the active bit array
    static inline unsigned int get_trace_id( size_t N, const BIT_WORD *trace );

    // Function to return the string of active timers
    std::string get_active_list( BIT_WORD *active, unsigned int myIndex, thread_info *head );

    // Function to get the hash index given a timer id
    unsigned int get_timer_hash( unsigned int id );

    // Function to get the hash index given a thread id
    unsigned int get_thread_hash( unsigned int id );

    // Handle to a mutex lock
    #ifdef USE_WINDOWS
        HANDLE lock;                // Handle to a mutex lock
    #elif defined(USE_LINUX)
        pthread_mutex_t lock;       // Handle to a mutex lock
    #endif
    
    // Misc variables
    bool store_trace_data;          // Do we want to store trace information
    char d_level;                   // Level of timing to use (default is 0, -1 is disabled)
    TIME_TYPE construct_time;       // Store when the constructor was called
    TIME_TYPE frequency;            // Clock frequency (only used for windows)
};


}


// The global profiler
extern AMP::ProfilerApp global_profiler;


#endif


