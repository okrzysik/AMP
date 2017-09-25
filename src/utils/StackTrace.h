#ifndef included_AMP_StackTrace
#define included_AMP_StackTrace

#include "utils/PIO.h"

#include <functional>
#include <iostream>
#include <set>
#include <stdio.h>
#include <stdlib.h>
#include <thread>
#include <vector>


namespace AMP {
namespace StackTrace {


struct stack_info {
    void *address;
    void *address2;
    std::string object;
    std::string function;
    std::string filename;
    int line;
    //! Default constructor
    stack_info() : address( nullptr ), address2( nullptr ), line( 0 ) {}
    //! Print the stack info
    std::string print() const;
};


struct multi_stack_info {
    int N;
    stack_info stack;
    std::vector<multi_stack_info> children;
    //! Default constructor
    multi_stack_info() : N( 0 ) {}
    //! Print the stack info
    std::vector<std::string> print( const std::string &prefix = std::string() ) const;
};


//! Function to return the current call stack
std::vector<stack_info> getCallStack();

//! Function to return the current call stack for the given thread
std::vector<stack_info> getCallStack( std::thread::native_handle_type id );

//! Function to return the current call stack for all threads
std::vector<multi_stack_info> getAllCallStacks();


//! Function to return the current call stack for the current thread
std::vector<void *> backtrace();

//! Function to return the current call stack for the given thread
std::vector<void *> backtrace( std::thread::native_handle_type id );

//! Function to return the stack info for a given address
stack_info getStackInfo( void *address );


//! Function to return the stack info for a given address
std::vector<stack_info> getStackInfo( const std::vector<void *> &address );


//! Function to return the signal name
std::string signalName( int signal );


/*!
 * Return the symbols from the current executable (not availible for all platforms)
 * @return              Returns 0 if sucessful
 */
int getSymbols( std::vector<void *> &address,
                std::vector<char> &type,
                std::vector<std::string> &obj );


/*!
 * Return the name of the executable
 * @return      Returns the name of the executable (usually the full path)
 */
std::string getExecutable();


/*!
 * Return the search path for the symbols
 * @return              Returns the search path for the symbols
 */
std::string getSymPaths();


//!< Terminate type
enum class terminateType { signal, exception };

/*!
 * Set the error handlers
 * @param[in] abort     Function to terminate the program: abort(msg,type)
 */
void setErrorHandlers( std::function<void( std::string, terminateType )> abort );


/*!
 * Set the given signals to the handler
 * @param[in] signals   List of signals to set
 * @param[in] handler   Function to handle signals
 */
void setSignals( const std::vector<int> &signals, void ( *handler )( int ) );


//! Clear a signal set by setSignals
void clearSignal( int signal );


//! Clear all signals set by setSignals
void clearSignals();


//! Return a list of all signals that can be caught
std::vector<int> allSignalsToCatch();

//! Return a default list of signals to catch
std::vector<int> defaultSignalsToCatch();


//! Get a list of the active threads
std::vector<std::thread::native_handle_type> activeThreads();

//! Get a handle to this thread
std::thread::native_handle_type thisThread();


} // namespace StackTrace
} // namespace AMP

#endif
