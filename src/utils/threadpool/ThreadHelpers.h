#ifndef included_AMP_ThreadPool_ThreadHelpers
#define included_AMP_ThreadPool_ThreadHelpers


#include <thread>
#include <vector>


namespace AMP::Thread {


//! Function to return the number of processors available
int getNumberOfProcessors();


//! Function to return the processor number that the current thread is running on
int getCurrentProcessor();


//! Return the native handle to the current thread
std::thread::native_handle_type getCurrentThread();


//! Function to return the affinity of the current process
std::vector<int> getProcessAffinity();


//! Function to set the affinity of the current process
void setProcessAffinity( const std::vector<int> &procs );


/*!
 *  Function to return the affinity of the given thread
 *  @param thread   The native thread handle
 */
std::vector<int> getThreadAffinity( std::thread::native_handle_type );


/*!
 *  Set the given thread to have the specified affinity
 *  @param thread   The native thread handle
 *  @param procs    The processors to use
 */
void setThreadAffinity( std::thread::native_handle_type, const std::vector<int> &procs );


} // namespace AMP::Thread

#endif
