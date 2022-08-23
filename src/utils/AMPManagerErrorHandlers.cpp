#include "AMP/AMP_TPLs.h"
#include "AMP/utils/AMP_MPI.I"
#include "AMP/utils/AMPManager.h"
#include "AMP/utils/AMP_MPI.I"
#include "AMP/utils/Utilities.h"

#include "StackTrace/ErrorHandlers.h"
#include "StackTrace/Utilities.h"

#include <cstring>
#include <sstream>


// Include external packages for startup/shutdown
// clang-format off
#undef NULL_USE
#ifdef USE_CUDA
    #include <cuda.h>
    #include <cuda_runtime_api.h>
#endif
#ifdef AMP_USE_PETSC
    #include "petsc.h"
    #include "petscerror.h"
#endif
#ifdef AMP_USE_TIMER
    #include "MemoryApp.h"
#endif
#ifdef AMP_USE_SAMRAI
    #include "SAMRAI/tbox/Logger.h"
    #include "SAMRAI/tbox/SAMRAIManager.h"
    #include "SAMRAI/tbox/StartupShutdownManager.h"
#endif
#ifdef AMP_USE_KOKKOS
    #include "AMP/utils/KokkosManager.h"
#endif
#ifdef AMP_USE_HDF5
    #include "hdf5.h"
#endif
// clang-format on


namespace AMP {


/****************************************************************************
 *  Function to terminate AMP with a message for exceptions                  *
 ****************************************************************************/
static int force_exit     = 0;
static bool printed_stack = false;
static void terminate_AMP2( StackTrace::abort_error &err )
{
    printed_stack = true;
    StackTrace::Utilities::terminate( err );
}
void AMPManager::terminate_AMP( std::string message )
{
    AMP_MPI comm( AMP_COMM_WORLD );
    if ( !printed_stack ) {
        // Print the call stack and memory usage
        std::stringstream msg;
        msg << message << std::endl;
        msg << "Bytes used = " << AMP::Utilities::getMemoryUsage() << std::endl;
        StackTrace::multi_stack_info stack;
        if ( abort_stackType == 1 ) {
            stack = StackTrace::getCallStack();
        } else if ( abort_stackType == 2 ) {
            stack = StackTrace::getAllCallStacks();
        } else if ( abort_stackType == 3 ) {
            stack = StackTrace::getGlobalCallStacks();
        }
        StackTrace::cleanupStackTrace( stack );
        auto data = stack.print();
        msg << std::endl;
        msg << "Stack Trace:\n";
        for ( const auto &i : data )
            msg << " " << i << std::endl;
        // Add a rank dependent wait to hopefully print the stack trace cleanly
        Utilities::sleep_ms( ( 100 * comm.getRank() ) / comm.getSize() );
        perr << msg.str();
        printed_stack = true;
        force_exit    = 1;
    }
    if ( force_exit > 1 ) {
        exit( -1 );
    } else if ( AMP::AMPManager::use_MPI_Abort == true ) {
        // Use MPI_abort (will terminate all processes)
        force_exit = 2;
        comm.abort();
    } else if ( force_exit > 0 ) {
        exit( -1 );
    } else {
        // Throw and standard exception (allows the use of try, catch)
        force_exit = 1;
        throw std::logic_error( message );
    }
}
void AMPManager::exitFun()
{
    if ( initialized != 1 || printed_stack )
        return;
    auto stack = StackTrace::getCallStack();
    for ( auto &elem : stack ) {
        if ( strcmp( elem.function.data(), "MPID_Abort" ) == 0 )
            return;
    }
    std::stringstream msg;
    msg << "Calling exit without calling shutdown\n";
    msg << "Bytes used = " << AMP::Utilities::getMemoryUsage() << std::endl;
    msg << "Stack Trace:\n";
    for ( auto &elem : stack )
        msg << "   " << elem.print() << std::endl;
    perr << msg.str();
}


/******************************************************************
 * Create custom error handler                                     *
 ******************************************************************/
#ifdef AMP_USE_HDF5
herr_t hdf5_error_handler( hid_t err_stack, void * )
{
    FILE *fid = tmpfile();
    H5Eprint2( err_stack, fid );
    H5Eclear2( err_stack );
    rewind( fid );
    char msg[1024];
    size_t N = fread( msg, 1, sizeof( msg ) - 1, fid );
    fclose( fid );
    msg[N]           = 0;
    std::string msg2 = "Error calling HDF5 routine:\n";
    AMP_ERROR( msg2 + msg );
    return 0;
}
#endif


/****************************************************************************
 *  Function to handle PETSc errors                                          *
 ****************************************************************************/
#ifdef AMP_USE_PETSC
static_assert( PETSC_VERSION_GE( 3, 7, 5 ), "AMP only supports PETSc 3.7.5 or greater" );
static PetscErrorCode petsc_err_handler( MPI_Comm,
                                         int line,
                                         const char *dir,
                                         const char *file,
                                         PetscErrorCode,
                                         PetscErrorType,
                                         const char *buf,
                                         void * )
{
    std::stringstream msg;
    msg << "PETSc error:" << std::endl;
    msg << "   File: " << dir << file << ", line: " << line << std::endl;
    msg << "   " << buf << std::endl;
    AMPManager::terminate_AMP( msg.str() );
    return 0;
}
#endif


/****************************************************************************
 *  Functions to handle MPI errors                                           *
 ****************************************************************************/
void AMPManager::setMPIErrorHandler()
{
#ifdef AMP_USE_MPI
    StackTrace::setMPIErrorHandler( getCommWorld().getCommunicator() );
    #ifdef AMP_USE_SAMRAI
    auto comm = SAMRAI::tbox::SAMRAI_MPI::getSAMRAIWorld().getCommunicator();
    StackTrace::setMPIErrorHandler( comm );
    #endif
#endif
}
void AMPManager::clearMPIErrorHandler()
{
#ifdef AMP_USE_MPI
    StackTrace::clearMPIErrorHandler( getCommWorld().getCommunicator() );
    #ifdef AMP_USE_SAMRAI
    auto comm = SAMRAI::tbox::SAMRAI_MPI::getSAMRAIWorld().getCommunicator();
    StackTrace::clearMPIErrorHandler( comm );
    #endif
#endif
}


/****************************************************************************
 *  Class to override the output appender for abort messages                 *
 ****************************************************************************/
#ifdef AMP_USE_SAMRAI
class SAMRAIAbortAppender : public SAMRAI::tbox::Logger::Appender
{
public:
    void logMessage( const std::string &msg, const std::string &file, const int line ) override
    {
        auto msg2 = "SAMRAIAbortAppender called from " + file + " at line " +
                    std::to_string( line ) + ":\n" + msg;
        StackTrace::Utilities::abort( msg2, SOURCE_LOCATION_CURRENT() );
    }
    SAMRAIAbortAppender()           = default;
    ~SAMRAIAbortAppender() override = default;
};
#endif


/****************************************************************************
 * Functions to set/clear the error handlers                                 *
 ****************************************************************************/
void AMPManager::setHandlers()
{
    // Set the MPI error handler for comm_world
    setMPIErrorHandler();
    // Set the error handlers for petsc
#ifdef AMP_USE_PETSC
    PetscPopSignalHandler();
    PetscPopErrorHandler();
    PetscPushErrorHandler( &petsc_err_handler, PETSC_NULL );
#endif
    // Set the error handlers for SAMRAI
#ifdef AMP_USE_SAMRAI
    SAMRAI::tbox::SAMRAI_MPI::setCallAbortInSerialInsteadOfExit( true );
    SAMRAI::tbox::SAMRAI_MPI::setCallAbortInParallelInsteadOfMPIAbort( true );
    auto appender = std::make_shared<SAMRAIAbortAppender>();
    SAMRAI::tbox::Logger::getInstance()->setAbortAppender( appender );
#endif
    // Set the error handlers for HDF5
#ifdef AMP_USE_HDF5
    hid_t error_stack = 0;
    H5E_auto2_t fun   = hdf5_error_handler;
    H5Eset_auto2( error_stack, fun, nullptr );
#endif
    // Set the terminate routine for runtime errors
    StackTrace::Utilities::setErrorHandlers( terminate_AMP2 );
    // Set atexit function
    std::atexit( exitFun );
    int err = std::at_quick_exit( exitFun );
    AMP_ASSERT( err == 0 );
}
void AMPManager::clearHandlers()
{
    // Don't call the global version of the call stack
    StackTrace::globalCallStackFinalize();
    // Clear the MPI error handler for comm_world
    clearMPIErrorHandler();
    // Clear error handlers for StackTrace
    StackTrace::Utilities::clearErrorHandlers();
    StackTrace::clearSignals();
    StackTrace::clearSymbols();
    // Clear the error handlers for petsc
#ifdef AMP_USE_PETSC
    PetscPopSignalHandler();
    PetscPopErrorHandler();
#endif
}


} // namespace AMP
