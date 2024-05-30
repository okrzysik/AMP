#include "AMP/AMP_TPLs.h"
#include "AMP/utils/AMPManager.h"

#ifdef AMP_USE_PETSC
    #include "petsc.h"
    #include "petscerror.h"
#endif

#include <chrono>


// Get the elapsed duration
[[maybe_unused]] static double
getDuration( const std::chrono::time_point<std::chrono::steady_clock> &start )
{
    auto stop  = std::chrono::steady_clock::now();
    int64_t ns = std::chrono::duration_cast<std::chrono::nanoseconds>( stop - start ).count();
    return 1e-9 * ns;
}


namespace AMP {


/****************************************************************************
 * Function to start/stop PETSc                                              *
 ****************************************************************************/
#ifdef AMP_USE_PETSC
static bool called_PetscInitialize = false;
double AMPManager::start_PETSc()
{
    auto start = std::chrono::steady_clock::now();
    if ( PetscInitializeCalled ) {
        called_PetscInitialize = false;
    } else {
        int nargsPetsc       = 1;
        const char *noMalloc = "-malloc no";
        char **petscArgs     = const_cast<char **>( &noMalloc );
        PetscInitialize( &nargsPetsc, &petscArgs, nullptr, nullptr );
        called_PetscInitialize = true;
    }
    #ifndef AMP_USE_MPI
    // Fix minor bug in petsc where first call to dup returns MPI_COMM_WORLD instead of a new comm
    AMP::AMP_MPI( MPI_COMM_WORLD ).dup();
    #endif
    return getDuration( start );
}
double AMPManager::stop_PETSc()
{
    double time = 0;
    if ( called_PetscInitialize ) {
        auto start = std::chrono::steady_clock::now();
        PetscPopSignalHandler();
        PetscPopErrorHandler();
        PetscFinalize();
        time = getDuration( start );
    }
    return time;
}
#else
double AMPManager::start_PETSc() { return 0; }
double AMPManager::stop_PETSc() { return 0; }
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
void AMPManager::set_PETSc_error_handler()
{
    PetscPopSignalHandler();
    PetscPopErrorHandler();
    PetscPushErrorHandler( &petsc_err_handler, nullptr );
}
void AMPManager::clear_PETSc_error_handler()
{
    PetscPopSignalHandler();
    PetscPopErrorHandler();
}
#else
void AMPManager::set_PETSc_error_handler() {}
void AMPManager::clear_PETSc_error_handler() {}
#endif

} // namespace AMP
