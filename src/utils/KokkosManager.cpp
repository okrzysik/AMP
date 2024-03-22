#include "AMP/utils/KokkosManager.h"
#include "AMP/AMP_TPLs.h"
#include "AMP/utils/Utilities.h"

#if defined( AMP_USE_KOKKOS ) || defined( AMP_USE_TRILINOS_KOKKOS )
    #define USE_KOKKOS
    #include <Kokkos_Core.hpp>
#endif


namespace AMP::Utilities {


#if defined( AMP_USE_KOKKOS ) || defined( AMP_USE_TRILINOS_KOKKOS )
void initializeKokkos( int &argc_in, char *argv_in[] )
{
    // Copy the input arguments, swap the kokkos arguments to the end,
    //    and remove them from the input argument list
    int argc = argc_in;
    for ( int i = argc_in - 1; i >= 0; i-- ) {
        if ( strncmp( argv_in[i], "--kokkos-", 9 ) == 0 )
            std::swap( argv_in[i], argv_in[--argc_in] );
    }
    char *argv[1024] = { nullptr };
    for ( int i = 0; i < argc; i++ )
        argv[i] = argv_in[i];
    // Set some basic environmental variables
    if ( getenv( "OMP_PROC_BIND" ).empty() )
        setenv( "OMP_PROC_BIND", "false" );
    // Check if we need to set the number of threads
    bool setThreads = !getenv( "OMP_NUM_THREADS" ).empty();
    for ( int i = 0; i < argc; i++ ) {
        setThreads = setThreads || strncmp( argv[i], "--threads", 9 ) == 0;
        setThreads = setThreads || strncmp( argv[i], "--kokkos-threads", 16 ) == 0;
        setThreads = setThreads || strncmp( argv[i], "--kokkos-num-threads", 20 ) == 0;
    }
    char defaultThreads[] = "--kokkos-num-threads=3";
    if ( !setThreads )
        argv[argc++] = defaultThreads;
    // Initialize kokkos
    Kokkos::initialize( argc, argv );
}
void finalizeKokkos() { Kokkos::finalize(); }
#else
void initializeKokkos( int &, char *[] ) {}
void finalizeKokkos() {}
#endif

} // namespace AMP::Utilities
