#include "AMP/utils/KokkosManager.h"
#include "AMP/AMP_TPLs.h"
#include "AMP/utils/Utilities.h"

#if defined( AMP_USE_KOKKOS ) || defined( AMP_USE_TRILINOS_KOKKOS )
    #define USE_KOKKOS
    #include <Kokkos_Core.hpp>
#endif


namespace AMP::Utilities {

#ifdef USE_KOKKOS
void initializeKokkos( std::vector<char *> args )
{
    // Set some basic environmental variables
    if ( getenv( "OMP_PROC_BIND" ).empty() )
        setenv( "OMP_PROC_BIND", "false" );
    // Check if we need to set the number of threads
    bool setThreads = !getenv( "OMP_NUM_THREADS" ).empty();
    for ( size_t i = 0; i < args.size(); i++ ) {
        setThreads = setThreads || strncmp( args[i], "--threads", 9 ) == 0;
        setThreads = setThreads || strncmp( args[i], "--kokkos-threads", 16 ) == 0;
        setThreads = setThreads || strncmp( args[i], "--kokkos-num-threads", 20 ) == 0;
    }
    char defaultThreads[] = "--kokkos-num-threads=3";
    if ( !setThreads )
        args.push_back( defaultThreads );
    // Initialize kokkos
    int argc = args.size();
    Kokkos::initialize( argc, args.data() );
}
void finalizeKokkos() { Kokkos::finalize(); }
#else
void initializeKokkos( const std::vector<std::string> & ) {}
void finalizeKokkos() {}
#endif

} // namespace AMP::Utilities
