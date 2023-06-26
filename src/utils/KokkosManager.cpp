#include "AMP/utils/KokkosManager.h"
#include "AMP/AMP_TPLs.h"

#if defined( AMP_USE_KOKKOS ) || defined( AMP_USE_TRILINOS_KOKKOS )
    #define USE_KOKKOS
    #include <Kokkos_Core.hpp>
#endif


namespace AMP::Utilities {

#ifdef USE_KOKKOS
void initializeKokkos( int &argc, char **argv ) { Kokkos::initialize( argc, argv ); }
void finalizeKokkos() { Kokkos::finalize(); }
#else
void initializeKokkos( int &, char ** ) {}
void finalizeKokkos() {}
#endif

} // namespace AMP::Utilities
