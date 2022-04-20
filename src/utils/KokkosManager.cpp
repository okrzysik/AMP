#include "AMP/AMP_TPLs.h"

#ifdef AMP_USE_KOKKOS

    #include "AMP/utils/KokkosManager.h"
    #include <Kokkos_Core.hpp>

namespace AMP::Utilities {

void initializeKokkos( int &argc, char **argv ) { Kokkos::initialize( argc, argv ); }

void finalizeKokkos() { Kokkos::finalize(); }

} // namespace AMP::Utilities
#endif
