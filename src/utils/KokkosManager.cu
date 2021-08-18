#ifdef USE_KOKKOS
#ifdef USE_CUDA

#include "AMP/utils/KokkosManager.h"
#include <Kokkos_Core.hpp>

namespace AMP::Utilities{

void initializeKokkos( int argc, char **argv)
{
  Kokkos::initialize( argc, argv );
}

void finalizeKokkos()
{
  Kokkos::finalize();
}

}
#endif
#endif
