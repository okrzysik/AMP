#ifndef included_AMP_KokkosManager
#define included_AMP_KokkosManager

#include <string>
#include <vector>

namespace AMP::Utilities {


void initializeKokkos( int &argc, char *argv[] );
void finalizeKokkos();


} // namespace AMP::Utilities

#endif
