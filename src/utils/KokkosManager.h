#ifndef included_AMP_KokkosManager
#define included_AMP_KokkosManager

#include <string>
#include <vector>


namespace AMP::Utilities {


void initializeKokkos( std::vector<char *> args );
void finalizeKokkos();


} // namespace AMP::Utilities

#endif
