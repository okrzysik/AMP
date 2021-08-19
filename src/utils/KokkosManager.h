#ifndef included_KokkosManager_H_
#define included_KokkosManager_H_
namespace AMP::Utilities {
void initializeKokkos( int &argc, char **argv );
void finalizeKokkos();
} // namespace AMP::Utilities
#endif
