#include "AMP/utils/kdtree2.h"
#include "AMP/utils/Utilities.hpp"
#include "AMP/utils/kdtree2.hpp"


/********************************************************
 *  Explicit instantiations                              *
 ********************************************************/
template class AMP::kdtree2<1, int>;
template class AMP::kdtree2<2, int>;
template class AMP::kdtree2<3, int>;
template class AMP::kdtree2<4, int>;
template class AMP::kdtree2<5, int>;
template void AMP::Utilities::quicksort<double, uint64_t>( size_t, double *, uint64_t * );
