#include "AMP/utils/kdtree2.h"
#include "AMP/utils/Utilities.hpp"
#include "AMP/utils/kdtree2.hpp"


/********************************************************
 *  Explicit instantiations                              *
 ********************************************************/
#define INSANTIATE( TYPE )                \
    template class AMP::kdtree2<1, TYPE>; \
    template class AMP::kdtree2<2, TYPE>; \
    template class AMP::kdtree2<3, TYPE>; \
    template class AMP::kdtree2<4, TYPE>; \
    template class AMP::kdtree2<5, TYPE>
INSANTIATE( bool );
INSANTIATE( int );
INSANTIATE( float );
INSANTIATE( double );
template void AMP::Utilities::quicksort<double, uint64_t>( size_t, double *, uint64_t * );
