#include "AMP/utils/kdtree2.h"
#include "AMP/utils/kdtree2.hpp"


namespace AMP {

/********************************************************
 *  Explicit instantiations of kdtree2                   *
 ********************************************************/
template class kdtree2<1, int>;
template class kdtree2<2, int>;
template class kdtree2<3, int>;

} // namespace AMP
