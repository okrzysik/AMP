#include "AMP/utils/Array.hpp"
#include "AMP/utils/hip/GPUFunctionTable.h"
#include "AMP/utils/hip/HipAllocator.h"

namespace AMP {

/********************************************************
 *  Explicit instantiations of Array                     *
 ********************************************************/
template class Array<double, AMP::GPUFunctionTable, HipDevAllocator<void>>;
template class Array<double, AMP::GPUFunctionTable, HipManagedAllocator<void>>;
template class Array<float, AMP::GPUFunctionTable, HipDevAllocator<void>>;
template class Array<float, AMP::GPUFunctionTable, HipManagedAllocator<void>>;

} // namespace AMP
