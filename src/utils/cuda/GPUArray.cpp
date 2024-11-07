#include "AMP/utils/Array.hpp"
#include "AMP/utils/cuda/CudaAllocator.h"
#include "AMP/utils/cuda/GPUFunctionTable.h"


namespace AMP {


/********************************************************
 *  Explicit instantiations of Array                     *
 ********************************************************/
template class Array<double, AMP::GPUFunctionTable, CudaDevAllocator<void>>;
template class Array<double, AMP::GPUFunctionTable, CudaManagedAllocator<void>>;
template class Array<float, AMP::GPUFunctionTable, CudaDevAllocator<void>>;
template class Array<float, AMP::GPUFunctionTable, CudaManagedAllocator<void>>;


} // namespace AMP
