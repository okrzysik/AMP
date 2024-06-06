#include "AMP/vectors/data/ArrayVectorData.h"
#ifdef USE_CUDA
    #include "AMP/utils/cuda/CudaAllocator.h"
#endif

// Explicit instantiations
template class AMP::LinearAlgebra::ArrayVectorData<double>;
template class AMP::LinearAlgebra::ArrayVectorData<float>;

#ifdef USE_CUDA
template class AMP::LinearAlgebra::
    ArrayVectorData<double, AMP::GPUFunctionTable, AMP::CudaDevAllocator<double>>;
template class AMP::LinearAlgebra::
    ArrayVectorData<float, AMP::GPUFunctionTable, AMP::CudaDevAllocator<float>>;
template class AMP::LinearAlgebra::
    ArrayVectorData<double, AMP::GPUFunctionTable, AMP::CudaManagedAllocator<double>>;
template class AMP::LinearAlgebra::
    ArrayVectorData<float, AMP::GPUFunctionTable, AMP::CudaManagedAllocator<float>>;
#endif
