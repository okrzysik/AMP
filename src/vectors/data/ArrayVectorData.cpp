#include "AMP/vectors/data/ArrayVectorData.h"
#ifdef USE_HIP
    #include "AMP/utils/hip/HipAllocator.h"
    #include "AMP/utils/hip/GPUFunctionTable.h"
#endif
#ifdef USE_CUDA
    #include "AMP/utils/cuda/CudaAllocator.h"
    #include "AMP/utils/cuda/GPUFunctionTable.h"
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
#ifdef USE_HIP
template class AMP::LinearAlgebra::
    ArrayVectorData<double, AMP::FunctionTable, AMP::HipDevAllocator<double>>;
template class AMP::LinearAlgebra::
    ArrayVectorData<float, AMP::FunctionTable, AMP::HipDevAllocator<float>>;
template class AMP::LinearAlgebra::
    ArrayVectorData<double, AMP::FunctionTable, AMP::HipManagedAllocator<double>>;
template class AMP::LinearAlgebra::
    ArrayVectorData<float, AMP::FunctionTable, AMP::HipManagedAllocator<float>>;
#endif
