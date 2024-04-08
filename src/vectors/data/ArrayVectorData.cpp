#include "AMP/vectors/data/ArrayVectorData.h"


// Explicit instantiations
template class AMP::LinearAlgebra::ArrayVectorData<double>;
template class AMP::LinearAlgebra::ArrayVectorData<float>;

#ifdef USE_CUDA
template class AMP::LinearAlgebra::
    ArrayVectorData<double, AMP::FunctionTable, AMP::CudaDevAllocator<double>>;
template class AMP::LinearAlgebra::
    ArrayVectorData<float, AMP::FunctionTable, AMP::CudaDevAllocator<float>>;
template class AMP::LinearAlgebra::
    ArrayVectorData<double, AMP::FunctionTable, AMP::CudaManagedAllocator<double>>;
template class AMP::LinearAlgebra::
    ArrayVectorData<float, AMP::FunctionTable, AMP::CudaManagedAllocator<float>>;
#endif
