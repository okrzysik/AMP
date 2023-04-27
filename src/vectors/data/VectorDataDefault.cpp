#include "AMP/vectors/data/VectorDataDefault.h"
#include "AMP/vectors/data/VectorDataDefault.hpp"


#ifdef USE_CUDA
#include "AMP/utils/cuda/CudaAllocator.h"
#endif


// Explicit instantiations
template class AMP::LinearAlgebra::VectorDataDefault<double>;
template class AMP::LinearAlgebra::VectorDataDefault<float>;
#ifdef USE_CUDA
template class AMP::LinearAlgebra::VectorDataDefault<double,AMP::CudaDevAllocator<double>>;
template class AMP::LinearAlgebra::VectorDataDefault<float,AMP::CudaDevAllocator<float>>;
template class AMP::LinearAlgebra::VectorDataDefault<double,AMP::CudaManagedAllocator<double>>;
template class AMP::LinearAlgebra::VectorDataDefault<float,AMP::CudaManagedAllocator<float>>;
#endif


