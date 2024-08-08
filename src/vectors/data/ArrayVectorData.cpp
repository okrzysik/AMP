#include "AMP/vectors/data/ArrayVectorData.h"
#ifdef USE_DEVICE
    #include "AMP/utils/device/GPUFunctionTable.h"
#endif

#include "AMP/utils/memory.h"


// Explicit instantiations
template class AMP::LinearAlgebra::ArrayVectorData<double>;
template class AMP::LinearAlgebra::ArrayVectorData<float>;

template class AMP::LinearAlgebra::
    ArrayVectorData<double, AMP::GPUFunctionTable, AMP::DeviceAllocator<double>>;
template class AMP::LinearAlgebra::
    ArrayVectorData<float, AMP::GPUFunctionTable, AMP::DeviceAllocator<float>>;
template class AMP::LinearAlgebra::
    ArrayVectorData<double, AMP::GPUFunctionTable, AMP::ManagedAllocator<double>>;
template class AMP::LinearAlgebra::
    ArrayVectorData<float, AMP::GPUFunctionTable, AMP::ManagedAllocator<float>>;
