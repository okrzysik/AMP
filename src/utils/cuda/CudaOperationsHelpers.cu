#include "AMP/utils/device/DeviceOperationsHelpers.h"
#include "AMP/utils/cuda/CudaOperationsHelpers.hpp"


// Explicit instantiations
template class AMP::LinearAlgebra::DeviceOperationsHelpers<double>;
template class AMP::LinearAlgebra::DeviceOperationsHelpers<float>;
