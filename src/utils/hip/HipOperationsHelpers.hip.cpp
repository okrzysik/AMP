#include "AMP/utils/device/DeviceOperationsHelpers.h"
#include "AMP/utils/hip/HipOperationsHelpers.hpp"


// Explicit instantiations
template class AMP::LinearAlgebra::DeviceOperationsHelpers<double>;
template class AMP::LinearAlgebra::DeviceOperationsHelpers<float>;
