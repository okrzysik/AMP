#include "AMP/utils/hip/HipOperationsHelpers.hpp"
#include "AMP/utils/device/DeviceOperationsHelpers.h"


// Explicit instantiations
template class AMP::LinearAlgebra::DeviceOperationsHelpers<double>;
template class AMP::LinearAlgebra::DeviceOperationsHelpers<float>;
