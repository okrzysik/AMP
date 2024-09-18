#include "AMP/matrices/data/DeviceDataHelpers.h"
#include "AMP/matrices/data/DeviceDataHelpers.hpp"

#include <cstddef>
#include <cstdint>
#include <cuda.h>

// explicit instantiations
template class AMP::LinearAlgebra::DeviceDataHelpers<double>;
template class AMP::LinearAlgebra::DeviceDataHelpers<float>;
template class AMP::LinearAlgebra::DeviceDataHelpers<size_t>;
template class AMP::LinearAlgebra::DeviceDataHelpers<int32_t>;
template class AMP::LinearAlgebra::DeviceDataHelpers<long long>;
