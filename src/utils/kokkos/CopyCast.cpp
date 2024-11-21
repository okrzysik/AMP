#include "CopyCast.h"
#include "AMP/utils/kokkos/CopyCast.hpp"

namespace AMP::Utilities {

template struct copyCast_<double, float, AMP::Utilities::MemoryType::host>;
template struct copyCast_<float, double, AMP::Utilities::MemoryType::host>;
template struct copyCast_<double, float, AMP::Utilities::MemoryType::unregistered>;
template struct copyCast_<float, double, AMP::Utilities::MemoryType::unregistered>;

} // namespace AMP::Utilities