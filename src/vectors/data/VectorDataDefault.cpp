#include "AMP/vectors/data/VectorDataDefault.h"
#include "AMP/utils/memory.h"
#include "AMP/vectors/data/VectorDataDefault.hpp"


// Explicit instantiations
template class AMP::LinearAlgebra::VectorDataDefault<double>;
template class AMP::LinearAlgebra::VectorDataDefault<float>;

template class AMP::LinearAlgebra::VectorDataDefault<double, AMP::DeviceAllocator<double>>;
template class AMP::LinearAlgebra::VectorDataDefault<float, AMP::DeviceAllocator<float>>;
template class AMP::LinearAlgebra::VectorDataDefault<double, AMP::ManagedAllocator<double>>;
template class AMP::LinearAlgebra::VectorDataDefault<float, AMP::ManagedAllocator<float>>;
