#include "AMP/vectors/data/cuda/VectorDataGPU.h"
#include "AMP/vectors/data/cuda/VectorDataGPU.hpp"


namespace AMP::LinearAlgebra {


// Explicit instantiations
template class AMP::LinearAlgebra::VectorDataGPU<double>;
template class AMP::LinearAlgebra::VectorDataGPU<float>;


} // namespace AMP::LinearAlgebra
