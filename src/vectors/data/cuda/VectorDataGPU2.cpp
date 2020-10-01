#include "AMP/vectors/data/cuda/VectorDataGPU.h"
#include "AMP/vectors/data/cuda/VectorDataGPU.hpp"


namespace AMP {
namespace LinearAlgebra {


// Explicit instantiations
template class AMP::LinearAlgebra::VectorDataGPU<double>;
template class AMP::LinearAlgebra::VectorDataGPU<float>;


} // namespace LinearAlgebra
} // namespace AMP
