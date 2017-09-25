#include "vectors/operations/cuda/VectorOperationsCuda.h"
#include "vectors/operations/cuda/VectorOperationsCuda.hpp"


// Explicit instantiations
template class AMP::LinearAlgebra::VectorOperationsCuda<double>;
template class AMP::LinearAlgebra::VectorOperationsCuda<float>;
