#ifdef USE_DEVICE
    #include "AMP/vectors/operations/VectorOperationsDevice.h"
    #include "AMP/vectors/operations/VectorOperationsDevice.hpp"
#endif

// Explicit instantiations
#ifdef USE_DEVICE
template class AMP::LinearAlgebra::VectorOperationsDevice<double>;
template class AMP::LinearAlgebra::VectorOperationsDevice<float>;
#endif
