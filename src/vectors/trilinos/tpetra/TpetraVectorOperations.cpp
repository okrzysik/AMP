#include "AMP/vectors/trilinos/tpetra/TpetraVectorOperations.hpp"

namespace AMP::LinearAlgebra {

template class TpetraVectorOperations<double, int32_t, long long>;
template class TpetraVectorOperations<float, int32_t, long long>;

} // namespace AMP::LinearAlgebra
