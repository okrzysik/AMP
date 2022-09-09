#include "AMP/vectors/data/VectorDataIterator.h"
#include "AMP/vectors/data/VectorDataIterator.hpp"


// Explicit instantiations
template class AMP::LinearAlgebra::VectorDataIterator<const double>;
template class AMP::LinearAlgebra::VectorDataIterator<const float>;
template class AMP::LinearAlgebra::VectorDataIterator<double>;
template class AMP::LinearAlgebra::VectorDataIterator<float>;
template class AMP::LinearAlgebra::VectorDataIterator<std::complex<const double>>;
template class AMP::LinearAlgebra::VectorDataIterator<std::complex<const float>>;
template class AMP::LinearAlgebra::VectorDataIterator<std::complex<double>>;
template class AMP::LinearAlgebra::VectorDataIterator<std::complex<float>>;
