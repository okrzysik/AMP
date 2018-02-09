#include "AMP/vectors/data/VectorDataCPU.h"
#include "AMP/vectors/data/VectorDataCPU.hpp"


// Specializations
template<>
std::string AMP::LinearAlgebra::VectorDataCPU<double>::VectorDataName() const
{
    return "VectorDataCPU<double>";
}
template<>
std::string AMP::LinearAlgebra::VectorDataCPU<float>::VectorDataName() const
{
    return "VectorDataCPU<float>";
}


// Explicit instantiations
template class AMP::LinearAlgebra::VectorDataCPU<double>;
template class AMP::LinearAlgebra::VectorDataCPU<float>;
