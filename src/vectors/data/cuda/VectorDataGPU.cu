#include "vectors/data/cuda/VectorDataGPU.h"
#include "vectors/data/cuda/VectorDataGPU.hpp"


// Specializations
template<>
std::string AMP::LinearAlgebra::VectorDataGPU<double>::VectorDataName() const
{
    return "VectorDataCPU<double>";
}
template<>
std::string AMP::LinearAlgebra::VectorDataGPU<float>::VectorDataName() const
{
    return "VectorDataCPU<float>";
}


// Explicit instantiations
template class AMP::LinearAlgebra::VectorDataGPU<double>;
template class AMP::LinearAlgebra::VectorDataGPU<float>;
