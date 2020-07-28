#include "AMP/vectors/data/cuda/VectorDataGPU.h"


namespace AMP {
namespace LinearAlgebra {


template<typename TYPE>
std::string VectorDataGPU<TYPE>::VectorDataName() const
{
    std::string type = typeid( TYPE ).name();
    return "VectorDataGPU<" + type + ">";
}
template std::string VectorDataGPU<double>::VectorDataName() const;
template std::string VectorDataGPU<float>::VectorDataName() const;


} // namespace LinearAlgebra
} // namespace AMP
