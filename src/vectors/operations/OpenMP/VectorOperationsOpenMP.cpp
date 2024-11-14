#include "AMP/vectors/operations/OpenMP/VectorOperationsOpenMP.h"
#include "AMP/vectors/VectorBuilder.hpp"
#include "AMP/vectors/operations/OpenMP/VectorOperationsOpenMP.hpp"


namespace AMP::LinearAlgebra {


/********************************************************
 *  Explicit instantiations                              *
 ********************************************************/
template class VectorOperationsOpenMP<double>;
template class VectorOperationsOpenMP<float>;
INSTANTIATE_SIMPLE_VECTOR( float, VectorOperationsOpenMP<float>, VectorDataDefault<double> );
INSTANTIATE_SIMPLE_VECTOR( double, VectorOperationsOpenMP<double>, VectorDataDefault<double> );


} // namespace AMP::LinearAlgebra
