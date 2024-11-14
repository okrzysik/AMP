#include "AMP/vectors/VectorBuilder.h"
#include "AMP/discretization/MultiDOF_Manager.h"
#include "AMP/utils/memory.h"
#include "AMP/vectors/MultiVariable.h"
#include "AMP/vectors/MultiVector.h"
#include "AMP/vectors/VectorBuilder.hpp"

namespace AMP::LinearAlgebra {


/********************************************************
 *  Explicit instantiations                              *
 ********************************************************/
INSTANTIATE_ARRAY_VECTOR( float );
INSTANTIATE_ARRAY_VECTOR( double );
INSTANTIATE_SIMPLE_VECTOR( float, VectorOperationsDefault<float>, VectorDataDefault<double> );
INSTANTIATE_SIMPLE_VECTOR( double, VectorOperationsDefault<double>, VectorDataDefault<double> );

#ifdef USE_DEVICE
using float_op  = VectorOperationsDevice<float>;
using double_op = VectorOperationsDevice<double>;
using float_data  = VectorDataDefault<float, AMP::ManagedAllocator<void>>;
using double_data = VectorDataDefault<double, AMP::ManagedAllocator<void>>;
INSTANTIATE_SIMPLE_VECTOR( float, float_op, float_data );
INSTANTIATE_SIMPLE_VECTOR( double, double_op, double_data );
#endif

} // namespace AMP::LinearAlgebra
