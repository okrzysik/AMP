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
INSTANTIATE_SIMPLE_VECTOR( float,
                           VectorOperationsDevice<float>,
                           VectorDataDefault<float, AMP::ManagedAllocator<void>> );
INSTANTIATE_SIMPLE_VECTOR( double,
                           VectorOperationsDevice<double>,
                           VectorDataDefault<double, AMP::ManagedAllocator<void>> );
#endif

} // namespace AMP::LinearAlgebra
