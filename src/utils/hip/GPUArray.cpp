#include "AMP/utils/Array.h"
#include "AMP/utils/Array.hpp"
#include "AMP/utils/FunctionTable.h"
#include "AMP/utils/FunctionTable.hpp"
#include "AMP/utils/hip/HipAllocator.h"
#include "AMP/utils/hip/GPUFunctionTable.h"
#include "AMP/utils/hip/GPUFunctionTable.hpp"


namespace AMP {


/********************************************************
 *  Explicit instantiations of Array                     *
 ********************************************************/
template class Array<double, AMP::GPUFunctionTable, HipDevAllocator<double>>;
template class Array<double, AMP::GPUFunctionTable, HipManagedAllocator<double>>;
template class Array<float, AMP::GPUFunctionTable, HipDevAllocator<float>>;
template class Array<float, AMP::GPUFunctionTable, HipManagedAllocator<float>>;


} // namespace AMP
