#include "AMP/utils/Array.h"
#include "AMP/utils/FunctionTable.h"
#include "AMP/utils/cuda/GPUFunctionTable.h"
#include "AMP/utils/cuda/GPUUmemAllocator.h"

#include "AMP/utils/Array.hpp"
#include "AMP/utils/FunctionTable.hpp"
#include "AMP/utils/cuda/GPUFunctionTable.hpp"
#include "AMP/utils/cuda/GPUUmemAllocator.hpp"


namespace AMP {


/********************************************************
 *  Explicit instantiations of Array                     *
 ********************************************************/
template class Array<double, AMP::GPUFunctionTable, GPUUmemAllocator<double>>;


} // namespace AMP
