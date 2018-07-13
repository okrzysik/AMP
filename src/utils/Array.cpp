#include "AMP/utils/Array.h"
#include "AMP/utils/Array.hpp"


namespace AMP {


/********************************************************
 *  ArraySize                                            *
 ********************************************************/
ArraySize::ArraySize( const std::vector<size_t> &N )
{
    d_ndim = N.size();
    d_N[0] = 0;
    d_N[1] = 1;
    d_N[2] = 1;
    d_N[3] = 1;
    d_N[4] = 1;
    for ( size_t i = 0; i < d_ndim; i++ )
        d_N[i] = N[i];
    d_length = 1;
    for ( unsigned long i : d_N )
        d_length *= i;
    if ( d_ndim == 0 )
        d_length = 0;
}


/********************************************************
 *  Explicit instantiations of Array                     *
 ********************************************************/
template class Array<uint8_t, FunctionTable>;
template class Array<uint16_t, FunctionTable>;
template class Array<uint32_t, FunctionTable>;
template class Array<uint64_t, FunctionTable>;
template class Array<int8_t, FunctionTable>;
template class Array<int16_t, FunctionTable>;
template class Array<int32_t, FunctionTable>;
template class Array<int64_t, FunctionTable>;
template class Array<float, FunctionTable>;
template class Array<double, FunctionTable>;


} // namespace AMP