#ifndef included_AMP_DeviceDataHelpers_hpp
#define included_AMP_DeviceDataHelpers_hpp

#include "AMP/matrices/data/DeviceDataHelpers.h"

#include <thrust/device_vector.h>
#include <thrust/execution_policy.h>
#include <thrust/extrema.h>
#include <thrust/for_each.h>
#include <thrust/inner_product.h>
#include <thrust/scan.h>

namespace AMP {
namespace LinearAlgebra {

template<typename TYPE>
void DeviceDataHelpers<TYPE>::fill_n( TYPE *x, const size_t N, const TYPE alpha )
{
    thrust::fill_n( thrust::device, x, N, alpha );
}

template<typename TYPE>
void DeviceDataHelpers<TYPE>::copy_n( TYPE *x, const size_t N, TYPE *y )
{
    thrust::copy_n( thrust::device, x, N, y );
}

template<typename TYPE>
void DeviceDataHelpers<TYPE>::exclusive_scan( TYPE *x, const size_t N, TYPE *y, TYPE init )
{
    thrust::exclusive_scan( thrust::device, x, x + N, y, init );
}

template<typename TYPE>
TYPE DeviceDataHelpers<TYPE>::max_element( TYPE *x, const size_t N )
{
    return *thrust::max_element( thrust::device, x, x + N );
}

template<typename TYPE>
TYPE DeviceDataHelpers<TYPE>::accumulate( TYPE *x, const size_t N, TYPE init )
{
    return thrust::reduce( thrust::device, x, x + N, init, thrust::plus<TYPE>() );
}

} // namespace LinearAlgebra
} // namespace AMP

#endif
