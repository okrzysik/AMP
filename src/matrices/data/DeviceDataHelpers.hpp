#ifndef included_AMP_DeviceDataHelpers_hpp
#define included_AMP_DeviceDataHelpers_hpp

#include "AMP/matrices/data/DeviceDataHelpers.h"

#include <thrust/device_vector.h>
#include <thrust/execution_policy.h>
#include <thrust/for_each.h>
#include <thrust/inner_product.h>

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

} // namespace LinearAlgebra
} // namespace AMP

#endif
