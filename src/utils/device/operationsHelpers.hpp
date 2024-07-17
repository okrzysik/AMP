#ifndef included_AMP_OperationsHelpers_hpp
#define included_AMP_OperationsHelpers_hpp

#ifdef USE_CUDA
#include <cuda.h>
#endif

#include <thrust/device_vector.h>
#include <thrust/execution_policy.h>
#include <thrust/for_each.h>
#include <thrust/inner_product.h>

namespace AMP {
  namespace LinearAlgebra {
    template<typename TYPE>
    void OperationsHelpers<TYPE>::copy_n(const TYPE *x, const size_t N, TYPE *y){
      thrust::copy_n( thrust::device, x, N, y);
    }
    template<typename TYPE>
    void OperationsHelpers<TYPE>::fill_n( TYPE *x, const size_t N, const TYPE data){
      thrust::fill_n( thrust::device, x, N, data);
    }
  }
}

#endif
