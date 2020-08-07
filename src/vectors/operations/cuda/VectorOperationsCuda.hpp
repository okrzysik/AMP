#ifndef included_AMP_VectorOperationsCuda_hpp
#define included_AMP_VectorOperationsCuda_hpp

#include "AMP/vectors/Vector.h"
#include "AMP/vectors/data/VectorData.h"
#include "AMP/vectors/operations/VectorOperationsDefault.hpp"
#include "AMP/vectors/operations/cuda/VectorOperationsCuda.h"

#include <thrust/device_vector.h>
#include <thrust/execution_policy.h>
#include <thrust/for_each.h>
#include <thrust/inner_product.h>


namespace AMP {
namespace LinearAlgebra {


extern template class VectorOperationsCuda<double>; // Suppresses implicit instantiation below --
extern template class VectorOperationsCuda<float>;  // Suppresses implicit instantiation below --


/****************************************************************
 * Constructors                                                  *
 ****************************************************************/
template<typename TYPE>
std::shared_ptr<VectorOperations> VectorOperationsCuda<TYPE>::cloneOperations() const
{
    auto ptr = std::make_shared<VectorOperationsCuda<TYPE>>();
    return ptr;
}


/****************************************************************
 * Check that all data can be passed to cuda                     *
 ****************************************************************/
template<typename TYPE>
bool VectorOperationsCuda<TYPE>::checkData( ) const
{
    return d_VectorData->numberOfDataBlocks() == 1;
}
template<typename TYPE>
bool VectorOperationsCuda<TYPE>::checkData( const VectorOperations &x ) const
{
    return x.getVectorData()->numberOfDataBlocks() == 1 && d_VectorData->numberOfDataBlocks() == 1;
}
template<typename TYPE>
bool VectorOperationsCuda<TYPE>::checkData( const VectorOperations &x,
                                            const VectorOperations &y ) const
{
    return x.getVectorData()->numberOfDataBlocks() == 1 &&
           y.getVectorData()->numberOfDataBlocks() == 1 && d_VectorData->numberOfDataBlocks() == 1;
}
template<typename TYPE>
bool VectorOperationsCuda<TYPE>::checkData( const VectorData &x ) const
{
    return x->numberOfDataBlocks() == 1 ;
}
template<typename TYPE>
bool VectorOperationsCuda<TYPE>::checkData( const VectorData &x,
                                            const VectorData &y ) const
{
    return x.numberOfDataBlocks() == 1 &&
      y.numberOfDataBlocks() == 1 ;
}
template<typename TYPE>
bool VectorOperationsCuda<TYPE>::checkData( const VectorData &x,
                                            const VectorData &y,
					    const VectorData &z ) const
{
    return x.numberOfDataBlocks() == 1
        && y.numberOfDataBlocks() == 1 ;
        && z.numberOfDataBlocks() == 1 ;
}

//**********************************************************************
// Static functions that operate on VectorData objects

template<typename TYPE>
void VectorOperationsCuda<TYPE>::zero( VectorData &x )
{
  VectorOperationsCuda<TYPE>::setToScalar( 0.0, x );
}

template<typename TYPE>
void VectorOperationsCuda<TYPE>::setToScalar( double alpha, VectorData &x )
{
    bool useGPU = checkData();
    if ( useGPU ) {
        TYPE *data = x.getRawDataBlock<TYPE>( 0 );
        size_t N   = x.sizeOfDataBlock( 0 );
        thrust::fill_n( thrust::device, data, N, static_cast<TYPE>( alpha ) );
    } else {
        // Default to cpu version
        auto curMe = x.begin<TYPE>();
        auto last  = x.end<TYPE>();
        while ( curMe != last ) {
            *curMe = alpha;
            ++curMe;
        }
    }
    if ( hasGhosts() ) {
        auto &ghosts = getGhosts();
        for ( size_t i = 0; i != ghosts.size(); i++ )
            ghosts[i] = alpha;
    }
    // Override the status state since we set the ghost values
    *( x.getUpdateStatusPtr() ) = VectorData::UpdateState::UNCHANGED;
    // Wait for cuda data to complete
    if ( useGPU )
        cudaDeviceSynchronize();
}

template<typename TYPE>
void VectorOperationsCuda<TYPE>::setRandomValues( VectorData &x )
{
  // Default to VectorOperationsDefault (on cpu)
  return VectorOperationsDefault<TYPE>::setRandomValues(x);
}

template<typename TYPE>
void VectorOperationsCuda<TYPE>::setRandomValues( RNG::shared_ptr rng, VectorData &x )
{
    // Default to VectorOperationsDefault (on cpu)
  return VectorOperationsDefault<TYPE>::setRandomValues( rng, x );
}

template<typename TYPE>
void VectorOperationsCuda<TYPE>::copy( const VectorData &x, VectorData &y )
{
  if ( checkData( x, y ) ) {
        TYPE *data        = y.getRawDataBlock<TYPE>( 0 );
        const TYPE *xdata = x.getRawDataBlock<TYPE>( 0 );
        size_t N          = y.sizeOfDataBlock( 0 );
        thrust::copy_n( thrust::device, xdata, N, data );
        cudaDeviceSynchronize();
    } else {
        // Default to VectorOperationsDefault (on cpu)
      return VectorOperationsDefault<TYPE>::copy( x, y );
    }
}

template<typename TYPE>
void VectorOperationsCuda<TYPE>::scale( double alpha_in, VectorData &x )
{
    if ( checkData(x) ) {
        TYPE *data  = x.getRawDataBlock<TYPE>( 0 );
        size_t N    = x.sizeOfDataBlock( 0 );
        TYPE alpha  = static_cast<TYPE>( alpha_in );
        auto lambda = [alpha] __device__( TYPE y ) { return y * alpha; };
        thrust::transform( thrust::device, data, data + N, data, lambda );
        cudaDeviceSynchronize();
    } else {
        // Default to VectorOperationsDefault (on cpu)
      VectorOperationsDefault<TYPE>::scale( alpha_in, x );
    }
}

template<typename TYPE>
void VectorOperationsCuda<TYPE>::scale( double alpha_in, const VectorData &x, VectorData &y )
{
  if ( checkData( x, y ) ) {
        TYPE *data        = y.getRawDataBlock<TYPE>( 0 );
        const TYPE *xdata = x.getRawDataBlock<TYPE>( 0 );
        size_t N          = y.sizeOfDataBlock( 0 );
        TYPE alpha        = static_cast<TYPE>( alpha_in );
        auto lambda       = [alpha] __device__( TYPE x ) { return x * alpha; };
        thrust::transform( thrust::device, xdata, xdata + N, data, lambda );
        cudaDeviceSynchronize();
    } else {
        // Default to VectorOperationsDefault (on cpu)
      VectorOperationsDefault<TYPE>::scale( alpha_in, x, y );
    }
}

template<typename TYPE>
void VectorOperationsCuda<TYPE>::add( const VectorData &x, const VectorData &y, VectorData &z )
{
    if ( checkData( x, y ) ) {
        TYPE *data        = z.getRawDataBlock<TYPE>( 0 );
        const TYPE *xdata = x.getRawDataBlock<TYPE>( 0 );
        const TYPE *ydata = y.getRawDataBlock<TYPE>( 0 );
        size_t N          = z.sizeOfDataBlock( 0 );
        thrust::transform( thrust::device, xdata, xdata + N, ydata, data, thrust::plus<TYPE>() );
        cudaDeviceSynchronize();
    } else {
        // Default to VectorOperationsDefault (on cpu)
      VectorOperationsDefault<TYPE>::add( x, y, z );
    }
}

template<typename TYPE>
void VectorOperationsCuda<TYPE>::subtract( const VectorData &x, const VectorData &y, VectorData &z  )
{
  if ( checkData( x, y, z ) ) {
        TYPE *data        = z.getRawDataBlock<TYPE>( 0 );
        const TYPE *xdata = x.getRawDataBlock<TYPE>( 0 );
        const TYPE *ydata = y.getRawDataBlock<TYPE>( 0 );
        size_t N          = z.sizeOfDataBlock( 0 );
        thrust::transform( thrust::device, xdata, xdata + N, ydata, data, thrust::minus<TYPE>() );
        cudaDeviceSynchronize();
    } else {
        // Default to VectorOperationsDefault (on cpu)
      VectorOperationsDefault<TYPE>::subtract( x, y, z );
    }
}

template<typename TYPE>
void VectorOperationsCuda<TYPE>::multiply( const VectorData &x, const VectorData &y, VectorData &z )
{
  if ( checkData( x, y, z ) ) {
        TYPE *data        = z.getRawDataBlock<TYPE>( 0 );
        const TYPE *xdata = x.getRawDataBlock<TYPE>( 0 );
        const TYPE *ydata = y.getRawDataBlock<TYPE>( 0 );
        size_t N          = z.sizeOfDataBlock( 0 );
        thrust::transform(
            thrust::device, xdata, xdata + N, ydata, data, thrust::multiplies<TYPE>() );
        cudaDeviceSynchronize();
    } else {
        // Default to VectorOperationsDefault (on cpu)
      VectorOperationsDefault<TYPE>::multiply( x, y, z );
    }
}

template<typename TYPE>
void VectorOperationsCuda<TYPE>::divide( const VectorData &x, const VectorData &y, VectorData &z )
{
  if ( checkData( x, y, z ) ) {
        TYPE *data        = z.getRawDataBlock<TYPE>( 0 );
        const TYPE *xdata = x.getRawDataBlock<TYPE>( 0 );
        const TYPE *ydata = y.getRawDataBlock<TYPE>( 0 );
        size_t N          = z.sizeOfDataBlock( 0 );
        thrust::transform( thrust::device, xdata, xdata + N, ydata, data, thrust::divides<TYPE>() );
        cudaDeviceSynchronize();
    } else {
        // Default to VectorOperationsDefault (on cpu)
      VectorOperationsDefault<TYPE>::divide( x, y, z );
    }
}


template<typename TYPE>
void VectorOperationsCuda<TYPE>::reciprocal( const VectorData &x, VectorData &y )
{
  if ( checkData( x, y ) ) {
        TYPE *data        = y.getRawDataBlock<TYPE>( 0 );
        const TYPE *xdata = x.getRawDataBlock<TYPE>( 0 );
        size_t N          = y.sizeOfDataBlock( 0 );
        auto lambda       = [] __device__( TYPE x ) { return (TYPE) 1 / x; };
        thrust::transform( thrust::device, xdata, xdata + N, data, lambda );
        cudaDeviceSynchronize();
    } else {
        // Default to VectorOperationsDefault (on cpu)
      VectorOperationsDefault<TYPE>::reciprocal( x, y );
    }
}


template<typename TYPE>
void VectorOperationsCuda<TYPE>::linearSum( double alpha_in,
                                   const VectorData &x,
                                   double beta_in,
                                   const VectorData &y,
				   VectorData &z)
{
  if ( checkData( x, y, z ) ) {
        TYPE alpha        = alpha_in;
        TYPE beta         = beta_in;
        TYPE *data        = z.getRawDataBlock<TYPE>( 0 );
        const TYPE *xdata = x.getRawDataBlock<TYPE>( 0 );
        const TYPE *ydata = y.getRawDataBlock<TYPE>( 0 );
        size_t N          = z.sizeOfDataBlock( 0 );
        auto lambda = [alpha, beta] __device__( TYPE x, TYPE y ) { return alpha * x + beta * y; };
        thrust::transform( thrust::device, xdata, xdata + N, ydata, data, lambda );
        cudaDeviceSynchronize();
    } else {
        // Default to VectorOperationsDefault (on cpu)
      VectorOperationsDefault<TYPE>::linearSum( alpha_in, x, beta_in, y, z );
    }
}

template<typename TYPE>
void VectorOperationsCuda<TYPE>::axpy( double alpha_in, const VectorData &x, const VectorData &y, VectorData &z )
{
  VectorOperationsCuda<TYPE>::linearSum( alpha_in, x, 1.0, y, z );
}

template<typename TYPE>
void VectorOperationsCuda<TYPE>::axpby( double alpha_in, double beta_in, const VectorData &x, VectorData &z )
{
  VectorOperationsCuda<TYPE>::linearSum( alpha_in, x, beta_in, z, z );
}

template<typename TYPE>
void VectorOperationsCuda<TYPE>::abs( const VectorData &x, VectorData &y )
{
  if ( checkData( x, y ) ) {
        TYPE *data        = y.getRawDataBlock<TYPE>( 0 );
        const TYPE *xdata = x.getRawDataBlock<TYPE>( 0 );
        size_t N          = y.sizeOfDataBlock( 0 );
        auto lambda       = [] __device__( TYPE x ) { return x < 0 ? -x : x; };
        thrust::transform( thrust::device, xdata, xdata + N, data, lambda );
        cudaDeviceSynchronize();
    } else {
        // Default to VectorOperationsDefault (on cpu)
      VectorOperationsDefault<TYPE>::abs( x, y );
    }
}

template<typename TYPE>
void VectorOperationsCuda<TYPE>::addScalar( const VectorData &x, double alpha_in, VectorData &y )
{
  if ( checkData( x, y ) ) {
        TYPE *data        = y.getRawDataBlock<TYPE>( 0 );
        const TYPE *xdata = x.getRawDataBlock<TYPE>( 0 );
        size_t N          = y.sizeOfDataBlock( 0 );
        TYPE alpha        = static_cast<TYPE>( alpha_in );
        auto lambda       = [alpha] __device__( TYPE x ) { return x + alpha; };
        thrust::transform( thrust::device, xdata, xdata + N, data, lambda );
        cudaDeviceSynchronize();
    } else {
        // Default to VectorOperationsDefault (on cpu)
        VectorOperationsDefault<TYPE>::addScalar( x, alpha_in );
    }
}

template<typename TYPE>
double VectorOperationsCuda<TYPE>::localMin( const VectorData &x )  const
{
    double result = 0;
    if ( checkData(x) ) {
        const TYPE *data = x.getRawDataBlock<TYPE>( 0 );
        size_t N         = x.sizeOfDataBlock( 0 );
        result           = thrust::reduce( thrust::device,
                                 data,
                                 data + N,
                                 std::numeric_limits<TYPE>::max(),
                                 thrust::minimum<TYPE>() );
    } else {
        // Default to VectorOperationsDefault (on cpu)
        result = VectorOperationsDefault<TYPE>::localMin(x);
    }
    return result;
}

template<typename TYPE>
double VectorOperationsCuda<TYPE>::localMax( const VectorData &x )  const
{
    double result = 0;
    if ( checkData(x) ) {
        const TYPE *data = x.getRawDataBlock<TYPE>( 0 );
        size_t N         = x.sizeOfDataBlock( 0 );
        result           = thrust::reduce( thrust::device,
                                 data,
                                 data + N,
                                 -std::numeric_limits<TYPE>::max(),
                                 thrust::maximum<TYPE>() );
    } else {
        // Default to VectorOperationsDefault (on cpu)
        result = VectorOperationsDefault<TYPE>::localMax(x);
    }
    return result;
}

template<typename TYPE>
double VectorOperationsCuda<TYPE>::localL1Norm( const VectorData &x )  const
{
    double result = 0;
    if ( checkData(x) ) {
        const TYPE *data = x.getRawDataBlock<TYPE>( 0 );
        size_t N         = x.sizeOfDataBlock( 0 );
        auto lambda      = [=] __device__( TYPE x ) { return x < 0 ? -x : x; };
        result           = thrust::transform_reduce(
            thrust::device, data, data + N, lambda, (TYPE) 0, thrust::plus<TYPE>() );
    } else {
        // Default to VectorOperationsDefault (on cpu)
        result = VectorOperationsDefault<TYPE>::localL1Norm(x);
    }
    return result;
}

template<typename TYPE>
double VectorOperationsCuda<TYPE>::localL2Norm( const VectorData &x )  const
{
    double result = 0;
    if ( checkData(x) ) {
        const TYPE *data = x.getRawDataBlock<TYPE>( 0 );
        size_t N         = x.sizeOfDataBlock( 0 );
        auto lambda      = [=] __device__( TYPE x ) { return x * x; };
        result           = thrust::transform_reduce(
            thrust::device, data, data + N, lambda, (TYPE) 0, thrust::plus<TYPE>() );
        result = sqrt( result );
    } else {
        // Default to VectorOperationsDefault (on cpu)
        result = VectorOperationsDefault<TYPE>::localL2Norm(x);
    }
    return result;
}

template<typename TYPE>
double VectorOperationsCuda<TYPE>::localMaxNorm( const VectorData &x )  const
{
    double result = 0;
    if ( checkData(x) ) {
        const TYPE *data = x.getRawDataBlock<TYPE>( 0 );
        size_t N         = x.sizeOfDataBlock( 0 );
        auto lambda      = [=] __device__( TYPE x ) { return x < 0 ? -x : x; };
        result           = thrust::transform_reduce(
            thrust::device, data, data + N, lambda, (TYPE) 0, thrust::maximum<TYPE>() );
    } else {
        // Default to VectorOperationsDefault (on cpu)
        result = VectorOperationsDefault<TYPE>::localMaxNorm();
    }
    return result;
}

template<typename TYPE>
double VectorOperationsCuda<TYPE>::localDot( const VectorData &x, const VectorData &y ) const
{
    double result = 0;
    if ( checkData( x, y ) ) {
        const TYPE *data1 = y.getRawDataBlock<TYPE>( 0 );
        const TYPE *data2 = x.getRawDataBlock<TYPE>( 0 );
        size_t N          = y.sizeOfDataBlock( 0 );
        result = thrust::inner_product( thrust::device, data1, data1 + N, data2, (TYPE) 0 );
    } else {
        // Default to VectorOperationsDefault (on cpu)
      result = VectorOperationsDefault<TYPE>::localDot( x, y );
    }
    return result;
}

template<typename TYPE>
double VectorOperationsCuda<TYPE>::localMinQuotient( const VectorData &x, const VectorData &y ) const
{
    double result = 0;
    if ( checkData(x,y) ) {
        const TYPE *data1 = y.getRawDataBlock<TYPE>( 0 );
        const TYPE *data2 = x.getRawDataBlock<TYPE>( 0 );
        size_t N          = y.sizeOfDataBlock( 0 );
        result            = thrust::inner_product( thrust::device,
                                        data1,
                                        data1 + N,
                                        data2,
                                        std::numeric_limits<TYPE>::max(),
                                        thrust::minimum<TYPE>(),
                                        thrust::divides<TYPE>() );
    } else {
        // Default to VectorOperationsDefault (on cpu)
      result = VectorOperationsDefault<TYPE>::localMinQuotient( x, y );
    }
    return result;
}
template<typename T>
struct thrust_wrs {
    typedef T first_argument_type;
    typedef T second_argument_type;
    typedef T result_type;
    __host__ __device__ T operator()( const T &x, const T &y ) const { return x * x * y * y; }
};

template<typename TYPE>
double VectorOperationsCuda<TYPE>::localWrmsNorm( const VectorData &x, const VectorData &y ) const
{
    double result = 0;
    if ( checkData(x, y) ) {
        const TYPE *data1 = y.getRawDataBlock<TYPE>( 0 );
        const TYPE *data2 = x.getRawDataBlock<TYPE>( 0 );
        size_t N          = y.sizeOfDataBlock( 0 );
        result            = thrust::inner_product(
            thrust::device, data1, data1 + N, data2, 0, thrust::plus<TYPE>(), thrust_wrs<TYPE>() );
    } else {
        // Default to VectorOperationsDefault (on cpu)
      result = VectorOperationsDefault<TYPE>::localWrmsNorm( x, y );
    }
    return result;
}

template<typename TYPE>
double VectorOperationsCuda<TYPE>::localWrmsNormMask( const VectorData &x, const VectorData &mask, const VectorData &y ) const
{
    // Default to VectorOperationsDefault (on cpu)
  return VectorOperationsDefault<TYPE>::localWrmsNormMask( x, mask, y );
}

template<typename TYPE>
bool VectorOperationsCuda<TYPE>::localEquals( const VectorData &x, const VectorData &y, double tol ) const
{
  if ( checkData( x, y ) ) {
        // Call Cuda
      return VectorOperationsDefault<TYPE>::localEquals( x, y, tol );
    } else {
        // Default to VectorOperationsDefault (on cpu)
      return VectorOperationsDefault<TYPE>::localEquals( x, y, tol );
    }
}

} // namespace LinearAlgebra
} // namespace AMP

#endif
