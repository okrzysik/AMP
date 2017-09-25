#ifndef included_AMP_VectorOperationsCuda_hpp
#define included_AMP_VectorOperationsCuda_hpp

#include "vectors/Vector.h"
#include "vectors/data/VectorData.h"
#include "vectors/operations/VectorOperationsDefault.hpp"
#include "vectors/operations/cuda/VectorOperationsCuda.h"

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
template <typename TYPE>
AMP::shared_ptr<VectorOperations> VectorOperationsCuda<TYPE>::cloneOperations() const
{
    auto ptr = AMP::make_shared<VectorOperationsCuda<TYPE>>();
    return ptr;
}


/****************************************************************
* Check that all data can be passed to cuda                     *
****************************************************************/
template <typename TYPE>
bool VectorOperationsCuda<TYPE>::checkData() const
{
    return d_VectorData->numberOfDataBlocks() == 1;
}
template <typename TYPE>
bool VectorOperationsCuda<TYPE>::checkData( const VectorOperations &x ) const
{
    return x.getVectorData()->numberOfDataBlocks() == 1 && d_VectorData->numberOfDataBlocks() == 1;
}
template <typename TYPE>
bool VectorOperationsCuda<TYPE>::checkData( const VectorOperations &x,
                                            const VectorOperations &y ) const
{
    return x.getVectorData()->numberOfDataBlocks() == 1 &&
           y.getVectorData()->numberOfDataBlocks() == 1 && d_VectorData->numberOfDataBlocks() == 1;
}


/****************************************************************
* min, max, norms, etc.                                         *
****************************************************************/
template <typename TYPE>
bool VectorOperationsCuda<TYPE>::localEquals( const VectorOperations &rhs, double tol ) const
{
    if ( checkData( rhs ) ) {
        // Call Cuda
        return VectorOperationsDefault<TYPE>::localEquals( rhs, tol );
    } else {
        // Default to VectorOperationsDefault (on cpu)
        return VectorOperationsDefault<TYPE>::localEquals( rhs, tol );
    }
}
template <typename TYPE>
double VectorOperationsCuda<TYPE>::localMin( void ) const
{
    double result = 0;
    if ( checkData() ) {
        const TYPE *data = d_VectorData->getRawDataBlock<TYPE>( 0 );
        size_t N         = d_VectorData->sizeOfDataBlock( 0 );
        result           = thrust::reduce( thrust::device,
                                 data,
                                 data + N,
                                 std::numeric_limits<TYPE>::max(),
                                 thrust::minimum<TYPE>() );
    } else {
        // Default to VectorOperationsDefault (on cpu)
        result = VectorOperationsDefault<TYPE>::localMin();
    }
    return result;
}
template <typename TYPE>
double VectorOperationsCuda<TYPE>::localMax( void ) const
{
    double result = 0;
    if ( checkData() ) {
        const TYPE *data = d_VectorData->getRawDataBlock<TYPE>( 0 );
        size_t N         = d_VectorData->sizeOfDataBlock( 0 );
        result           = thrust::reduce( thrust::device,
                                 data,
                                 data + N,
                                 -std::numeric_limits<TYPE>::max(),
                                 thrust::maximum<TYPE>() );
    } else {
        // Default to VectorOperationsDefault (on cpu)
        result = VectorOperationsDefault<TYPE>::localMax();
    }
    return result;
}
template <typename TYPE>
double VectorOperationsCuda<TYPE>::localL1Norm( void ) const
{
    double result = 0;
    if ( checkData() ) {
        const TYPE *data = d_VectorData->getRawDataBlock<TYPE>( 0 );
        size_t N         = d_VectorData->sizeOfDataBlock( 0 );
        auto lambda = [=] __device__( TYPE x ) { return x < 0 ? -x : x; };
        result         = thrust::transform_reduce(
            thrust::device, data, data + N, lambda, (TYPE) 0, thrust::plus<TYPE>() );
    } else {
        // Default to VectorOperationsDefault (on cpu)
        result = VectorOperationsDefault<TYPE>::localL1Norm();
    }
    return result;
}
template <typename TYPE>
double VectorOperationsCuda<TYPE>::localL2Norm( void ) const
{
    double result = 0;
    if ( checkData() ) {
        const TYPE *data = d_VectorData->getRawDataBlock<TYPE>( 0 );
        size_t N         = d_VectorData->sizeOfDataBlock( 0 );
        auto lambda = [=] __device__( TYPE x ) { return x * x; };
        result         = thrust::transform_reduce(
            thrust::device, data, data + N, lambda, (TYPE) 0, thrust::plus<TYPE>() );
        result = sqrt( result );
    } else {
        // Default to VectorOperationsDefault (on cpu)
        result = VectorOperationsDefault<TYPE>::localL2Norm();
    }
    return result;
}
template <typename TYPE>
double VectorOperationsCuda<TYPE>::localMaxNorm( void ) const
{
    double result = 0;
    if ( checkData() ) {
        const TYPE *data = d_VectorData->getRawDataBlock<TYPE>( 0 );
        size_t N         = d_VectorData->sizeOfDataBlock( 0 );
        auto lambda = [=] __device__( TYPE x ) { return x < 0 ? -x : x; };
        result         = thrust::transform_reduce(
            thrust::device, data, data + N, lambda, (TYPE) 0, thrust::maximum<TYPE>() );
    } else {
        // Default to VectorOperationsDefault (on cpu)
        result = VectorOperationsDefault<TYPE>::localMaxNorm();
    }
    return result;
}
template <typename TYPE>
double VectorOperationsCuda<TYPE>::localDot( const VectorOperations &x ) const
{
    double result = 0;
    if ( checkData( x ) ) {
        const TYPE *data1 = d_VectorData->getRawDataBlock<TYPE>( 0 );
        const TYPE *data2 = x.getVectorData()->getRawDataBlock<TYPE>( 0 );
        size_t N          = d_VectorData->sizeOfDataBlock( 0 );
        result = thrust::inner_product( thrust::device, data1, data1 + N, data2, (TYPE) 0 );
    } else {
        // Default to VectorOperationsDefault (on cpu)
        result = VectorOperationsDefault<TYPE>::localDot( x );
    }
    return result;
}
template <typename TYPE>
double VectorOperationsCuda<TYPE>::localMinQuotient( const VectorOperations &x ) const
{
    double result = 0;
    if ( checkData() ) {
        const TYPE *data1 = d_VectorData->getRawDataBlock<TYPE>( 0 );
        const TYPE *data2 = x.getVectorData()->getRawDataBlock<TYPE>( 0 );
        size_t N          = d_VectorData->sizeOfDataBlock( 0 );
        result            = thrust::inner_product( thrust::device,
                                        data1,
                                        data1 + N,
                                        data2,
                                        std::numeric_limits<TYPE>::max(),
                                        thrust::minimum<TYPE>(),
                                        thrust::divides<TYPE>() );
    } else {
        // Default to VectorOperationsDefault (on cpu)
        result = VectorOperationsDefault<TYPE>::localMinQuotient( x );
    }
    return result;
}
template <typename T>
struct thrust_wrs {
    typedef T first_argument_type;
    typedef T second_argument_type;
    typedef T result_type;
    __host__ __device__ T operator()( const T &x, const T &y ) const { return x * x * y * y; }
};
template <typename TYPE>
double VectorOperationsCuda<TYPE>::localWrmsNorm( const VectorOperations &x ) const
{
    double result = 0;
    if ( checkData() ) {
        const TYPE *data1 = d_VectorData->getRawDataBlock<TYPE>( 0 );
        const TYPE *data2 = x.getVectorData()->getRawDataBlock<TYPE>( 0 );
        size_t N          = d_VectorData->sizeOfDataBlock( 0 );
        result            = thrust::inner_product(
            thrust::device, data1, data1 + N, data2, 0, thrust::plus<TYPE>(), thrust_wrs<TYPE>() );
    } else {
        // Default to VectorOperationsDefault (on cpu)
        result = VectorOperationsDefault<TYPE>::localWrmsNorm( x );
    }
    return result;
}
template <typename TYPE>
double VectorOperationsCuda<TYPE>::localWrmsNormMask( const VectorOperations &x,
                                                      const VectorOperations &mask ) const
{
    // Default to VectorOperationsDefault (on cpu)
    return VectorOperationsDefault<TYPE>::localWrmsNormMask( x, mask );
}


/****************************************************************
* Functions to initalize the data                               *
****************************************************************/
template <typename TYPE>
void VectorOperationsCuda<TYPE>::zero()
{
    VectorOperationsCuda<TYPE>::setToScalar( 0.0 );
}
template <typename TYPE>
void VectorOperationsCuda<TYPE>::setToScalar( double alpha )
{
    bool useGPU = checkData();
    if ( useGPU ) {
        const TYPE *data = d_VectorData->getRawDataBlock<TYPE>( 0 );
        size_t N         = d_VectorData->sizeOfDataBlock( 0 );
        thrust::fill_n( thrust::device, data, N, static_cast<TYPE>( alpha ) );
    } else {
        // Default to cpu version
        auto curMe = d_VectorData->begin<TYPE>();
        auto last  = d_VectorData->end<TYPE>();
        while ( curMe != last ) {
            *curMe = alpha;
            ++curMe;
        }
    }
    if ( hasGhosts() ) {
        auto &ghosts = getGhosts();
        for ( size_t i = 0; i != ghosts.size(); i++ )
            ghosts[i]  = alpha;
    }
    // Override the status state since we set the ghost values
    *( d_VectorData->getUpdateStatusPtr() ) = VectorData::UpdateState::UNCHANGED;
    // Wait for cuda data to complete
    if ( useGPU )
        cudaDeviceSynchronize();
}
template <typename TYPE>
void VectorOperationsCuda<TYPE>::setRandomValues()
{
    // Default to VectorOperationsDefault (on cpu)
    return VectorOperationsDefault<TYPE>::setRandomValues();
}
template <typename TYPE>
void VectorOperationsCuda<TYPE>::setRandomValues( RNG::shared_ptr rng )
{
    // Default to VectorOperationsDefault (on cpu)
    return VectorOperationsDefault<TYPE>::setRandomValues( rng );
}


/****************************************************************
* Basic linear algebra                                          *
****************************************************************/
template <typename TYPE>
void VectorOperationsCuda<TYPE>::copy( const VectorOperations &x )
{
    if ( checkData( x ) ) {
        TYPE *data        = d_VectorData->getRawDataBlock<TYPE>( 0 );
        const TYPE *xdata = x.getVectorData()->getRawDataBlock<TYPE>( 0 );
        size_t N          = d_VectorData->sizeOfDataBlock( 0 );
        thrust::copy_n( thrust::device, xdata, N, data );
        cudaDeviceSynchronize();
    } else {
        // Default to VectorOperationsDefault (on cpu)
        return VectorOperationsDefault<TYPE>::copy( x );
    }
}
template <typename TYPE>
void VectorOperationsCuda<TYPE>::scale( double alpha_in )
{
    if ( checkData() ) {
        TYPE *data  = d_VectorData->getRawDataBlock<TYPE>( 0 );
        size_t N    = d_VectorData->sizeOfDataBlock( 0 );
        TYPE alpha  = static_cast<TYPE>( alpha_in );
        auto lambda = [alpha] __device__( TYPE x ) { return x * alpha; };
        thrust::transform( thrust::device, data, data + N, data, lambda );
        cudaDeviceSynchronize();
    } else {
        // Default to VectorOperationsDefault (on cpu)
        VectorOperationsDefault<TYPE>::scale( alpha_in );
    }
}
template <typename TYPE>
void VectorOperationsCuda<TYPE>::scale( double alpha_in, const VectorOperations &x )
{
    if ( checkData( x ) ) {
        TYPE *data        = d_VectorData->getRawDataBlock<TYPE>( 0 );
        const TYPE *xdata = x.getVectorData()->getRawDataBlock<TYPE>( 0 );
        size_t N          = d_VectorData->sizeOfDataBlock( 0 );
        TYPE alpha        = static_cast<TYPE>( alpha_in );
        auto lambda       = [alpha] __device__( TYPE x ) { return x * alpha; };
        thrust::transform( thrust::device, xdata, xdata + N, data, lambda );
        cudaDeviceSynchronize();
    } else {
        // Default to VectorOperationsDefault (on cpu)
        VectorOperationsDefault<TYPE>::scale( alpha_in, x );
    }
}
template <typename TYPE>
void VectorOperationsCuda<TYPE>::add( const VectorOperations &x, const VectorOperations &y )
{
    if ( checkData( x, y ) ) {
        TYPE *data        = d_VectorData->getRawDataBlock<TYPE>( 0 );
        const TYPE *xdata = x.getVectorData()->getRawDataBlock<TYPE>( 0 );
        const TYPE *ydata = y.getVectorData()->getRawDataBlock<TYPE>( 0 );
        size_t N          = d_VectorData->sizeOfDataBlock( 0 );
        thrust::transform( thrust::device, xdata, xdata + N, ydata, data, thrust::plus<TYPE>() );
        cudaDeviceSynchronize();
    } else {
        // Default to VectorOperationsDefault (on cpu)
        VectorOperationsDefault<TYPE>::add( x, y );
    }
}
template <typename TYPE>
void VectorOperationsCuda<TYPE>::subtract( const VectorOperations &x, const VectorOperations &y )
{
    if ( checkData( x, y ) ) {
        TYPE *data        = d_VectorData->getRawDataBlock<TYPE>( 0 );
        const TYPE *xdata = x.getVectorData()->getRawDataBlock<TYPE>( 0 );
        const TYPE *ydata = y.getVectorData()->getRawDataBlock<TYPE>( 0 );
        size_t N          = d_VectorData->sizeOfDataBlock( 0 );
        thrust::transform( thrust::device, xdata, xdata + N, ydata, data, thrust::minus<TYPE>() );
        cudaDeviceSynchronize();
    } else {
        // Default to VectorOperationsDefault (on cpu)
        VectorOperationsDefault<TYPE>::subtract( x, y );
    }
}
template <typename TYPE>
void VectorOperationsCuda<TYPE>::multiply( const VectorOperations &x, const VectorOperations &y )
{
    if ( checkData( x, y ) ) {
        TYPE *data        = d_VectorData->getRawDataBlock<TYPE>( 0 );
        const TYPE *xdata = x.getVectorData()->getRawDataBlock<TYPE>( 0 );
        const TYPE *ydata = y.getVectorData()->getRawDataBlock<TYPE>( 0 );
        size_t N          = d_VectorData->sizeOfDataBlock( 0 );
        thrust::transform(
            thrust::device, xdata, xdata + N, ydata, data, thrust::multiplies<TYPE>() );
        cudaDeviceSynchronize();
    } else {
        // Default to VectorOperationsDefault (on cpu)
        VectorOperationsDefault<TYPE>::multiply( x, y );
    }
}
template <typename TYPE>
void VectorOperationsCuda<TYPE>::divide( const VectorOperations &x, const VectorOperations &y )
{
    if ( checkData( x, y ) ) {
        TYPE *data        = d_VectorData->getRawDataBlock<TYPE>( 0 );
        const TYPE *xdata = x.getVectorData()->getRawDataBlock<TYPE>( 0 );
        const TYPE *ydata = y.getVectorData()->getRawDataBlock<TYPE>( 0 );
        size_t N          = d_VectorData->sizeOfDataBlock( 0 );
        thrust::transform( thrust::device, xdata, xdata + N, ydata, data, thrust::divides<TYPE>() );
        cudaDeviceSynchronize();
    } else {
        // Default to VectorOperationsDefault (on cpu)
        VectorOperationsDefault<TYPE>::divide( x, y );
    }
}
template <typename TYPE>
void VectorOperationsCuda<TYPE>::reciprocal( const VectorOperations &x )
{
    if ( checkData( x ) ) {
        TYPE *data        = d_VectorData->getRawDataBlock<TYPE>( 0 );
        const TYPE *xdata = x.getVectorData()->getRawDataBlock<TYPE>( 0 );
        size_t N          = d_VectorData->sizeOfDataBlock( 0 );
        auto lambda       = [] __device__( TYPE x ) { return (TYPE) 1 / x; };
        thrust::transform( thrust::device, xdata, xdata + N, data, lambda );
        cudaDeviceSynchronize();
    } else {
        // Default to VectorOperationsDefault (on cpu)
        VectorOperationsDefault<TYPE>::reciprocal( x );
    }
}
template <typename TYPE>
void VectorOperationsCuda<TYPE>::linearSum( double alpha_in,
                                            const VectorOperations &x,
                                            double beta_in,
                                            const VectorOperations &y )
{
    if ( checkData( x, y ) ) {
        TYPE alpha        = alpha_in;
        TYPE beta         = beta_in;
        TYPE *data        = d_VectorData->getRawDataBlock<TYPE>( 0 );
        const TYPE *xdata = x.getVectorData()->getRawDataBlock<TYPE>( 0 );
        const TYPE *ydata = y.getVectorData()->getRawDataBlock<TYPE>( 0 );
        size_t N          = d_VectorData->sizeOfDataBlock( 0 );
        auto lambda = [alpha, beta] __device__( TYPE x, TYPE y ) { return alpha * x + beta * y; };
        thrust::transform( thrust::device, xdata, xdata + N, ydata, data, lambda );
        cudaDeviceSynchronize();
    } else {
        // Default to VectorOperationsDefault (on cpu)
        VectorOperationsDefault<TYPE>::linearSum( alpha_in, x, beta_in, y );
    }
}
template <typename TYPE>
void VectorOperationsCuda<TYPE>::axpy( double alpha_in,
                                       const VectorOperations &x,
                                       const VectorOperations &y )
{
    VectorOperationsCuda<TYPE>::linearSum( alpha_in, x, 1, y );
}
template <typename TYPE>
void VectorOperationsCuda<TYPE>::axpby( double alpha_in, double beta_in, const VectorOperations &x )
{
    VectorOperationsCuda<TYPE>::linearSum( alpha_in, x, beta_in, *this );
}
template <typename TYPE>
void VectorOperationsCuda<TYPE>::abs( const VectorOperations &x )
{
    if ( checkData( x ) ) {
        TYPE *data        = d_VectorData->getRawDataBlock<TYPE>( 0 );
        const TYPE *xdata = x.getVectorData()->getRawDataBlock<TYPE>( 0 );
        size_t N          = d_VectorData->sizeOfDataBlock( 0 );
        auto lambda       = [] __device__( TYPE x ) { return x < 0 ? -x : x; };
        thrust::transform( thrust::device, xdata, xdata + N, data, lambda );
        cudaDeviceSynchronize();
    } else {
        // Default to VectorOperationsDefault (on cpu)
        VectorOperationsDefault<TYPE>::abs( x );
    }
}
template <typename TYPE>
void VectorOperationsCuda<TYPE>::addScalar( const VectorOperations &x, double alpha_in )
{
    if ( checkData( x ) ) {
        TYPE *data        = d_VectorData->getRawDataBlock<TYPE>( 0 );
        const TYPE *xdata = x.getVectorData()->getRawDataBlock<TYPE>( 0 );
        size_t N          = d_VectorData->sizeOfDataBlock( 0 );
        TYPE alpha        = static_cast<TYPE>( alpha_in );
        auto lambda       = [alpha] __device__( TYPE x ) { return x + alpha; };
        thrust::transform( thrust::device, xdata, xdata + N, data, lambda );
        cudaDeviceSynchronize();
    } else {
        // Default to VectorOperationsDefault (on cpu)
        VectorOperationsDefault<TYPE>::addScalar( x, alpha_in );
    }
}


} // LinearAlgebra namespace
} // AMP namespace

#endif
