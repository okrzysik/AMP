#ifndef included_AMP_VectorOperationsCuda_hpp
#define included_AMP_VectorOperationsCuda_hpp

#include "AMP/utils/UtilityMacros.h"
#include "AMP/vectors/Vector.h"
#include "AMP/vectors/data/VectorData.h"
#include "AMP/vectors/operations/VectorOperationsDefault.hpp"
#include "AMP/vectors/operations/cuda/CudaOperationsHelpers.h"
#include "AMP/vectors/operations/cuda/VectorOperationsCuda.h"

#include <cuda.h>
#include <cuda_runtime_api.h>


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
    return std::make_shared<VectorOperationsCuda<TYPE>>();
}

template<typename TYPE>
VectorOperationsCuda<TYPE>::~VectorOperationsCuda<TYPE>()
{
    if ( d_default_ops )
        delete d_default_ops;
}


/****************************************************************
 * Check that all data can be passed to cuda                     *
 ****************************************************************/
inline bool checkData( const VectorData &x ) { return x.numberOfDataBlocks() == 1; }
template<class TYPE>
inline bool checkData( const VectorData &x, const VectorData &y )
{
    constexpr auto type = getTypeID<TYPE>();
    return x.numberOfDataBlocks() == 1 && y.numberOfDataBlocks() == 1 && y.getType( 0 ) == type;
}
template<class TYPE>
inline bool checkData( const VectorData &x, const VectorData &y, const VectorData &z )
{
    constexpr auto type = getTypeID<TYPE>();
    return x.numberOfDataBlocks() == 1 && y.numberOfDataBlocks() == 1 &&
           z.numberOfDataBlocks() == 1 && y.getType( 0 ) == type && z.getType( 0 ) == type;
}

template<typename TYPE>
inline VectorOperationsDefault<TYPE> *VectorOperationsCuda<TYPE>::getDefaultOps( void )
{
    if ( !d_default_ops )
        d_default_ops = new VectorOperationsDefault<TYPE>();
    return d_default_ops;
}

template<typename TYPE>
inline const VectorOperationsDefault<TYPE> *VectorOperationsCuda<TYPE>::getDefaultOps( void ) const
{
    if ( !d_default_ops )
        d_default_ops = new VectorOperationsDefault<TYPE>();
    return d_default_ops;
}

//**********************************************************************
// Static functions that operate on VectorData objects

template<typename TYPE>
void VectorOperationsCuda<TYPE>::zero( VectorData &x )
{
    VectorOperationsCuda<TYPE>::setToScalar( 0.0, x );
    x.setUpdateStatus( UpdateState::UNCHANGED );
}

template<typename TYPE>
void VectorOperationsCuda<TYPE>::setToScalar( const Scalar &alpha_in, VectorData &x )
{
    bool useGPU = checkData( x );
    TYPE alpha  = alpha_in.get<TYPE>();
    if ( useGPU ) {
        TYPE *data = x.getRawDataBlock<TYPE>( 0 );
        size_t N   = x.sizeOfDataBlock( 0 );
        CudaOperationsHelpers<TYPE>::setToScalar( alpha, N, data );
    } else {
        // Default to cpu version
        auto curMe = x.begin<TYPE>();
        auto last  = x.end<TYPE>();
        while ( curMe != last ) {
            *curMe = alpha;
            ++curMe;
        }
    }
    if ( x.hasGhosts() ) {
        auto &ghosts = x.getGhosts();
        for ( size_t i = 0; i != ghosts.size(); i++ )
            ghosts[i] = alpha;
    }
    // Override the status state since we set the ghost values
    x.setUpdateStatus( UpdateState::UNCHANGED );
    // Wait for cuda data to complete
    if ( useGPU )
        cudaDeviceSynchronize();
}

template<typename TYPE>
void VectorOperationsCuda<TYPE>::setRandomValues( VectorData &x )
{
    // Default to VectorOperationsDefault (on cpu)
    getDefaultOps()->setRandomValues( x );
}

template<typename TYPE>
void VectorOperationsCuda<TYPE>::copy( const VectorData &x, VectorData &y )
{
    if ( checkData<TYPE>( x, y ) ) {
        auto ydata = y.getRawDataBlock<TYPE>( 0 );
        auto xdata = x.getRawDataBlock<TYPE>( 0 );
        auto N     = y.sizeOfDataBlock( 0 );
        CudaOperationsHelpers<TYPE>::copy( N, xdata, ydata );
        cudaDeviceSynchronize();
    } else {
        // Default to VectorOperationsDefault (on cpu)
        getDefaultOps()->copy( x, y );
    }
}

template<typename TYPE>
void VectorOperationsCuda<TYPE>::scale( const Scalar &alpha_in, VectorData &x )
{
    if ( checkData( x ) ) {
        TYPE *data = x.getRawDataBlock<TYPE>( 0 );
        size_t N   = x.sizeOfDataBlock( 0 );
        TYPE alpha = alpha_in.get<TYPE>();
        CudaOperationsHelpers<TYPE>::scale( alpha, N, data );
        cudaDeviceSynchronize();
    } else {
        // Default to VectorOperationsDefault (on cpu)
        getDefaultOps()->scale( alpha_in, x );
    }
}

template<typename TYPE>
void VectorOperationsCuda<TYPE>::scale( const Scalar &alpha_in, const VectorData &x, VectorData &y )
{
    if ( checkData<TYPE>( x, y ) ) {
        auto xdata = x.getRawDataBlock<TYPE>( 0 );
        auto ydata = y.getRawDataBlock<TYPE>( 0 );
        auto N     = y.sizeOfDataBlock( 0 );
        auto alpha = alpha_in.get<TYPE>();
        CudaOperationsHelpers<TYPE>::scale( alpha, N, xdata, ydata );
        cudaDeviceSynchronize();
    } else {
        // Default to VectorOperationsDefault (on cpu)
        getDefaultOps()->scale( alpha_in, x, y );
    }
}

template<typename TYPE>
void VectorOperationsCuda<TYPE>::add( const VectorData &x, const VectorData &y, VectorData &z )
{
    if ( checkData<TYPE>( x, y ) ) {
        auto xdata = x.getRawDataBlock<TYPE>( 0 );
        auto ydata = y.getRawDataBlock<TYPE>( 0 );
        auto zdata = z.getRawDataBlock<TYPE>( 0 );
        auto N     = z.sizeOfDataBlock( 0 );
        CudaOperationsHelpers<TYPE>::add( N, xdata, ydata, zdata );
        cudaDeviceSynchronize();
    } else {
        // Default to VectorOperationsDefault (on cpu)
        getDefaultOps()->add( x, y, z );
    }
}

template<typename TYPE>
void VectorOperationsCuda<TYPE>::subtract( const VectorData &x, const VectorData &y, VectorData &z )
{
    if ( checkData<TYPE>( x, y, z ) ) {
        auto xdata = x.getRawDataBlock<TYPE>( 0 );
        auto ydata = y.getRawDataBlock<TYPE>( 0 );
        auto zdata = z.getRawDataBlock<TYPE>( 0 );
        size_t N   = z.sizeOfDataBlock( 0 );
        CudaOperationsHelpers<TYPE>::subtract( N, xdata, ydata, zdata );
        cudaDeviceSynchronize();
    } else {
        // Default to VectorOperationsDefault (on cpu)
        getDefaultOps()->subtract( x, y, z );
    }
}

template<typename TYPE>
void VectorOperationsCuda<TYPE>::multiply( const VectorData &x, const VectorData &y, VectorData &z )
{
    if ( checkData<TYPE>( x, y, z ) ) {
        auto xdata = x.getRawDataBlock<TYPE>( 0 );
        auto ydata = y.getRawDataBlock<TYPE>( 0 );
        auto zdata = z.getRawDataBlock<TYPE>( 0 );
        size_t N   = z.sizeOfDataBlock( 0 );
        CudaOperationsHelpers<TYPE>::multiply( N, xdata, ydata, zdata );
        cudaDeviceSynchronize();
    } else {
        // Default to VectorOperationsDefault (on cpu)
        getDefaultOps()->multiply( x, y, z );
    }
}

template<typename TYPE>
void VectorOperationsCuda<TYPE>::divide( const VectorData &x, const VectorData &y, VectorData &z )
{
    if ( checkData<TYPE>( x, y, z ) ) {
        auto xdata = x.getRawDataBlock<TYPE>( 0 );
        auto ydata = y.getRawDataBlock<TYPE>( 0 );
        auto zdata = z.getRawDataBlock<TYPE>( 0 );
        size_t N   = z.sizeOfDataBlock( 0 );
        CudaOperationsHelpers<TYPE>::divide( N, xdata, ydata, zdata );
        cudaDeviceSynchronize();
    } else {
        // Default to VectorOperationsDefault (on cpu)
        getDefaultOps()->divide( x, y, z );
    }
}


template<typename TYPE>
void VectorOperationsCuda<TYPE>::reciprocal( const VectorData &x, VectorData &y )
{
    if ( checkData<TYPE>( x, y ) ) {
        auto xdata = x.getRawDataBlock<TYPE>( 0 );
        auto ydata = y.getRawDataBlock<TYPE>( 0 );
        size_t N   = y.sizeOfDataBlock( 0 );
        CudaOperationsHelpers<TYPE>::reciprocal( N, xdata, ydata );
        cudaDeviceSynchronize();
    } else {
        // Default to VectorOperationsDefault (on cpu)
        getDefaultOps()->reciprocal( x, y );
    }
}


template<typename TYPE>
void VectorOperationsCuda<TYPE>::linearSum( const Scalar &alpha_in,
                                            const VectorData &x,
                                            const Scalar &beta_in,
                                            const VectorData &y,
                                            VectorData &z )
{
    if ( checkData<TYPE>( x, y, z ) ) {
        TYPE alpha = alpha_in.get<TYPE>();
        TYPE beta  = beta_in.get<TYPE>();
        auto xdata = x.getRawDataBlock<TYPE>( 0 );
        auto ydata = y.getRawDataBlock<TYPE>( 0 );
        auto zdata = z.getRawDataBlock<TYPE>( 0 );
        size_t N   = z.sizeOfDataBlock( 0 );
        CudaOperationsHelpers<TYPE>::linearSum( alpha, N, xdata, beta, ydata, zdata );
        cudaDeviceSynchronize();
    } else {
        // Default to VectorOperationsDefault (on cpu)
        getDefaultOps()->linearSum( alpha_in, x, beta_in, y, z );
    }
}

template<typename TYPE>
void VectorOperationsCuda<TYPE>::axpy( const Scalar &alpha_in,
                                       const VectorData &x,
                                       const VectorData &y,
                                       VectorData &z )
{
    VectorOperationsCuda<TYPE>::linearSum( alpha_in, x, 1.0, y, z );
}

template<typename TYPE>
void VectorOperationsCuda<TYPE>::axpby( const Scalar &alpha_in,
                                        const Scalar &beta_in,
                                        const VectorData &x,
                                        VectorData &z )
{
    VectorOperationsCuda<TYPE>::linearSum( alpha_in, x, beta_in, z, z );
}

template<typename TYPE>
void VectorOperationsCuda<TYPE>::abs( const VectorData &x, VectorData &y )
{
    if ( checkData<TYPE>( x, y ) ) {
        auto xdata = x.getRawDataBlock<TYPE>( 0 );
        auto ydata = y.getRawDataBlock<TYPE>( 0 );
        size_t N   = y.sizeOfDataBlock( 0 );
        CudaOperationsHelpers<TYPE>::abs( N, xdata, ydata );
        cudaDeviceSynchronize();
    } else {
        // Default to VectorOperationsDefault (on cpu)
        getDefaultOps()->abs( x, y );
    }
}

template<typename TYPE>
void VectorOperationsCuda<TYPE>::addScalar( const VectorData &x,
                                            const Scalar &alpha_in,
                                            VectorData &y )
{
    if ( checkData<TYPE>( x, y ) ) {
        auto xdata = x.getRawDataBlock<TYPE>( 0 );
        auto ydata = y.getRawDataBlock<TYPE>( 0 );
        size_t N   = y.sizeOfDataBlock( 0 );
        TYPE alpha = alpha_in.get<TYPE>();
        CudaOperationsHelpers<TYPE>::addScalar( N, xdata, alpha, ydata );
        cudaDeviceSynchronize();
    } else {
        // Default to VectorOperationsDefault (on cpu)
        getDefaultOps()->addScalar( x, alpha_in, y );
    }
}

template<typename TYPE>
Scalar VectorOperationsCuda<TYPE>::localMin( const VectorData &x ) const
{
    if ( checkData( x ) ) {
        auto xdata = x.getRawDataBlock<TYPE>( 0 );
        size_t N   = x.sizeOfDataBlock( 0 );
        return CudaOperationsHelpers<TYPE>::localMin( N, xdata );
    } else {
        // Default to VectorOperationsDefault (on cpu)
        return getDefaultOps()->localMin( x );
    }
}

template<typename TYPE>
Scalar VectorOperationsCuda<TYPE>::localMax( const VectorData &x ) const
{
    if ( checkData( x ) ) {
        auto xdata = x.getRawDataBlock<TYPE>( 0 );
        size_t N   = x.sizeOfDataBlock( 0 );
        return CudaOperationsHelpers<TYPE>::localMax( N, xdata );
    } else {
        // Default to VectorOperationsDefault (on cpu)
        return getDefaultOps()->localMax( x );
    }
}

template<typename TYPE>
Scalar VectorOperationsCuda<TYPE>::localSum( const VectorData &x ) const
{
    if ( checkData( x ) ) {
        auto xdata = x.getRawDataBlock<TYPE>( 0 );
        size_t N   = x.sizeOfDataBlock( 0 );
        return CudaOperationsHelpers<TYPE>::localSum( N, xdata );
    } else {
        // Default to VectorOperationsDefault (on cpu)
        return getDefaultOps()->localSum( x );
    }
}

template<typename TYPE>
Scalar VectorOperationsCuda<TYPE>::localL1Norm( const VectorData &x ) const
{
    if ( checkData( x ) ) {
        auto xdata = x.getRawDataBlock<TYPE>( 0 );
        size_t N   = x.sizeOfDataBlock( 0 );
        return CudaOperationsHelpers<TYPE>::localL1Norm( N, xdata );
    } else {
        // Default to VectorOperationsDefault (on cpu)
        return getDefaultOps()->localL1Norm( x );
    }
}

template<typename TYPE>
Scalar VectorOperationsCuda<TYPE>::localL2Norm( const VectorData &x ) const
{
    if ( checkData( x ) ) {
        auto xdata = x.getRawDataBlock<TYPE>( 0 );
        size_t N   = x.sizeOfDataBlock( 0 );
        return CudaOperationsHelpers<TYPE>::localL2Norm( N, xdata );
    } else {
        // Default to VectorOperationsDefault (on cpu)
        return getDefaultOps()->localL2Norm( x );
    }
}

template<typename TYPE>
Scalar VectorOperationsCuda<TYPE>::localMaxNorm( const VectorData &x ) const
{
    if ( checkData( x ) ) {
        auto xdata = x.getRawDataBlock<TYPE>( 0 );
        size_t N   = x.sizeOfDataBlock( 0 );
        return CudaOperationsHelpers<TYPE>::localMaxNorm( N, xdata );
    } else {
        // Default to VectorOperationsDefault (on cpu)
        return getDefaultOps()->localMaxNorm( x );
    }
}

template<typename TYPE>
Scalar VectorOperationsCuda<TYPE>::localDot( const VectorData &x, const VectorData &y ) const
{
    if ( checkData<TYPE>( x, y ) ) {
        auto xdata = x.getRawDataBlock<TYPE>( 0 );
        auto ydata = y.getRawDataBlock<TYPE>( 0 );
        size_t N   = x.sizeOfDataBlock( 0 );
        return CudaOperationsHelpers<TYPE>::localDot( N, xdata, ydata );
    } else {
        // Default to VectorOperationsDefault (on cpu)
        return getDefaultOps()->localDot( x, y );
    }
}

template<typename TYPE>
Scalar VectorOperationsCuda<TYPE>::localMinQuotient( const VectorData &x,
                                                     const VectorData &y ) const
{
    if ( checkData<TYPE>( x, y ) ) {
        auto xdata = x.getRawDataBlock<TYPE>( 0 );
        auto ydata = y.getRawDataBlock<TYPE>( 0 );
        size_t N   = x.sizeOfDataBlock( 0 );
        return CudaOperationsHelpers<TYPE>::localMinQuotient( N, xdata, ydata );
    } else {
        // Default to VectorOperationsDefault (on cpu)
        return getDefaultOps()->localMinQuotient( x, y );
    }
}

template<typename TYPE>
Scalar VectorOperationsCuda<TYPE>::localWrmsNorm( const VectorData &x, const VectorData &y ) const
{
    if ( checkData<TYPE>( x, y ) ) {
        auto xdata = x.getRawDataBlock<TYPE>( 0 );
        auto ydata = y.getRawDataBlock<TYPE>( 0 );
        size_t N   = x.sizeOfDataBlock( 0 );
        return CudaOperationsHelpers<TYPE>::localWrmsNorm( N, xdata, ydata );
    } else {
        // Default to VectorOperationsDefault (on cpu)
        return getDefaultOps()->localWrmsNorm( x, y );
    }
}

template<typename TYPE>
Scalar VectorOperationsCuda<TYPE>::localWrmsNormMask( const VectorData &x,
                                                      const VectorData &mask,
                                                      const VectorData &y ) const
{
    // Default to VectorOperationsDefault (on cpu)
    return getDefaultOps()->localWrmsNormMask( x, mask, y );
}

template<typename TYPE>
bool VectorOperationsCuda<TYPE>::localEquals( const VectorData &x,
                                              const VectorData &y,
                                              const Scalar &tol_in ) const
{
    TYPE tol = tol_in.get<TYPE>();
    if ( checkData<TYPE>( x, y ) ) {
        // Call Cuda
        return getDefaultOps()->localEquals( x, y, tol );
    } else {
        // Default to VectorOperationsDefault (on cpu)
        return getDefaultOps()->localEquals( x, y, tol );
    }
}

} // namespace LinearAlgebra
} // namespace AMP

#endif
