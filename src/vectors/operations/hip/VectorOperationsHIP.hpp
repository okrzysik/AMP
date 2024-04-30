#ifndef included_AMP_VectorOperationsHIP_hpp
#define included_AMP_VectorOperationsHIP_hpp

#include "AMP/utils/UtilityMacros.h"
#include "AMP/vectors/Vector.h"
#include "AMP/vectors/data/VectorData.h"
#include "AMP/vectors/operations/VectorOperationsDefault.hpp"
#include "AMP/vectors/operations/hip/HIPOperationsHelpers.h"
#include "AMP/vectors/operations/hip/VectorOperationsHIP.h"

#include <hip/hip_runtime_api.h>


namespace AMP {
namespace LinearAlgebra {


extern template class VectorOperationsHip<double>; // Suppresses implicit instantiation below --
extern template class VectorOperationsHip<float>;  // Suppresses implicit instantiation below --


/****************************************************************
 * Constructors                                                  *
 ****************************************************************/
template<typename TYPE>
std::shared_ptr<VectorOperations> VectorOperationsHip<TYPE>::cloneOperations() const
{
    return std::make_shared<VectorOperationsHip<TYPE>>();
}

template<typename TYPE>
VectorOperationsHip<TYPE>::~VectorOperationsHip<TYPE>()
{
    if ( d_default_ops )
        delete d_default_ops;
}


/****************************************************************
 * Check that all data can be passed to hip                     *
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
inline VectorOperationsDefault<TYPE> *VectorOperationsHip<TYPE>::getDefaultOps( void )
{
    if ( !d_default_ops )
        d_default_ops = new VectorOperationsDefault<TYPE>();
    return d_default_ops;
}

template<typename TYPE>
inline const VectorOperationsDefault<TYPE> *VectorOperationsHip<TYPE>::getDefaultOps( void ) const
{
    if ( !d_default_ops )
        d_default_ops = new VectorOperationsDefault<TYPE>();
    return d_default_ops;
}

//**********************************************************************
// Static functions that operate on VectorData objects

template<typename TYPE>
void VectorOperationsHip<TYPE>::zero( VectorData &x )
{
    VectorOperationsHip<TYPE>::setToScalar( 0.0, x );
    x.setUpdateStatus( UpdateState::UNCHANGED );
}

template<typename TYPE>
void VectorOperationsHip<TYPE>::setToScalar( const Scalar &alpha_in, VectorData &x )
{
    bool useGPU = checkData( x );
    TYPE alpha  = alpha_in.get<TYPE>();
    if ( useGPU ) {
        TYPE *data = x.getRawDataBlock<TYPE>( 0 );
        size_t N   = x.sizeOfDataBlock( 0 );
        HipOperationsHelpers<TYPE>::setToScalar( alpha, N, data );
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
    // Wait for hip data to complete
    if ( useGPU )
        hipDeviceSynchronize();
}

template<typename TYPE>
void VectorOperationsHip<TYPE>::setRandomValues( VectorData &x )
{
    // Default to VectorOperationsDefault (on cpu)
    getDefaultOps()->setRandomValues( x );
}

template<typename TYPE>
void VectorOperationsHip<TYPE>::copy( const VectorData &x, VectorData &y )
{
    if ( checkData<TYPE>( x, y ) ) {
        auto ydata = y.getRawDataBlock<TYPE>( 0 );
        auto xdata = x.getRawDataBlock<TYPE>( 0 );
        auto N     = y.sizeOfDataBlock( 0 );
        HipOperationsHelpers<TYPE>::copy( N, xdata, ydata );
        hipDeviceSynchronize();
    } else {
        // Default to VectorOperationsDefault (on cpu)
        getDefaultOps()->copy( x, y );
    }
}

template<typename TYPE>
void VectorOperationsHip<TYPE>::scale( const Scalar &alpha_in, VectorData &x )
{
    if ( checkData( x ) ) {
        TYPE *data = x.getRawDataBlock<TYPE>( 0 );
        size_t N   = x.sizeOfDataBlock( 0 );
        TYPE alpha = alpha_in.get<TYPE>();
        HipOperationsHelpers<TYPE>::scale( alpha, N, data );
        hipDeviceSynchronize();
    } else {
        // Default to VectorOperationsDefault (on cpu)
        getDefaultOps()->scale( alpha_in, x );
    }
}

template<typename TYPE>
void VectorOperationsHip<TYPE>::scale( const Scalar &alpha_in, const VectorData &x, VectorData &y )
{
    if ( checkData<TYPE>( x, y ) ) {
        auto xdata = x.getRawDataBlock<TYPE>( 0 );
        auto ydata = y.getRawDataBlock<TYPE>( 0 );
        auto N     = y.sizeOfDataBlock( 0 );
        auto alpha = alpha_in.get<TYPE>();
        HipOperationsHelpers<TYPE>::scale( alpha, N, xdata, ydata );
        hipDeviceSynchronize();
    } else {
        // Default to VectorOperationsDefault (on cpu)
        getDefaultOps()->scale( alpha_in, x, y );
    }
}

template<typename TYPE>
void VectorOperationsHip<TYPE>::add( const VectorData &x, const VectorData &y, VectorData &z )
{
    if ( checkData<TYPE>( x, y ) ) {
        auto xdata = x.getRawDataBlock<TYPE>( 0 );
        auto ydata = y.getRawDataBlock<TYPE>( 0 );
        auto zdata = z.getRawDataBlock<TYPE>( 0 );
        auto N     = z.sizeOfDataBlock( 0 );
        HipOperationsHelpers<TYPE>::add( N, xdata, ydata, zdata );
        hipDeviceSynchronize();
    } else {
        // Default to VectorOperationsDefault (on cpu)
        getDefaultOps()->add( x, y, z );
    }
}

template<typename TYPE>
void VectorOperationsHip<TYPE>::subtract( const VectorData &x, const VectorData &y, VectorData &z )
{
    if ( checkData<TYPE>( x, y, z ) ) {
        auto xdata = x.getRawDataBlock<TYPE>( 0 );
        auto ydata = y.getRawDataBlock<TYPE>( 0 );
        auto zdata = z.getRawDataBlock<TYPE>( 0 );
        size_t N   = z.sizeOfDataBlock( 0 );
        HipOperationsHelpers<TYPE>::subtract( N, xdata, ydata, zdata );
        hipDeviceSynchronize();
    } else {
        // Default to VectorOperationsDefault (on cpu)
        getDefaultOps()->subtract( x, y, z );
    }
}

template<typename TYPE>
void VectorOperationsHip<TYPE>::multiply( const VectorData &x, const VectorData &y, VectorData &z )
{
    if ( checkData<TYPE>( x, y, z ) ) {
        auto xdata = x.getRawDataBlock<TYPE>( 0 );
        auto ydata = y.getRawDataBlock<TYPE>( 0 );
        auto zdata = z.getRawDataBlock<TYPE>( 0 );
        size_t N   = z.sizeOfDataBlock( 0 );
        HipOperationsHelpers<TYPE>::multiply( N, xdata, ydata, zdata );
        hipDeviceSynchronize();
    } else {
        // Default to VectorOperationsDefault (on cpu)
        getDefaultOps()->multiply( x, y, z );
    }
}

template<typename TYPE>
void VectorOperationsHip<TYPE>::divide( const VectorData &x, const VectorData &y, VectorData &z )
{
    if ( checkData<TYPE>( x, y, z ) ) {
        auto xdata = x.getRawDataBlock<TYPE>( 0 );
        auto ydata = y.getRawDataBlock<TYPE>( 0 );
        auto zdata = z.getRawDataBlock<TYPE>( 0 );
        size_t N   = z.sizeOfDataBlock( 0 );
        HipOperationsHelpers<TYPE>::divide( N, xdata, ydata, zdata );
        hipDeviceSynchronize();
    } else {
        // Default to VectorOperationsDefault (on cpu)
        getDefaultOps()->divide( x, y, z );
    }
}


template<typename TYPE>
void VectorOperationsHip<TYPE>::reciprocal( const VectorData &x, VectorData &y )
{
    if ( checkData<TYPE>( x, y ) ) {
        auto xdata = x.getRawDataBlock<TYPE>( 0 );
        auto ydata = y.getRawDataBlock<TYPE>( 0 );
        size_t N   = y.sizeOfDataBlock( 0 );
        HipOperationsHelpers<TYPE>::reciprocal( N, xdata, ydata );
        hipDeviceSynchronize();
    } else {
        // Default to VectorOperationsDefault (on cpu)
        getDefaultOps()->reciprocal( x, y );
    }
}


template<typename TYPE>
void VectorOperationsHip<TYPE>::linearSum( const Scalar &alpha_in,
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
        HipOperationsHelpers<TYPE>::linearSum( alpha, N, xdata, beta, ydata, zdata );
        hipDeviceSynchronize();
    } else {
        // Default to VectorOperationsDefault (on cpu)
        getDefaultOps()->linearSum( alpha_in, x, beta_in, y, z );
    }
}

template<typename TYPE>
void VectorOperationsHip<TYPE>::axpy( const Scalar &alpha_in,
                                       const VectorData &x,
                                       const VectorData &y,
                                       VectorData &z )
{
    VectorOperationsHip<TYPE>::linearSum( alpha_in, x, 1.0, y, z );
}

template<typename TYPE>
void VectorOperationsHip<TYPE>::axpby( const Scalar &alpha_in,
                                        const Scalar &beta_in,
                                        const VectorData &x,
                                        VectorData &z )
{
    VectorOperationsHip<TYPE>::linearSum( alpha_in, x, beta_in, z, z );
}

template<typename TYPE>
void VectorOperationsHip<TYPE>::abs( const VectorData &x, VectorData &y )
{
    if ( checkData<TYPE>( x, y ) ) {
        auto xdata = x.getRawDataBlock<TYPE>( 0 );
        auto ydata = y.getRawDataBlock<TYPE>( 0 );
        size_t N   = y.sizeOfDataBlock( 0 );
        HipOperationsHelpers<TYPE>::abs( N, xdata, ydata );
        hipDeviceSynchronize();
    } else {
        // Default to VectorOperationsDefault (on cpu)
        getDefaultOps()->abs( x, y );
    }
}

template<typename TYPE>
void VectorOperationsHip<TYPE>::addScalar( const VectorData &x,
                                            const Scalar &alpha_in,
                                            VectorData &y )
{
    if ( checkData<TYPE>( x, y ) ) {
        auto xdata = x.getRawDataBlock<TYPE>( 0 );
        auto ydata = y.getRawDataBlock<TYPE>( 0 );
        size_t N   = y.sizeOfDataBlock( 0 );
        TYPE alpha = alpha_in.get<TYPE>();
        HipOperationsHelpers<TYPE>::addScalar( N, xdata, alpha, ydata );
        hipDeviceSynchronize();
    } else {
        // Default to VectorOperationsDefault (on cpu)
        getDefaultOps()->addScalar( x, alpha_in, y );
    }
}

template<typename TYPE>
Scalar VectorOperationsHip<TYPE>::localMin( const VectorData &x ) const
{
    if ( checkData( x ) ) {
        auto xdata = x.getRawDataBlock<TYPE>( 0 );
        size_t N   = x.sizeOfDataBlock( 0 );
        return HipOperationsHelpers<TYPE>::localMin( N, xdata );
    } else {
        // Default to VectorOperationsDefault (on cpu)
        return getDefaultOps()->localMin( x );
    }
}

template<typename TYPE>
Scalar VectorOperationsHip<TYPE>::localMax( const VectorData &x ) const
{
    if ( checkData( x ) ) {
        auto xdata = x.getRawDataBlock<TYPE>( 0 );
        size_t N   = x.sizeOfDataBlock( 0 );
        return HipOperationsHelpers<TYPE>::localMax( N, xdata );
    } else {
        // Default to VectorOperationsDefault (on cpu)
        return getDefaultOps()->localMax( x );
    }
}

template<typename TYPE>
Scalar VectorOperationsHip<TYPE>::localSum( const VectorData &x ) const
{
    if ( checkData( x ) ) {
        auto xdata = x.getRawDataBlock<TYPE>( 0 );
        size_t N   = x.sizeOfDataBlock( 0 );
        return HipOperationsHelpers<TYPE>::localSum( N, xdata );
    } else {
        // Default to VectorOperationsDefault (on cpu)
        return getDefaultOps()->localSum( x );
    }
}

template<typename TYPE>
Scalar VectorOperationsHip<TYPE>::localL1Norm( const VectorData &x ) const
{
    if ( checkData( x ) ) {
        auto xdata = x.getRawDataBlock<TYPE>( 0 );
        size_t N   = x.sizeOfDataBlock( 0 );
        return HipOperationsHelpers<TYPE>::localL1Norm( N, xdata );
    } else {
        // Default to VectorOperationsDefault (on cpu)
        return getDefaultOps()->localL1Norm( x );
    }
}

template<typename TYPE>
Scalar VectorOperationsHip<TYPE>::localL2Norm( const VectorData &x ) const
{
    if ( checkData( x ) ) {
        auto xdata = x.getRawDataBlock<TYPE>( 0 );
        size_t N   = x.sizeOfDataBlock( 0 );
        return HipOperationsHelpers<TYPE>::localL2Norm( N, xdata );
    } else {
        // Default to VectorOperationsDefault (on cpu)
        return getDefaultOps()->localL2Norm( x );
    }
}

template<typename TYPE>
Scalar VectorOperationsHip<TYPE>::localMaxNorm( const VectorData &x ) const
{
    if ( checkData( x ) ) {
        auto xdata = x.getRawDataBlock<TYPE>( 0 );
        size_t N   = x.sizeOfDataBlock( 0 );
        return HipOperationsHelpers<TYPE>::localMaxNorm( N, xdata );
    } else {
        // Default to VectorOperationsDefault (on cpu)
        return getDefaultOps()->localMaxNorm( x );
    }
}

template<typename TYPE>
Scalar VectorOperationsHip<TYPE>::localDot( const VectorData &x, const VectorData &y ) const
{
    if ( checkData<TYPE>( x, y ) ) {
        auto xdata = x.getRawDataBlock<TYPE>( 0 );
        auto ydata = y.getRawDataBlock<TYPE>( 0 );
        size_t N   = x.sizeOfDataBlock( 0 );
        return HipOperationsHelpers<TYPE>::localDot( N, xdata, ydata );
    } else {
        // Default to VectorOperationsDefault (on cpu)
        return getDefaultOps()->localDot( x, y );
    }
}

template<typename TYPE>
Scalar VectorOperationsHip<TYPE>::localMinQuotient( const VectorData &x,
                                                     const VectorData &y ) const
{
    if ( checkData<TYPE>( x, y ) ) {
        auto xdata = x.getRawDataBlock<TYPE>( 0 );
        auto ydata = y.getRawDataBlock<TYPE>( 0 );
        size_t N   = x.sizeOfDataBlock( 0 );
        return HipOperationsHelpers<TYPE>::localMinQuotient( N, xdata, ydata );
    } else {
        // Default to VectorOperationsDefault (on cpu)
        return getDefaultOps()->localMinQuotient( x, y );
    }
}

template<typename TYPE>
Scalar VectorOperationsHip<TYPE>::localWrmsNorm( const VectorData &x, const VectorData &y ) const
{
    if ( checkData<TYPE>( x, y ) ) {
        auto xdata = x.getRawDataBlock<TYPE>( 0 );
        auto ydata = y.getRawDataBlock<TYPE>( 0 );
        size_t N   = x.sizeOfDataBlock( 0 );
        return HipOperationsHelpers<TYPE>::localWrmsNorm( N, xdata, ydata );
    } else {
        // Default to VectorOperationsDefault (on cpu)
        return getDefaultOps()->localWrmsNorm( x, y );
    }
}

template<typename TYPE>
Scalar VectorOperationsHip<TYPE>::localWrmsNormMask( const VectorData &x,
                                                      const VectorData &mask,
                                                      const VectorData &y ) const
{
    // Default to VectorOperationsDefault (on cpu)
    return getDefaultOps()->localWrmsNormMask( x, mask, y );
}

template<typename TYPE>
bool VectorOperationsHip<TYPE>::localEquals( const VectorData &x,
                                              const VectorData &y,
                                              const Scalar &tol_in ) const
{
    TYPE tol = tol_in.get<TYPE>();
    if ( checkData<TYPE>( x, y ) ) {
        // Call Hip
        return getDefaultOps()->localEquals( x, y, tol );
    } else {
        // Default to VectorOperationsDefault (on cpu)
        return getDefaultOps()->localEquals( x, y, tol );
    }
}

} // namespace LinearAlgebra
} // namespace AMP

#endif
