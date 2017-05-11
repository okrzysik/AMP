#ifndef included_AMP_VectorOperationsDefault_hpp
#define included_AMP_VectorOperationsDefault_hpp

#include "vectors/operations/VectorOperationsDefault.h"
#include "vectors/VectorData.h"
#include "vectors/Vector.h"


namespace AMP {
namespace LinearAlgebra {


extern template class VectorOperationsDefault<double>; // Suppresses implicit instantiation below --


/****************************************************************
* min, max, norms, etc.                                         *
****************************************************************/
template<typename TYPE>
bool VectorOperationsDefault<TYPE>::localEquals( const VectorOperations &rhs, double tol ) const
{
    const auto& x = *d_VectorData;
    const auto& y = *rhs.getVectorData();
    if ( ( x.getGlobalSize() != y.getGlobalSize() ) || ( x.getLocalSize() != y.getLocalSize() ) )
        return false;
    bool equal = true;
    auto cur1 = x.begin<TYPE>();
    auto cur2 = y.begin<TYPE>();
    auto last = x.end<TYPE>();
    while ( cur1 != last ) {
        if ( fabs( *cur1 - *cur2 ) > tol ) {
            equal = false;
            break;
        }
        ++cur1;
        ++cur2;
    }
    return equal;
}
template<typename TYPE>
double VectorOperationsDefault<TYPE>::localMin( void ) const
{
    size_t N_blocks = d_VectorData->numberOfDataBlocks();
    TYPE ans        = std::numeric_limits<TYPE>::max();
    for ( size_t i = 0; i < N_blocks; i++ ) {
        size_t size      = d_VectorData->sizeOfDataBlock( i );
        const TYPE *data = d_VectorData->getRawDataBlock<TYPE>( i );
        for ( size_t j = 0; j < size; j++ )
            ans = std::min( data[j], ans );
    }
    return static_cast<double>(ans);
}
template<typename TYPE>
double VectorOperationsDefault<TYPE>::localMax( void ) const
{
    size_t N_blocks = d_VectorData->numberOfDataBlocks();
    TYPE ans        = std::numeric_limits<TYPE>::lowest();
    for ( size_t i = 0; i < N_blocks; i++ ) {
        size_t size      = d_VectorData->sizeOfDataBlock( i );
        const TYPE *data = d_VectorData->getRawDataBlock<TYPE>( i );
        for ( size_t j = 0; j < size; j++ )
            ans = std::max( data[j], ans );
    }
    return static_cast<double>(ans);
}
template<typename TYPE>
double VectorOperationsDefault<TYPE>::localL1Norm( void ) const
{
    size_t N_blocks = d_VectorData->numberOfDataBlocks();
    double ans      = 0;
    for ( size_t i = 0; i < N_blocks; i++ ) {
        size_t size      = d_VectorData->sizeOfDataBlock( i );
        const TYPE *data = d_VectorData->getRawDataBlock<TYPE>( i );
        for ( size_t j = 0; j < size; j++ )
            ans += static_cast<double>( std::abs( data[j] ) );
    }
    return ans;
}
template<typename TYPE>
double VectorOperationsDefault<TYPE>::localL2Norm( void ) const
{
    size_t N_blocks = d_VectorData->numberOfDataBlocks();
    double ans      = 0;
    for ( size_t i = 0; i < N_blocks; i++ ) {
        size_t size      = d_VectorData->sizeOfDataBlock( i );
        const TYPE *data = d_VectorData->getRawDataBlock<TYPE>( i );
        for ( size_t j = 0; j < size; j++ ) {
            double tmp = static_cast<double>(data[j]);
            ans += tmp*tmp;
        }
    }
    return sqrt( ans );
}
template<typename TYPE>
double VectorOperationsDefault<TYPE>::localMaxNorm( void ) const
{
    size_t N_blocks = d_VectorData->numberOfDataBlocks();
    TYPE ans        = 0;
    for ( size_t i = 0; i < N_blocks; i++ ) {
        size_t size      = d_VectorData->sizeOfDataBlock( i );
        const TYPE *data = d_VectorData->getRawDataBlock<TYPE>( i );
        for ( size_t j = 0; j < size; j++ )
            ans = std::max( std::abs( data[j] ), ans );
    }
    return static_cast<double>(ans);
}
template<typename TYPE>
double VectorOperationsDefault<TYPE>::localDot( const VectorOperations& x ) const
{
    AMP_ASSERT( d_VectorData->getLocalSize() == x.getVectorData()->getLocalSize() );
    auto curMe   = d_VectorData->begin<TYPE>();
    auto last    = d_VectorData->end<TYPE>();
    auto curXRhs = x.getVectorData()->begin<TYPE>();
    double ans   = 0;
    while ( curMe != last ) {
        double v1 = static_cast<double>( *curMe );
        double v2 = static_cast<double>( *curXRhs );
        ans += v1*v2;
        ++curXRhs;
        ++curMe;
    }
    return ans;
}
template<typename TYPE>
double VectorOperationsDefault<TYPE>::localMinQuotient( const VectorOperations &x ) const
{
    auto curx = x.getVectorData()->begin<TYPE>();
    auto endx = x.getVectorData()->end<TYPE>();
    auto cury = d_VectorData->begin<TYPE>();
    double ans = std::numeric_limits<double>::max();
    while ( curx != endx ) {
        if ( *cury != 0 ) {
            double v1 = static_cast<double>( *curx );
            double v2 = static_cast<double>( *cury );
            ans = std::min( ans, v1 / v2 );
        }
        ++curx;
        ++cury;
    }
    return ans;
}
template<typename TYPE>
double VectorOperationsDefault<TYPE>::localWrmsNorm( const VectorOperations &x ) const
{
    auto curx = x.getVectorData()->begin<TYPE>();
    auto endx = x.getVectorData()->end<TYPE>();
    auto cury = d_VectorData->begin<TYPE>();
    double ans = 0;
    size_t N = 0;
    while ( curx != endx ) {
        double v1 = static_cast<double>( *curx );
        double v2 = static_cast<double>( *cury );
        ans += v1*v1 * v2*v2;
        ++curx;
        ++cury;
        ++N;
    }
    return sqrt(ans/N);
}
template<typename TYPE>
double VectorOperationsDefault<TYPE>::localWrmsNormMask( const VectorOperations &x,
                             const VectorOperations &mask ) const
{
    auto curx = x.getVectorData()->begin<TYPE>();
    auto endx = x.getVectorData()->end<TYPE>();
    auto cury = d_VectorData->begin<TYPE>();
    auto curm = mask.getVectorData()->begin<TYPE>();
    double ans = 0;
    size_t N = 0;
    while ( curx != endx ) {
        if ( *curm > 0 ) {
            double v1 = static_cast<double>( *curx );
            double v2 = static_cast<double>( *cury );
            ans += v1*v1 * v2*v2;
        }
        ++curx;
        ++cury;
        ++curm;
        ++N;
    }
    return sqrt(ans/N);
}


/****************************************************************
* Functions to initalize the data                               *
****************************************************************/
template<typename TYPE>
void VectorOperationsDefault<TYPE>::zero()
{
    auto curMe = d_VectorData->begin<TYPE>();
    auto last  = d_VectorData->end<TYPE>();
    while ( curMe != last ) {
        *curMe = 0;
        ++curMe;
    }
    if ( haGhosts() ) {
        auto& ghosts = getGhosts();
        for ( size_t i = 0; i != ghosts.size(); i++ )
            ghosts[i] = 0;
    }
    *( d_VectorData->getUpdateStatusPtr() ) = VectorData::UpdateState::UNCHANGED;
}
template<typename TYPE>
void VectorOperationsDefault<TYPE>::setToScalar( double alpha )
{
    auto curMe = d_VectorData->begin<TYPE>();
    auto last  = d_VectorData->end<TYPE>();
    while ( curMe != last ) {
        *curMe = alpha;
        ++curMe;
    }
    if ( haGhosts() ) {
        auto& ghosts = getGhosts();
        for ( size_t i = 0; i != ghosts.size(); i++ )
            ghosts[i] = alpha;
    }
    *( d_VectorData->getUpdateStatusPtr() ) = VectorData::UpdateState::UNCHANGED;
}
template<typename TYPE>
void VectorOperationsDefault<TYPE>::setRandomValues()
{
    RandomVariable<double> r( 0, 1, Vector::getDefaultRNG() );
    auto curMe = d_VectorData->begin<TYPE>();
    auto last  = d_VectorData->end<TYPE>();
    while ( curMe != last ) {
        double curRand = r;
        *curMe         = curRand;
        ++curMe;
    }
    d_VectorData->dataChanged();
    d_VectorData->makeConsistent( VectorData::ScatterType::CONSISTENT_SET );
}
template<typename TYPE>
void VectorOperationsDefault<TYPE>::setRandomValues( RNG::shared_ptr rng )
{
    RandomVariable<double> r( 0, 1, rng );
    auto curMe = d_VectorData->begin<TYPE>();
    auto last  = d_VectorData->end<TYPE>();
    while ( curMe != last ) {
        *curMe = r;
        ++curMe;
    }
    d_VectorData->dataChanged();
    d_VectorData->makeConsistent( VectorData::ScatterType::CONSISTENT_SET );
}


/****************************************************************
* Basic linear algebra                                          *
****************************************************************/
template<typename TYPE>
void VectorOperationsDefault<TYPE>::scale( double alpha )
{
    auto curMe = d_VectorData->begin<TYPE>();
    auto last  = d_VectorData->end<TYPE>();
    while ( curMe != last ) {
        *curMe *= alpha;
        ++curMe;
    }
    d_VectorData->dataChanged();
}
template<typename TYPE>
void VectorOperationsDefault<TYPE>::scale( double alpha, const VectorOperations &x )
{
    AMP_ASSERT( d_VectorData->getLocalSize() == x.getVectorData()->getLocalSize() );
    auto curMe = d_VectorData->begin<TYPE>();
    auto last  = d_VectorData->end<TYPE>();
    auto curRhs = x.getVectorData()->begin<TYPE>();
    while ( curMe != last ) {
        *curMe = alpha * *curRhs;
        ++curRhs;
        ++curMe;
    }
    d_VectorData->dataChanged();
}
template<typename TYPE>
void VectorOperationsDefault<TYPE>::add( const VectorOperations &x, const VectorOperations &y )
{
    AMP_ASSERT( d_VectorData->getLocalSize() == x.getVectorData()->getLocalSize() );
    AMP_ASSERT( d_VectorData->getLocalSize() == y.getVectorData()->getLocalSize() );
    auto curMe = d_VectorData->begin<TYPE>();
    auto last  = d_VectorData->end<TYPE>();
    auto curXRhs = x.getVectorData()->begin<TYPE>();
    auto curYRhs = y.getVectorData()->begin<TYPE>();
    while ( curMe != last ) {
        *curMe = *curXRhs + *curYRhs;
        ++curXRhs;
        ++curYRhs;
        ++curMe;
    }
    d_VectorData->dataChanged();
}
template<typename TYPE>
void VectorOperationsDefault<TYPE>::subtract( const VectorOperations &x, const VectorOperations &y )
{
    AMP_ASSERT( d_VectorData->getLocalSize() == x.getVectorData()->getLocalSize() );
    AMP_ASSERT( d_VectorData->getLocalSize() == y.getVectorData()->getLocalSize() );
    auto curMe = d_VectorData->begin<TYPE>();
    auto last  = d_VectorData->end<TYPE>();
    auto curXRhs = x.getVectorData()->begin<TYPE>();
    auto curYRhs = y.getVectorData()->begin<TYPE>();
    while ( curMe != last ) {
        *curMe = *curXRhs - *curYRhs;
        ++curXRhs;
        ++curYRhs;
        ++curMe;
    }
    d_VectorData->dataChanged();
}
template<typename TYPE>
void VectorOperationsDefault<TYPE>::multiply( const VectorOperations &x, const VectorOperations &y )
{
    AMP_ASSERT( d_VectorData->getLocalSize() == x.getVectorData()->getLocalSize() );
    AMP_ASSERT( d_VectorData->getLocalSize() == y.getVectorData()->getLocalSize() );
    auto curMe = d_VectorData->begin<TYPE>();
    auto last  = d_VectorData->end<TYPE>();
    auto curXRhs = x.getVectorData()->begin<TYPE>();
    auto curYRhs = y.getVectorData()->begin<TYPE>();
    while ( curMe != last ) {
        *curMe = *curXRhs * *curYRhs;
        ++curXRhs;
        ++curYRhs;
        ++curMe;
    }
    d_VectorData->dataChanged();
}
template<typename TYPE>
void VectorOperationsDefault<TYPE>::divide( const VectorOperations &x, const VectorOperations &y )
{
    AMP_ASSERT( d_VectorData->getLocalSize() == x.getVectorData()->getLocalSize() );
    AMP_ASSERT( d_VectorData->getLocalSize() == y.getVectorData()->getLocalSize() );
    auto curMe = d_VectorData->begin<TYPE>();
    auto last  = d_VectorData->end<TYPE>();
    auto curXRhs = x.getVectorData()->begin<TYPE>();
    auto curYRhs = y.getVectorData()->begin<TYPE>();
    while ( curMe != last ) {
        *curMe = *curXRhs / *curYRhs;
        ++curXRhs;
        ++curYRhs;
        ++curMe;
    }
    d_VectorData->dataChanged();
}
template<typename TYPE>
void VectorOperationsDefault<TYPE>::reciprocal( const VectorOperations &x )
{
    AMP_ASSERT( d_VectorData->getLocalSize() == x.getVectorData()->getLocalSize() );
    auto curMe = d_VectorData->begin<TYPE>();
    auto last  = d_VectorData->end<TYPE>();
    auto curRhs = x.getVectorData()->begin<TYPE>();
    while ( curMe != last ) {
        *curMe = 1. / *curRhs;
        ++curRhs;
        ++curMe;
    }
    d_VectorData->dataChanged();
}
template<typename TYPE>
void VectorOperationsDefault<TYPE>::linearSum( double alpha_in,
                        const VectorOperations &x,
                        double beta_in,
                        const VectorOperations &y )
{
    const TYPE alpha = static_cast<TYPE>( alpha_in );
    const TYPE beta  = static_cast<TYPE>( beta_in );
    AMP_ASSERT( d_VectorData->getLocalSize() == x.getVectorData()->getLocalSize() );
    AMP_ASSERT( d_VectorData->getLocalSize() == y.getVectorData()->getLocalSize() );
    auto curMe = d_VectorData->begin<TYPE>();
    auto last  = d_VectorData->end<TYPE>();
    auto curXRhs = x.getVectorData()->begin<TYPE>();
    auto curYRhs = y.getVectorData()->begin<TYPE>();
    while ( curMe != last ) {
        *curMe = alpha * *curXRhs + beta * *curYRhs;
        ++curXRhs;
        ++curYRhs;
        ++curMe;
    }
    d_VectorData->dataChanged();
}
template<typename TYPE>
void VectorOperationsDefault<TYPE>::axpy( double alpha_in, const VectorOperations &x, const VectorOperations &y )
{
    const TYPE alpha = static_cast<TYPE>( alpha_in );
    AMP_ASSERT( d_VectorData->getLocalSize() == x.getVectorData()->getLocalSize() );
    AMP_ASSERT( d_VectorData->getLocalSize() == y.getVectorData()->getLocalSize() );
    auto curMe = d_VectorData->begin<TYPE>();
    auto last  = d_VectorData->end<TYPE>();
    auto curXRhs = x.getVectorData()->begin<TYPE>();
    auto curYRhs = y.getVectorData()->begin<TYPE>();
    while ( curMe != last ) {
        *curMe = alpha * *curXRhs + *curYRhs;
        ++curXRhs;
        ++curYRhs;
        ++curMe;
    }
    d_VectorData->dataChanged();
}
template<typename TYPE>
void VectorOperationsDefault<TYPE>::axpby( double alpha_in, double beta_in, const VectorOperations &x )
{
    const TYPE alpha = static_cast<TYPE>( alpha_in );
    const TYPE beta  = static_cast<TYPE>( beta_in );
    AMP_ASSERT( d_VectorData->getLocalSize() == x.getVectorData()->getLocalSize() );
    auto curMe = d_VectorData->begin<TYPE>();
    auto last  = d_VectorData->end<TYPE>();
    auto curXRhs = x.getVectorData()->begin<TYPE>();
    while ( curMe != last ) {
        *curMe = alpha * *curXRhs + beta * *curMe;
        ++curXRhs;
        ++curMe;
    }
    d_VectorData->dataChanged();
}
template<typename TYPE>
void VectorOperationsDefault<TYPE>::abs( const VectorOperations &x )
{
    AMP_ASSERT( d_VectorData->getLocalSize() == x.getVectorData()->getLocalSize() );
    auto curMe = d_VectorData->begin<TYPE>();
    auto last  = d_VectorData->end<TYPE>();
    auto curXRhs = x.getVectorData()->begin<TYPE>();
    while ( curMe != last ) {
        *curMe = fabs( *curXRhs );
        ++curXRhs;
        ++curMe;
    }
    d_VectorData->dataChanged();
}
template<typename TYPE>
void VectorOperationsDefault<TYPE>::addScalar( const VectorOperations &x, double alpha_in )
{
    const TYPE alpha = static_cast<TYPE>( alpha_in );
    AMP_ASSERT( d_VectorData->getLocalSize() == x.getVectorData()->getLocalSize() );
    auto curMe = d_VectorData->begin<TYPE>();
    auto last  = d_VectorData->end<TYPE>();
    auto curXRhs = x.getVectorData()->begin<TYPE>();
    while ( curMe != last ) {
        *curMe = *curXRhs + alpha;
        ++curXRhs;
        ++curMe;
    }
    d_VectorData->dataChanged();
}


} // LinearAlgebra namespace
} // AMP namespace

#endif
