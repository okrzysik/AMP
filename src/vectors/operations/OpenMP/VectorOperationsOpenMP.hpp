#ifndef included_AMP_VectorOperationsOpenMP_hpp
#define included_AMP_VectorOperationsOpenMP_hpp

#include "AMP/vectors/Vector.h"
#include "AMP/vectors/data/VectorData.h"
#include "AMP/vectors/operations/OpenMP/VectorOperationsOpenMP.h"


namespace AMP {
namespace LinearAlgebra {


extern template class VectorOperationsOpenMP<double>; // Suppresses implicit instantiation below --
extern template class VectorOperationsOpenMP<float>;  // Suppresses implicit instantiation below --


/****************************************************************
 * Constructors                                                  *
 ****************************************************************/
template<typename TYPE>
AMP::shared_ptr<VectorOperations> VectorOperationsOpenMP<TYPE>::cloneOperations() const
{
    auto ptr = AMP::make_shared<VectorOperationsOpenMP<TYPE>>();
    return ptr;
}


/****************************************************************
 * min, max, norms, etc.                                         *
 ****************************************************************/
template<typename TYPE>
bool VectorOperationsOpenMP<TYPE>::localEquals( const VectorOperations &rhs, double tol ) const
{
    const auto &x = *d_VectorData;
    const auto &y = *rhs.getVectorData();
    if ( ( x.getGlobalSize() != y.getGlobalSize() ) || ( x.getLocalSize() != y.getLocalSize() ) )
        return false;
    bool equal      = true;
    const auto last = x.constEnd<TYPE>();
    //#pragma omp parallel for reduction(&&:equal)
    for ( auto it1 = x.constBegin<TYPE>(), it2 = y.constBegin<TYPE>(); it1 < last; ++it1, ++it2 )
        equal = equal && fabs( *it1 - *it2 ) <= tol;
    return equal;
}
template<typename TYPE>
double VectorOperationsOpenMP<TYPE>::localMin( void ) const
{
    TYPE min        = std::numeric_limits<TYPE>::max();
    const auto last = d_VectorData->end<TYPE>();
#pragma omp parallel for reduction( min : min )
    for ( auto it = d_VectorData->begin<TYPE>(); it < last; ++it )
        min = std::min( *it, min );
    return static_cast<double>( min );
}
template<typename TYPE>
double VectorOperationsOpenMP<TYPE>::localMax( void ) const
{
    TYPE max        = -std::numeric_limits<TYPE>::max();
    const auto last = d_VectorData->end<TYPE>();
#pragma omp parallel for reduction( max : max )
    for ( auto it = d_VectorData->begin<TYPE>(); it < last; ++it )
        max = std::max( *it, max );
    return static_cast<double>( max );
}
template<typename TYPE>
double VectorOperationsOpenMP<TYPE>::localL1Norm( void ) const
{
    TYPE sum        = 0;
    const auto last = d_VectorData->end<TYPE>();
#pragma omp parallel for reduction( + : sum )
    for ( auto it = d_VectorData->begin<TYPE>(); it < last; ++it )
        sum += std::abs( *it );
    return static_cast<double>( sum );
}
template<typename TYPE>
double VectorOperationsOpenMP<TYPE>::localL2Norm( void ) const
{
    TYPE sum        = 0;
    const auto last = d_VectorData->end<TYPE>();
#pragma omp parallel for reduction( + : sum )
    for ( auto it = d_VectorData->begin<TYPE>(); it < last; ++it )
        sum += ( *it ) * ( *it );
    return sqrt( sum );
}
template<typename TYPE>
double VectorOperationsOpenMP<TYPE>::localMaxNorm( void ) const
{
    TYPE max        = -std::numeric_limits<TYPE>::max();
    const auto last = d_VectorData->end<TYPE>();
#pragma omp parallel for reduction( max : max )
    for ( auto it = d_VectorData->begin<TYPE>(); it < last; ++it )
        max = std::max( std::abs( *it ), max );
    return static_cast<double>( max );
}
template<typename TYPE>
double VectorOperationsOpenMP<TYPE>::localDot( const VectorOperations &x ) const
{
    AMP_ASSERT( d_VectorData->getLocalSize() == x.getVectorData()->getLocalSize() );
    auto curMe   = d_VectorData->constBegin<TYPE>();
    auto last    = d_VectorData->constEnd<TYPE>();
    auto curXRhs = x.getVectorData()->constBegin<TYPE>();
    double ans   = 0;
    while ( curMe != last ) {
        double v1 = static_cast<double>( *curMe );
        double v2 = static_cast<double>( *curXRhs );
        ans += v1 * v2;
        ++curXRhs;
        ++curMe;
    }
    return ans;
}
template<typename TYPE>
double VectorOperationsOpenMP<TYPE>::localMinQuotient( const VectorOperations &x ) const
{
    auto curx  = x.getVectorData()->constBegin<TYPE>();
    auto endx  = x.getVectorData()->constEnd<TYPE>();
    auto cury  = d_VectorData->constBegin<TYPE>();
    double ans = std::numeric_limits<double>::max();
    while ( curx != endx ) {
        if ( *cury != 0 ) {
            double v1 = static_cast<double>( *curx );
            double v2 = static_cast<double>( *cury );
            ans       = std::min( ans, v1 / v2 );
        }
        ++curx;
        ++cury;
    }
    return ans;
}
template<typename TYPE>
double VectorOperationsOpenMP<TYPE>::localWrmsNorm( const VectorOperations &x ) const
{
    auto curx  = x.getVectorData()->constBegin<TYPE>();
    auto endx  = x.getVectorData()->constEnd<TYPE>();
    auto cury  = d_VectorData->constBegin<TYPE>();
    double ans = 0;
    size_t N   = 0;
    while ( curx != endx ) {
        double v1 = static_cast<double>( *curx );
        double v2 = static_cast<double>( *cury );
        ans += v1 * v1 * v2 * v2;
        ++curx;
        ++cury;
        ++N;
    }
    return sqrt( ans / N );
}
template<typename TYPE>
double VectorOperationsOpenMP<TYPE>::localWrmsNormMask( const VectorOperations &x,
                                                        const VectorOperations &mask ) const
{
    auto curx  = x.getVectorData()->constBegin<TYPE>();
    auto endx  = x.getVectorData()->constEnd<TYPE>();
    auto cury  = d_VectorData->constBegin<TYPE>();
    auto curm  = mask.getVectorData()->constBegin<TYPE>();
    double ans = 0;
    size_t N   = 0;
    while ( curx != endx ) {
        if ( *curm > 0 ) {
            double v1 = static_cast<double>( *curx );
            double v2 = static_cast<double>( *cury );
            ans += v1 * v1 * v2 * v2;
        }
        ++curx;
        ++cury;
        ++curm;
        ++N;
    }
    return sqrt( ans / N );
}


/****************************************************************
 * Functions to initalize the data                               *
 ****************************************************************/
template<typename TYPE>
void VectorOperationsOpenMP<TYPE>::zero()
{
    const auto last = d_VectorData->end<TYPE>();
#pragma omp parallel for
    for ( auto it = d_VectorData->begin<TYPE>(); it < last; ++it )
        *it = 0;
    if ( hasGhosts() ) {
        auto &ghosts = getGhosts();
#pragma omp parallel for
        for ( size_t i = 0; i < ghosts.size(); i++ )
            ghosts[i] = 0;
    }
    // Override the status state since we set the ghost values
    *( d_VectorData->getUpdateStatusPtr() ) = VectorData::UpdateState::UNCHANGED;
}
template<typename TYPE>
void VectorOperationsOpenMP<TYPE>::setToScalar( double alpha )
{
    const auto last = d_VectorData->end<TYPE>();
#pragma omp parallel for
    for ( auto it = d_VectorData->begin<TYPE>(); it < last; ++it )
        *it = alpha;
    if ( hasGhosts() ) {
        auto &ghosts = getGhosts();
#pragma omp parallel for
        for ( size_t i = 0; i < ghosts.size(); i++ )
            ghosts[i] = alpha;
    }
    // Override the status state since we set the ghost values
    *( d_VectorData->getUpdateStatusPtr() ) = VectorData::UpdateState::UNCHANGED;
}
template<typename TYPE>
void VectorOperationsOpenMP<TYPE>::setRandomValues()
{
    RandomVariable<double> r( 0, 1, Vector::getDefaultRNG() );
    auto curMe = d_VectorData->begin<TYPE>();
    auto last  = d_VectorData->end<TYPE>();
    while ( curMe != last ) {
        double curRand = r;
        *curMe         = curRand;
        ++curMe;
    }
    // Call makeConsistent to leave the vector in a consistent state
    d_VectorData->makeConsistent( VectorData::ScatterType::CONSISTENT_SET );
}
template<typename TYPE>
void VectorOperationsOpenMP<TYPE>::setRandomValues( RNG::shared_ptr rng )
{
    RandomVariable<double> r( 0, 1, rng );
    auto curMe = d_VectorData->begin<TYPE>();
    auto last  = d_VectorData->end<TYPE>();
    while ( curMe != last ) {
        *curMe = r;
        ++curMe;
    }
    // Call makeConsistent to leave the vector in a consistent state
    d_VectorData->makeConsistent( VectorData::ScatterType::CONSISTENT_SET );
}


/****************************************************************
 * Basic linear algebra                                          *
 ****************************************************************/
template<typename TYPE>
void VectorOperationsOpenMP<TYPE>::copy( const VectorOperations &x )
{
    auto x_data = x.getVectorData();
    auto y_data = d_VectorData;
    AMP_ASSERT( x_data->getLocalSize() == y_data->getLocalSize() );
    std::copy( x_data->begin<TYPE>(), x_data->end<TYPE>(), y_data->template begin<TYPE>() );
    y_data->copyGhostValues( *x_data );
}
template<typename TYPE>
void VectorOperationsOpenMP<TYPE>::scale( double alpha )
{
    const auto last = d_VectorData->end<TYPE>();
#pragma omp parallel for
    for ( auto it = d_VectorData->begin<TYPE>(); it < last; ++it )
        *it *= alpha;
}
template<typename TYPE>
void VectorOperationsOpenMP<TYPE>::scale( double alpha, const VectorOperations &x )
{
    AMP_ASSERT( d_VectorData->getLocalSize() == x.getVectorData()->getLocalSize() );
    auto curMe  = d_VectorData->begin<TYPE>();
    auto last   = d_VectorData->end<TYPE>();
    auto curRhs = x.getVectorData()->begin<TYPE>();
    while ( curMe != last ) {
        *curMe = alpha * *curRhs;
        ++curRhs;
        ++curMe;
    }
}
template<typename TYPE>
void VectorOperationsOpenMP<TYPE>::add( const VectorOperations &x, const VectorOperations &y )
{
    AMP_ASSERT( d_VectorData->getLocalSize() == x.getVectorData()->getLocalSize() );
    AMP_ASSERT( d_VectorData->getLocalSize() == y.getVectorData()->getLocalSize() );
    auto curMe   = d_VectorData->begin<TYPE>();
    auto last    = d_VectorData->end<TYPE>();
    auto curXRhs = x.getVectorData()->begin<TYPE>();
    auto curYRhs = y.getVectorData()->begin<TYPE>();
    while ( curMe != last ) {
        *curMe = *curXRhs + *curYRhs;
        ++curXRhs;
        ++curYRhs;
        ++curMe;
    }
}
template<typename TYPE>
void VectorOperationsOpenMP<TYPE>::subtract( const VectorOperations &x, const VectorOperations &y )
{
    AMP_ASSERT( d_VectorData->getLocalSize() == x.getVectorData()->getLocalSize() );
    AMP_ASSERT( d_VectorData->getLocalSize() == y.getVectorData()->getLocalSize() );
    auto curMe   = d_VectorData->begin<TYPE>();
    auto last    = d_VectorData->end<TYPE>();
    auto curXRhs = x.getVectorData()->begin<TYPE>();
    auto curYRhs = y.getVectorData()->begin<TYPE>();
    while ( curMe != last ) {
        *curMe = *curXRhs - *curYRhs;
        ++curXRhs;
        ++curYRhs;
        ++curMe;
    }
}
template<typename TYPE>
void VectorOperationsOpenMP<TYPE>::multiply( const VectorOperations &x, const VectorOperations &y )
{
    AMP_ASSERT( d_VectorData->getLocalSize() == x.getVectorData()->getLocalSize() );
    AMP_ASSERT( d_VectorData->getLocalSize() == y.getVectorData()->getLocalSize() );
    auto curMe   = d_VectorData->begin<TYPE>();
    auto last    = d_VectorData->end<TYPE>();
    auto curXRhs = x.getVectorData()->begin<TYPE>();
    auto curYRhs = y.getVectorData()->begin<TYPE>();
    while ( curMe != last ) {
        *curMe = *curXRhs * *curYRhs;
        ++curXRhs;
        ++curYRhs;
        ++curMe;
    }
}
template<typename TYPE>
void VectorOperationsOpenMP<TYPE>::divide( const VectorOperations &x, const VectorOperations &y )
{
    AMP_ASSERT( d_VectorData->getLocalSize() == x.getVectorData()->getLocalSize() );
    AMP_ASSERT( d_VectorData->getLocalSize() == y.getVectorData()->getLocalSize() );
    auto curMe   = d_VectorData->begin<TYPE>();
    auto last    = d_VectorData->end<TYPE>();
    auto curXRhs = x.getVectorData()->begin<TYPE>();
    auto curYRhs = y.getVectorData()->begin<TYPE>();
    while ( curMe != last ) {
        *curMe = *curXRhs / *curYRhs;
        ++curXRhs;
        ++curYRhs;
        ++curMe;
    }
}
template<typename TYPE>
void VectorOperationsOpenMP<TYPE>::reciprocal( const VectorOperations &x )
{
    AMP_ASSERT( d_VectorData->getLocalSize() == x.getVectorData()->getLocalSize() );
    auto curMe  = d_VectorData->begin<TYPE>();
    auto last   = d_VectorData->end<TYPE>();
    auto curRhs = x.getVectorData()->begin<TYPE>();
    while ( curMe != last ) {
        *curMe = 1. / *curRhs;
        ++curRhs;
        ++curMe;
    }
}
template<typename TYPE>
void VectorOperationsOpenMP<TYPE>::linearSum( double alpha_in,
                                              const VectorOperations &x,
                                              double beta_in,
                                              const VectorOperations &y )
{
    const TYPE alpha = static_cast<TYPE>( alpha_in );
    const TYPE beta  = static_cast<TYPE>( beta_in );
    AMP_ASSERT( d_VectorData->getLocalSize() == x.getVectorData()->getLocalSize() );
    AMP_ASSERT( d_VectorData->getLocalSize() == y.getVectorData()->getLocalSize() );
    auto curMe   = d_VectorData->begin<TYPE>();
    auto last    = d_VectorData->end<TYPE>();
    auto curXRhs = x.getVectorData()->begin<TYPE>();
    auto curYRhs = y.getVectorData()->begin<TYPE>();
    while ( curMe != last ) {
        *curMe = alpha * *curXRhs + beta * *curYRhs;
        ++curXRhs;
        ++curYRhs;
        ++curMe;
    }
}
template<typename TYPE>
void VectorOperationsOpenMP<TYPE>::axpy( double alpha_in,
                                         const VectorOperations &x,
                                         const VectorOperations &y )
{
    const TYPE alpha = static_cast<TYPE>( alpha_in );
    AMP_ASSERT( d_VectorData->getLocalSize() == x.getVectorData()->getLocalSize() );
    AMP_ASSERT( d_VectorData->getLocalSize() == y.getVectorData()->getLocalSize() );
    auto curMe   = d_VectorData->begin<TYPE>();
    auto last    = d_VectorData->end<TYPE>();
    auto curXRhs = x.getVectorData()->begin<TYPE>();
    auto curYRhs = y.getVectorData()->begin<TYPE>();
    while ( curMe != last ) {
        *curMe = alpha * *curXRhs + *curYRhs;
        ++curXRhs;
        ++curYRhs;
        ++curMe;
    }
}
template<typename TYPE>
void VectorOperationsOpenMP<TYPE>::axpby( double alpha_in,
                                          double beta_in,
                                          const VectorOperations &x )
{
    const TYPE alpha = static_cast<TYPE>( alpha_in );
    const TYPE beta  = static_cast<TYPE>( beta_in );
    AMP_ASSERT( d_VectorData->getLocalSize() == x.getVectorData()->getLocalSize() );
    auto curMe   = d_VectorData->begin<TYPE>();
    auto last    = d_VectorData->end<TYPE>();
    auto curXRhs = x.getVectorData()->begin<TYPE>();
    while ( curMe != last ) {
        *curMe = alpha * *curXRhs + beta * *curMe;
        ++curXRhs;
        ++curMe;
    }
}
template<typename TYPE>
void VectorOperationsOpenMP<TYPE>::abs( const VectorOperations &x )
{
    AMP_ASSERT( d_VectorData->getLocalSize() == x.getVectorData()->getLocalSize() );
    auto curMe   = d_VectorData->begin<TYPE>();
    auto last    = d_VectorData->end<TYPE>();
    auto curXRhs = x.getVectorData()->begin<TYPE>();
    while ( curMe != last ) {
        *curMe = fabs( *curXRhs );
        ++curXRhs;
        ++curMe;
    }
}
template<typename TYPE>
void VectorOperationsOpenMP<TYPE>::addScalar( const VectorOperations &x, double alpha_in )
{
    const TYPE alpha = static_cast<TYPE>( alpha_in );
    AMP_ASSERT( d_VectorData->getLocalSize() == x.getVectorData()->getLocalSize() );
    auto curMe   = d_VectorData->begin<TYPE>();
    auto last    = d_VectorData->end<TYPE>();
    auto curXRhs = x.getVectorData()->begin<TYPE>();
    while ( curMe != last ) {
        *curMe = *curXRhs + alpha;
        ++curXRhs;
        ++curMe;
    }
}


} // namespace LinearAlgebra
} // namespace AMP

#endif
