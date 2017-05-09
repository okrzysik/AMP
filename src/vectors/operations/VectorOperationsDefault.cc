#include "vectors/operations/VectorOperationsDefault.h"
#include "vectors/VectorData.h"
#include "vectors/Vector.h"


namespace AMP {
namespace LinearAlgebra {


/****************************************************************
* min, max, norms, etc.                                         *
****************************************************************/
bool VectorOperationsDefault::localEquals( const VectorOperations &rhs, double tol ) const
{
    const auto& x = *d_VectorData;
    const auto& y = *rhs.getVectorData();
    if ( ( x.getGlobalSize() != y.getGlobalSize() ) || ( x.getLocalSize() != y.getLocalSize() ) )
        return false;
    bool equal = true;
    auto cur1 = x.begin();
    auto cur2 = y.begin();
    auto last = x.end();
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
double VectorOperationsDefault::localMin( void ) const
{
    size_t N_blocks = d_VectorData->numberOfDataBlocks();
    double ans      = 1e300;
    for ( size_t i = 0; i < N_blocks; i++ ) {
        size_t size        = d_VectorData->sizeOfDataBlock( i );
        const double *data = d_VectorData->getRawDataBlock<double>( i );
        for ( size_t j = 0; j < size; j++ )
            ans = std::min( data[j], ans );
    }
    return ans;
}
double VectorOperationsDefault::localMax( void ) const
{
    size_t N_blocks = d_VectorData->numberOfDataBlocks();
    double ans      = -1e300;
    for ( size_t i = 0; i < N_blocks; i++ ) {
        size_t size        = d_VectorData->sizeOfDataBlock( i );
        const double *data = d_VectorData->getRawDataBlock<double>( i );
        for ( size_t j = 0; j < size; j++ )
            ans = std::max( data[j], ans );
    }
    return ans;
}
double VectorOperationsDefault::localL1Norm( void ) const
{
    size_t N_blocks = d_VectorData->numberOfDataBlocks();
    double ans      = 0.0;
    for ( size_t i = 0; i < N_blocks; i++ ) {
        size_t size        = d_VectorData->sizeOfDataBlock( i );
        const double *data = d_VectorData->getRawDataBlock<double>( i );
        for ( size_t j = 0; j < size; j++ )
            ans += fabs( data[j] );
    }
    return ans;
}
double VectorOperationsDefault::localL2Norm( void ) const
{
    size_t N_blocks = d_VectorData->numberOfDataBlocks();
    double ans      = 0.0;
    for ( size_t i = 0; i < N_blocks; i++ ) {
        size_t size        = d_VectorData->sizeOfDataBlock( i );
        const double *data = d_VectorData->getRawDataBlock<double>( i );
        for ( size_t j = 0; j < size; j++ )
            ans += data[j] * data[j];
    }
    return sqrt( ans );
}
double VectorOperationsDefault::localMaxNorm( void ) const
{
    size_t N_blocks = d_VectorData->numberOfDataBlocks();
    double ans      = 0.0;
    for ( size_t i = 0; i < N_blocks; i++ ) {
        size_t size        = d_VectorData->sizeOfDataBlock( i );
        const double *data = d_VectorData->getRawDataBlock<double>( i );
        for ( size_t j = 0; j < size; j++ )
            ans = std::max( fabs( data[j] ), ans );
    }
    return ans;
}
double VectorOperationsDefault::localDot( const VectorOperations& x ) const
{
    AMP_ASSERT( d_VectorData->getLocalSize() == x.getVectorData()->getLocalSize() );
    auto curMe   = d_VectorData->begin();
    auto last    = d_VectorData->end();
    auto curXRhs = x.getVectorData()->begin();
    double ans             = 0.0;
    while ( curMe != last ) {
        ans += *curMe * *curXRhs;
        ++curXRhs;
        ++curMe;
    }
    return ans;
}
double VectorOperationsDefault::localMinQuotient( const VectorOperations &x ) const
{
    auto curx = x.getVectorData()->begin();
    auto endx = x.getVectorData()->end();
    auto cury = d_VectorData->begin();
    double ans = std::numeric_limits<double>::max();
    while ( curx != endx ) {
        if ( *cury != 0.0 )
            ans = std::min( ans, ( *curx ) / ( *cury ) );
        ++curx;
        ++cury;
    }
    return ans;
}
double VectorOperationsDefault::localWrmsNorm( const VectorOperations &x ) const
{
    auto curx = x.getVectorData()->begin();
    auto endx = x.getVectorData()->end();
    auto cury = d_VectorData->begin();
    double ans = 0;
    size_t N = 0;
    while ( curx != endx ) {
        ans += (*curx)*(*curx) * (*cury)*(*cury);
        ++curx;
        ++cury;
        ++N;
    }
    return sqrt(ans/N);
}
double VectorOperationsDefault::localWrmsNormMask( const VectorOperations &x,
                             const VectorOperations &mask ) const
{
    auto curx = x.getVectorData()->begin();
    auto endx = x.getVectorData()->end();
    auto cury = d_VectorData->begin();
    auto curm = mask.getVectorData()->begin();
    double ans = 0;
    size_t N = 0;
    while ( curx != endx ) {
        if ( *curm > 0.0 )
            ans += (*curx)*(*curx) * (*cury)*(*cury);
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
void VectorOperationsDefault::zero()
{
    auto curMe = d_VectorData->begin();
    auto last  = d_VectorData->end();
    while ( curMe != last ) {
        *curMe = 0;
        ++curMe;
    }
    if ( haGhosts() ) {
        auto& ghosts = getGhosts();
        for ( size_t i = 0; i != ghosts.size(); i++ )
            ghosts[i] = 0.0;
    }
    *( d_VectorData->getUpdateStatusPtr() ) = VectorData::UpdateState::UNCHANGED;
}
void VectorOperationsDefault::setToScalar( double alpha )
{
    auto curMe = d_VectorData->begin();
    auto last  = d_VectorData->end();
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
void VectorOperationsDefault::setRandomValues()
{
    RandomVariable<double> r( 0., 1., Vector::getDefaultRNG() );
    auto curMe = d_VectorData->begin();
    auto last  = d_VectorData->end();
    while ( curMe != last ) {
        double curRand = r;
        *curMe         = curRand;
        ++curMe;
    }
    d_VectorData->dataChanged();
    d_VectorData->makeConsistent( VectorData::ScatterType::CONSISTENT_SET );
}
void VectorOperationsDefault::setRandomValues( RNG::shared_ptr rng )
{
    RandomVariable<double> r( 0., 1., rng );
    auto curMe = d_VectorData->begin();
    auto last  = d_VectorData->end();
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
void VectorOperationsDefault::scale( double alpha )
{
    auto curMe = d_VectorData->begin();
    auto last  = d_VectorData->end();
    while ( curMe != last ) {
        *curMe *= alpha;
        ++curMe;
    }
    d_VectorData->dataChanged();
}
void VectorOperationsDefault::scale( double alpha, const VectorOperations &x )
{
    AMP_ASSERT( d_VectorData->getLocalSize() == x.getVectorData()->getLocalSize() );
    auto curMe = d_VectorData->begin();
    auto last  = d_VectorData->end();
    auto curRhs = x.getVectorData()->begin();
    while ( curMe != last ) {
        *curMe = alpha * *curRhs;
        ++curRhs;
        ++curMe;
    }
    d_VectorData->dataChanged();
}
void VectorOperationsDefault::add( const VectorOperations &x, const VectorOperations &y )
{
    AMP_ASSERT( d_VectorData->getLocalSize() == x.getVectorData()->getLocalSize() );
    AMP_ASSERT( d_VectorData->getLocalSize() == y.getVectorData()->getLocalSize() );
    auto curMe = d_VectorData->begin();
    auto last  = d_VectorData->end();
    auto curXRhs = x.getVectorData()->begin();
    auto curYRhs = y.getVectorData()->begin();
    while ( curMe != last ) {
        *curMe = *curXRhs + *curYRhs;
        ++curXRhs;
        ++curYRhs;
        ++curMe;
    }
    d_VectorData->dataChanged();
}
void VectorOperationsDefault::subtract( const VectorOperations &x, const VectorOperations &y )
{
    AMP_ASSERT( d_VectorData->getLocalSize() == x.getVectorData()->getLocalSize() );
    AMP_ASSERT( d_VectorData->getLocalSize() == y.getVectorData()->getLocalSize() );
    auto curMe = d_VectorData->begin();
    auto last  = d_VectorData->end();
    auto curXRhs = x.getVectorData()->begin();
    auto curYRhs = y.getVectorData()->begin();
    while ( curMe != last ) {
        *curMe = *curXRhs - *curYRhs;
        ++curXRhs;
        ++curYRhs;
        ++curMe;
    }
    d_VectorData->dataChanged();
}
void VectorOperationsDefault::multiply( const VectorOperations &x, const VectorOperations &y )
{
    AMP_ASSERT( d_VectorData->getLocalSize() == x.getVectorData()->getLocalSize() );
    AMP_ASSERT( d_VectorData->getLocalSize() == y.getVectorData()->getLocalSize() );
    auto curMe = d_VectorData->begin();
    auto last  = d_VectorData->end();
    auto curXRhs = x.getVectorData()->begin();
    auto curYRhs = y.getVectorData()->begin();
    while ( curMe != last ) {
        *curMe = *curXRhs * *curYRhs;
        ++curXRhs;
        ++curYRhs;
        ++curMe;
    }
    d_VectorData->dataChanged();
}
void VectorOperationsDefault::divide( const VectorOperations &x, const VectorOperations &y )
{
    AMP_ASSERT( d_VectorData->getLocalSize() == x.getVectorData()->getLocalSize() );
    AMP_ASSERT( d_VectorData->getLocalSize() == y.getVectorData()->getLocalSize() );
    auto curMe = d_VectorData->begin();
    auto last  = d_VectorData->end();
    auto curXRhs = x.getVectorData()->begin();
    auto curYRhs = y.getVectorData()->begin();
    while ( curMe != last ) {
        *curMe = *curXRhs / *curYRhs;
        ++curXRhs;
        ++curYRhs;
        ++curMe;
    }
    d_VectorData->dataChanged();
}
void VectorOperationsDefault::reciprocal( const VectorOperations &x )
{
    AMP_ASSERT( d_VectorData->getLocalSize() == x.getVectorData()->getLocalSize() );
    auto curMe = d_VectorData->begin();
    auto last  = d_VectorData->end();
    auto curRhs = x.getVectorData()->begin();
    while ( curMe != last ) {
        *curMe = 1. / *curRhs;
        ++curRhs;
        ++curMe;
    }
    d_VectorData->dataChanged();
}
void VectorOperationsDefault::linearSum( double alpha,
                        const VectorOperations &x,
                        double beta,
                        const VectorOperations &y )
{
    AMP_ASSERT( d_VectorData->getLocalSize() == x.getVectorData()->getLocalSize() );
    AMP_ASSERT( d_VectorData->getLocalSize() == y.getVectorData()->getLocalSize() );
    auto curMe = d_VectorData->begin();
    auto last  = d_VectorData->end();
    auto curXRhs = x.getVectorData()->begin();
    auto curYRhs = y.getVectorData()->begin();
    while ( curMe != last ) {
        *curMe = alpha * *curXRhs + beta * *curYRhs;
        ++curXRhs;
        ++curYRhs;
        ++curMe;
    }
    d_VectorData->dataChanged();
}
void VectorOperationsDefault::axpy( double alpha, const VectorOperations &x, const VectorOperations &y )
{
    AMP_ASSERT( d_VectorData->getLocalSize() == x.getVectorData()->getLocalSize() );
    AMP_ASSERT( d_VectorData->getLocalSize() == y.getVectorData()->getLocalSize() );
    auto curMe = d_VectorData->begin();
    auto last  = d_VectorData->end();
    auto curXRhs = x.getVectorData()->begin();
    auto curYRhs = y.getVectorData()->begin();
    while ( curMe != last ) {
        *curMe = alpha * *curXRhs + *curYRhs;
        ++curXRhs;
        ++curYRhs;
        ++curMe;
    }
    d_VectorData->dataChanged();
}
void VectorOperationsDefault::axpby( double alpha, double beta, const VectorOperations &x )
{
    AMP_ASSERT( d_VectorData->getLocalSize() == x.getVectorData()->getLocalSize() );
    auto curMe = d_VectorData->begin();
    auto last  = d_VectorData->end();
    auto curXRhs = x.getVectorData()->begin();
    while ( curMe != last ) {
        *curMe = alpha * *curXRhs + beta * *curMe;
        ++curXRhs;
        ++curMe;
    }
    d_VectorData->dataChanged();
}
void VectorOperationsDefault::abs( const VectorOperations &x )
{
    AMP_ASSERT( d_VectorData->getLocalSize() == x.getVectorData()->getLocalSize() );
    auto curMe = d_VectorData->begin();
    auto last  = d_VectorData->end();
    auto curXRhs = x.getVectorData()->begin();
    while ( curMe != last ) {
        *curMe = fabs( *curXRhs );
        ++curXRhs;
        ++curMe;
    }
    d_VectorData->dataChanged();
}
void VectorOperationsDefault::addScalar( const VectorOperations &x, double alpha )
{
    AMP_ASSERT( d_VectorData->getLocalSize() == x.getVectorData()->getLocalSize() );
    auto curMe = d_VectorData->begin();
    auto last  = d_VectorData->end();
    auto curXRhs = x.getVectorData()->begin();
    while ( curMe != last ) {
        *curMe = *curXRhs + alpha;
        ++curXRhs;
        ++curMe;
    }
    d_VectorData->dataChanged();
}


} // LinearAlgebra namespace
} // AMP namespace

