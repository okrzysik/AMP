#include "vectors/VectorOperationsDefault.h"
#include "vectors/VectorData.h"


namespace AMP {
namespace LinearAlgebra {


/****************************************************************
* min, max, norms, etc.                                         *
****************************************************************/
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


} // LinearAlgebra namespace
} // AMP namespace

