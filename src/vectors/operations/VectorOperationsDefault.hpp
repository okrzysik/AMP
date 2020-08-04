#ifndef included_AMP_VectorOperationsDefault_hpp
#define included_AMP_VectorOperationsDefault_hpp

#include "AMP/vectors/Vector.h"
#include "AMP/vectors/data/VectorData.h"
#include "AMP/vectors/operations/VectorOperationsDefault.h"


namespace AMP {
namespace LinearAlgebra {


extern template class VectorOperationsDefault<double>; // Suppresses implicit instantiation below --
extern template class VectorOperationsDefault<float>;  // Suppresses implicit instantiation below --


/****************************************************************
 * Constructors                                                  *
 ****************************************************************/
template<typename TYPE>
std::shared_ptr<VectorOperations> VectorOperationsDefault<TYPE>::cloneOperations() const
{
    auto ptr = std::make_shared<VectorOperationsDefault<TYPE>>();
    return ptr;
}

template<typename TYPE>
void VectorOperationsDefault<TYPE>::zero( )
{
  zero( *getVectorData() );
}

template<typename TYPE>
void VectorOperationsDefault<TYPE>::setToScalar( double alpha )
{
  setToScalar(alpha, *getVectorData() );
}

template<typename TYPE>
void VectorOperationsDefault<TYPE>::setRandomValues( void )
{
   setRandomValues( *getVectorData() );
}

template<typename TYPE>
void VectorOperationsDefault<TYPE>::setRandomValues( RNG::shared_ptr rng  )
{
  setRandomValues( rng, *getVectorData() );
}

template<typename TYPE>
void VectorOperationsDefault<TYPE>::copy( const VectorOperations &x )
{
  copy( *(x.getVectorData()), *getVectorData());
}
  
template<typename TYPE>
void VectorOperationsDefault<TYPE>::scale( double alpha, const VectorOperations &x )
{
  scale(alpha, *(x.getVectorData()), *getVectorData());
}

template<typename TYPE>
void VectorOperationsDefault<TYPE>::scale( double alpha )
{
  scale(alpha, *getVectorData());
}

template<typename TYPE>
void VectorOperationsDefault<TYPE>::add( const VectorOperations &x, const VectorOperations &y )
{
  add( *(x.getVectorData()), *(y.getVectorData()), *getVectorData() );
}

template<typename TYPE>
void VectorOperationsDefault<TYPE>::subtract( const VectorOperations &x, const VectorOperations &y )
{
  subtract( *(x.getVectorData()), *(y.getVectorData()), *getVectorData() );
}

template<typename TYPE>
void VectorOperationsDefault<TYPE>::multiply( const VectorOperations &x, const VectorOperations &y )
{
  multiply( *(x.getVectorData()), *(y.getVectorData()), *getVectorData() );
}

template<typename TYPE>
void VectorOperationsDefault<TYPE>::divide( const VectorOperations &x, const VectorOperations &y )
{
  divide( *(x.getVectorData()), *(y.getVectorData()), *getVectorData() );
}

template<typename TYPE>
void VectorOperationsDefault<TYPE>::reciprocal( const VectorOperations &x )
{
  reciprocal( *(x.getVectorData()), *getVectorData() );
}

template<typename TYPE>
void VectorOperationsDefault<TYPE>::linearSum( double alpha,
                                   const VectorOperations &x,
                                   double beta,
                                   const VectorOperations &y )
{
  linearSum( alpha,
	     *(x.getVectorData()),
	     beta,
	     *(y.getVectorData()),
	     *getVectorData() );
}

template<typename TYPE>
void VectorOperationsDefault<TYPE>::axpy( double alpha, const VectorOperations &x, const VectorOperations &y )
{
  axpy( alpha,
	*(x.getVectorData()),
	*(y.getVectorData()),
	*getVectorData() );
}

template<typename TYPE>
void VectorOperationsDefault<TYPE>::axpby( double alpha, double beta, const VectorOperations &x )
{
  axpby( alpha,
	 beta,
	 *(x.getVectorData()),
	 *getVectorData() );
}

template<typename TYPE>
void VectorOperationsDefault<TYPE>::abs( const VectorOperations &x )
{
    abs( *(x.getVectorData()), *getVectorData() );
}

template<typename TYPE>
void VectorOperationsDefault<TYPE>::addScalar( const VectorOperations &x, double alpha_in )
{
  addScalar( *(x.getVectorData()), alpha_in, *getVectorData() );
}

template<typename TYPE>
double VectorOperationsDefault<TYPE>::localMin( void ) const
{
  return localMin( *getVectorData() );
}

template<typename TYPE>
double VectorOperationsDefault<TYPE>::localMax( void ) const
{
  return localMax( *getVectorData() );
}

template<typename TYPE>
double VectorOperationsDefault<TYPE>::localL1Norm( void ) const
{
  return localL1Norm( *getVectorData() );
}

template<typename TYPE>
double VectorOperationsDefault<TYPE>::localL2Norm( void ) const
{
  return localL2Norm( *getVectorData() );
}

template<typename TYPE>
double VectorOperationsDefault<TYPE>::localMaxNorm( void ) const
{
  return localMaxNorm( *getVectorData() );
}

template<typename TYPE>
double VectorOperationsDefault<TYPE>::localDot( const VectorOperations &x ) const
{
    return localDot( *(x.getVectorData()), *getVectorData() );
}

template<typename TYPE>
bool VectorOperationsDefault<TYPE>::localEquals( const VectorOperations &rhs, double tol ) const
{
  return localEquals( *(rhs.getVectorData()), *getVectorData(), tol );
}

template<typename TYPE>
double VectorOperationsDefault<TYPE>::localMinQuotient( const VectorOperations &x ) const
{
  return localMinQuotient(*(x.getVectorData()), *getVectorData());
}

template<typename TYPE>
double VectorOperationsDefault<TYPE>::localWrmsNorm( const VectorOperations &x ) const
{
  return localWrmsNorm(*(x.getVectorData()), *getVectorData());
}

template<typename TYPE>
double VectorOperationsDefault<TYPE>::localWrmsNormMask( const VectorOperations &x,
                                                         const VectorOperations &mask ) const
{
  return localWrmsNormMask( *(x.getVectorData()), *(mask.getVectorData()), *getVectorData());
}
  
//**********************************************************************
// Static functions that operate on VectorData objects

template<typename TYPE>
void VectorOperationsDefault<TYPE>::zero( VectorData &x )
{
    auto curMe = x.begin<TYPE>();
    auto last  = x.end<TYPE>();
    while ( curMe != last ) {
        *curMe = 0;
        ++curMe;
    }
    if ( x.hasGhosts() ) {
      auto &ghosts = x.getGhosts();
        for ( size_t i = 0; i != ghosts.size(); i++ )
            ghosts[i] = 0;
    }
    // Override the status state since we set the ghost values
    *( x.getUpdateStatusPtr() ) = VectorData::UpdateState::UNCHANGED;
}

template<typename TYPE>
void VectorOperationsDefault<TYPE>::setToScalar( double alpha, VectorData &x )
{
    auto curMe = x.begin<TYPE>();
    auto last  = x.end<TYPE>();
    while ( curMe != last ) {
        *curMe = alpha;
        ++curMe;
    }
    if ( x.hasGhosts() ) {
      auto &ghosts = x.getGhosts();
        for ( size_t i = 0; i != ghosts.size(); i++ )
            ghosts[i] = alpha;
    }
    // Override the status state since we set the ghost values
    *( x.getUpdateStatusPtr() ) = VectorData::UpdateState::UNCHANGED;
}

template<typename TYPE>
void VectorOperationsDefault<TYPE>::setRandomValues( VectorData &x )
{
    RandomVariable<double> r( 0, 1, Vector::getDefaultRNG() );
    auto curMe = x.begin<TYPE>();
    auto last  = x.end<TYPE>();
    while ( curMe != last ) {
        double curRand = r;
        *curMe         = curRand;
        ++curMe;
    }
    // Call makeConsistent to leave the vector in a consistent state
    x.makeConsistent( VectorData::ScatterType::CONSISTENT_SET );
}

template<typename TYPE>
void VectorOperationsDefault<TYPE>::setRandomValues( RNG::shared_ptr rng, VectorData &x )
{
    RandomVariable<double> r( 0, 1, rng );
    auto curMe = x.begin<TYPE>();
    auto last  = x.end<TYPE>();
    while ( curMe != last ) {
        double curRand = r;
        *curMe         = curRand;
        ++curMe;
    }
    // Call makeConsistent to leave the vector in a consistent state
    x.makeConsistent( VectorData::ScatterType::CONSISTENT_SET );
}

template<typename TYPE>
void VectorOperationsDefault<TYPE>::copy( const VectorData &x, VectorData &y )
{
    AMP_ASSERT( y.getLocalSize() == x.getLocalSize() );
    std::copy( x.begin<TYPE>(), x.end<TYPE>(), y.begin<TYPE>() );
    y.copyGhostValues( x );
}

template<typename TYPE>
void VectorOperationsDefault<TYPE>::scale( double alpha, VectorData &x )
{
    auto curMe = x.begin<TYPE>();
    auto last  = x.end<TYPE>();
    while ( curMe != last ) {
        *curMe *= alpha;
        ++curMe;
    }
}

template<typename TYPE>
void VectorOperationsDefault<TYPE>::scale( double alpha, const VectorData &x, VectorData &y )
{
    AMP_ASSERT( y.getLocalSize() == x.getLocalSize() );
    auto curMe  = y.begin<TYPE>();
    auto last   = y.end<TYPE>();
    auto curRhs = x.begin<TYPE>();
    while ( curMe != last ) {
        *curMe = alpha * *curRhs;
        ++curRhs;
        ++curMe;
    }
}

template<typename TYPE>
void VectorOperationsDefault<TYPE>::add( const VectorData &x, const VectorData &y, VectorData &z )
{
    AMP_ASSERT( z.getLocalSize() == x.getLocalSize() );
    AMP_ASSERT( z.getLocalSize() == y.getLocalSize() );
    auto curMe   = z.begin<TYPE>();
    auto last    = z.end<TYPE>();
    auto curXRhs = x.begin<TYPE>();
    auto curYRhs = y.begin<TYPE>();
    while ( curMe != last ) {
        *curMe = *curXRhs + *curYRhs;
        ++curXRhs;
        ++curYRhs;
        ++curMe;
    }
}

template<typename TYPE>
void VectorOperationsDefault<TYPE>::subtract( const VectorData &x, const VectorData &y, VectorData &z  )
{
    AMP_ASSERT( z.getLocalSize() == x.getLocalSize() );
    AMP_ASSERT( z.getLocalSize() == y.getLocalSize() );
    auto curMe   = z.begin<TYPE>();
    auto last    = z.end<TYPE>();
    auto curXRhs = x.begin<TYPE>();
    auto curYRhs = y.begin<TYPE>();
    while ( curMe != last ) {
        *curMe = *curXRhs - *curYRhs;
        ++curXRhs;
        ++curYRhs;
        ++curMe;
    }
}

template<typename TYPE>
void VectorOperationsDefault<TYPE>::multiply( const VectorData &x, const VectorData &y, VectorData &z )
{
    AMP_ASSERT( z.getLocalSize() == x.getLocalSize() );
    AMP_ASSERT( z.getLocalSize() == y.getLocalSize() );
    auto curMe   = z.begin<TYPE>();
    auto last    = z.end<TYPE>();
    auto curXRhs = x.begin<TYPE>();
    auto curYRhs = y.begin<TYPE>();
    while ( curMe != last ) {
        *curMe = *curXRhs * *curYRhs;
        ++curXRhs;
        ++curYRhs;
        ++curMe;
    }
}

template<typename TYPE>
void VectorOperationsDefault<TYPE>::divide( const VectorData &x, const VectorData &y, VectorData &z )
{
    AMP_ASSERT( z.getLocalSize() == x.getLocalSize() );
    AMP_ASSERT( z.getLocalSize() == y.getLocalSize() );
    auto curMe   = z.begin<TYPE>();
    auto last    = z.end<TYPE>();
    auto curXRhs = x.begin<TYPE>();
    auto curYRhs = y.begin<TYPE>();
    while ( curMe != last ) {
        *curMe = *curXRhs / *curYRhs;
        ++curXRhs;
        ++curYRhs;
        ++curMe;
    }
}


template<typename TYPE>
void VectorOperationsDefault<TYPE>::reciprocal( const VectorData &x, VectorData &y )
{
    AMP_ASSERT( y.getLocalSize() == x.getLocalSize() );
    auto curMe  = y.begin<TYPE>();
    auto last   = y.end<TYPE>();
    auto curRhs = x.begin<TYPE>();
    while ( curMe != last ) {
        *curMe = 1.0 / *curRhs;
        ++curRhs;
        ++curMe;
    }
}


template<typename TYPE>
void VectorOperationsDefault<TYPE>::linearSum( double alpha_in,
                                   const VectorData &x,
                                   double beta_in,
                                   const VectorData &y,
				   VectorData &z)
{
    const TYPE alpha = static_cast<TYPE>( alpha_in );
    const TYPE beta  = static_cast<TYPE>( beta_in );
    AMP_ASSERT( z.getLocalSize() == x.getLocalSize() );
    AMP_ASSERT( z.getLocalSize() == y.getLocalSize() );
    auto curMe   = z.begin<TYPE>();
    auto last    = z.end<TYPE>();
    auto curXRhs = x.begin<TYPE>();
    auto curYRhs = y.begin<TYPE>();
    while ( curMe != last ) {
        *curMe = alpha * *curXRhs + beta * *curYRhs;
        ++curXRhs;
        ++curYRhs;
        ++curMe;
    }
}

template<typename TYPE>
void VectorOperationsDefault<TYPE>::axpy( double alpha_in, const VectorData &x, const VectorData &y, VectorData &z )
{
    const TYPE alpha = static_cast<TYPE>( alpha_in );
    AMP_ASSERT( z.getLocalSize() == x.getLocalSize() );
    AMP_ASSERT( z.getLocalSize() == y.getLocalSize() );
    auto curMe   = z.begin<TYPE>();
    auto last    = z.end<TYPE>();
    auto curXRhs = x.begin<TYPE>();
    auto curYRhs = y.begin<TYPE>();
    while ( curMe != last ) {
        *curMe = alpha * *curXRhs + *curYRhs;
        ++curXRhs;
        ++curYRhs;
        ++curMe;
    }
}

template<typename TYPE>
void VectorOperationsDefault<TYPE>::axpby( double alpha_in, double beta_in, const VectorData &x, VectorData &z )
{
    const TYPE alpha = static_cast<TYPE>( alpha_in );
    const TYPE beta  = static_cast<TYPE>( beta_in );
    AMP_ASSERT( z.getLocalSize() == x.getLocalSize() );
    auto curMe   = z.begin<TYPE>();
    auto last    = z.end<TYPE>();
    auto curXRhs = x.begin<TYPE>();
    while ( curMe != last ) {
        *curMe = alpha * *curXRhs + beta * *curMe;
        ++curXRhs;
        ++curMe;
    }
}

template<typename TYPE>
void VectorOperationsDefault<TYPE>::abs( const VectorData &x, VectorData &y )
{
    AMP_ASSERT( y.getLocalSize() == x.getLocalSize() );
    auto curMe  = y.begin<TYPE>();
    auto last   = y.end<TYPE>();
    auto curRhs = x.begin<TYPE>();
    while ( curMe != last ) {
        *curMe = fabs( *curRhs );
        ++curRhs;
        ++curMe;
    }
}

template<typename TYPE>
void VectorOperationsDefault<TYPE>::addScalar( const VectorData &x, double alpha_in, VectorData &y )
{
    const TYPE alpha = static_cast<TYPE>( alpha_in );
    AMP_ASSERT( y.getLocalSize() == x.getLocalSize() );
    auto curMe   = y.begin<TYPE>();
    auto last    = y.end<TYPE>();
    auto curXRhs = x.begin<TYPE>();
    while ( curMe != last ) {
        *curMe = *curXRhs + alpha;
        ++curXRhs;
        ++curMe;
    }
}

template<typename TYPE>
double VectorOperationsDefault<TYPE>::localMin( const VectorData &x ) 
{
    size_t N_blocks = x.numberOfDataBlocks();
    TYPE ans        = std::numeric_limits<TYPE>::max();
    for ( size_t i = 0; i < N_blocks; i++ ) {
        size_t size      = x.sizeOfDataBlock( i );
        const TYPE *data = x.getRawDataBlock<TYPE>( i );
        for ( size_t j = 0; j < size; j++ )
            ans = std::min( data[j], ans );
    }
    return static_cast<double>( ans );
}

template<typename TYPE>
double VectorOperationsDefault<TYPE>::localMax( const VectorData &x ) 
{
    size_t N_blocks = x.numberOfDataBlocks();
    TYPE ans        = std::numeric_limits<TYPE>::lowest();
    for ( size_t i = 0; i < N_blocks; i++ ) {
        size_t size      = x.sizeOfDataBlock( i );
        const TYPE *data = x.getRawDataBlock<TYPE>( i );
        for ( size_t j = 0; j < size; j++ )
            ans = std::max( data[j], ans );
    }
    return static_cast<double>( ans );
}

template<typename TYPE>
double VectorOperationsDefault<TYPE>::localL1Norm( const VectorData &x ) 
{
    size_t N_blocks = x.numberOfDataBlocks();
    double ans      = 0;
    for ( size_t i = 0; i < N_blocks; i++ ) {
        size_t size      = x.sizeOfDataBlock( i );
        const TYPE *data = x.getRawDataBlock<TYPE>( i );
        for ( size_t j = 0; j < size; j++ )
            ans += static_cast<double>( std::abs( data[j] ) );
    }
    return ans;
}

template<typename TYPE>
double VectorOperationsDefault<TYPE>::localL2Norm( const VectorData &x ) 
{
    size_t N_blocks = x.numberOfDataBlocks();
    double ans      = 0;
    for ( size_t i = 0; i < N_blocks; i++ ) {
        size_t size      = x.sizeOfDataBlock( i );
        const TYPE *data = x.getRawDataBlock<TYPE>( i );
        for ( size_t j = 0; j < size; j++ ) {
            double tmp = static_cast<double>( data[j] );
            ans += tmp * tmp;
        }
    }
    return sqrt( ans );
}

template<typename TYPE>
double VectorOperationsDefault<TYPE>::localMaxNorm( const VectorData &x ) 
{
    size_t N_blocks = x.numberOfDataBlocks();
    TYPE ans        = 0;
    for ( size_t i = 0; i < N_blocks; i++ ) {
        size_t size      = x.sizeOfDataBlock( i );
        const TYPE *data = x.getRawDataBlock<TYPE>( i );
        for ( size_t j = 0; j < size; j++ )
            ans = std::max( std::abs( data[j] ), ans );
    }
    return static_cast<double>( ans );
}

template<typename TYPE>
double VectorOperationsDefault<TYPE>::localDot( const VectorData &x, const VectorData &y )
{
    AMP_ASSERT( y.getLocalSize() == x.getLocalSize() );
    auto curMe   = y.constBegin<TYPE>();
    auto last    = y.constEnd<TYPE>();
    auto curXRhs = x.constBegin<TYPE>();
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
double VectorOperationsDefault<TYPE>::localMinQuotient( const VectorData &x, const VectorData &y )
{
    auto curx  = x.constBegin<TYPE>();
    auto endx  = x.constEnd<TYPE>();
    auto cury  = y.constBegin<TYPE>();
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
double VectorOperationsDefault<TYPE>::localWrmsNorm( const VectorData &x, const VectorData &y )
{
    auto curx  = x.constBegin<TYPE>();
    auto endx  = x.constEnd<TYPE>();
    auto cury  = y.constBegin<TYPE>();
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
double VectorOperationsDefault<TYPE>::localWrmsNormMask( const VectorData &x, const VectorData &mask, const VectorData &y )
{
    auto curx  = x.constBegin<TYPE>();
    auto endx  = x.constEnd<TYPE>();
    auto cury  = y.constBegin<TYPE>();
    auto curm  = mask.constBegin<TYPE>();
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

template<typename TYPE>
bool VectorOperationsDefault<TYPE>::localEquals( const VectorData &x, const VectorData &y, double tol )
{
    if ( ( x.getGlobalSize() != y.getGlobalSize() ) || ( x.getLocalSize() != y.getLocalSize() ) )
        return false;
    bool equal = true;
    auto cur1  = x.constBegin<TYPE>();
    auto cur2  = y.constBegin<TYPE>();
    auto last  = x.constEnd<TYPE>();
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


} // namespace LinearAlgebra
} // namespace AMP

#endif
