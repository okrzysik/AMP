#include "AMP/vectors/operations/MultiVectorOperations.h"
#include "AMP/vectors/data/MultiVectorData.h"


namespace AMP {
namespace LinearAlgebra {


/****************************************************************
 * Constructors                                                  *
 ****************************************************************/
std::shared_ptr<VectorOperations> MultiVectorOperations::cloneOperations() const
{
    auto ptr = std::make_shared<MultiVectorOperations>();
    return ptr;
}

//**********************************************************************
// Static functions that operate on VectorData objects

VectorData *MultiVectorOperations::getVectorDataComponent( VectorData &x, size_t i )
{
    auto x2 = dynamic_cast<MultiVectorData *>( &x );
    AMP_ASSERT( x2 && ( i < x2->numberOfComponents() ) );
    return x2->getVectorData( i );
}
const VectorData *MultiVectorOperations::getVectorDataComponent( const VectorData &x, size_t i )
{
    auto x2 = dynamic_cast<const MultiVectorData *>( &x );
    AMP_ASSERT( x2 && ( i < x2->numberOfComponents() ) );
    return x2->getVectorData( i );
}

const MultiVectorData *MultiVectorOperations::getMultiVectorData( const VectorData &x )
{
    return dynamic_cast<const MultiVectorData *>( &x );
}

MultiVectorData *MultiVectorOperations::getMultiVectorData( VectorData &x )
{
    return dynamic_cast<MultiVectorData *>( &x );
}

void MultiVectorOperations::zero( VectorData &x )
{
    auto mData = getMultiVectorData( x );

    for ( size_t i = 0; i != mData->numberOfComponents(); ++i ) {

        d_operations[i]->zero( *getVectorDataComponent( x, i ) );
    }
}

void MultiVectorOperations::setToScalar( double alpha, VectorData &x )
{
    for ( size_t i = 0; i != d_operations.size(); i++ ) {
        d_operations[i]->setToScalar( alpha, *getVectorDataComponent( x, i ) );
    }
}

void MultiVectorOperations::setRandomValues( VectorData &x )
{
    for ( size_t i = 0; i != d_operations.size(); i++ ) {
        d_operations[i]->setRandomValues( *getVectorDataComponent( x, i ) );
    }
}

void MultiVectorOperations::setRandomValues( RNG::shared_ptr rng, VectorData &x )
{
    for ( size_t i = 0; i != d_operations.size(); i++ ) {
        d_operations[i]->setRandomValues( rng, *getVectorDataComponent( x, i ) );
    }
}

void MultiVectorOperations::copy( const VectorData &x, VectorData &y )
{

    auto xc = getMultiVectorData( x );
    auto yc = getMultiVectorData( y );

    if ( xc && yc ) {
        // Both this and x are multivectors
        for ( size_t i = 0; i != d_operations.size(); i++ )
            d_operations[i]->copy( *getVectorDataComponent( x, i ),
                                   *getVectorDataComponent( y, i ) );
    } else {
        // x is not a multivector, try to call a default implementation
        AMP_ASSERT( x.getLocalSize() == y.getLocalSize() );
        if ( x.isType<double>() && y.isType<double>() ) {
            std::copy( x.begin<double>(), x.end<double>(), y.begin<double>() );
        } else if ( x.isType<float>() && y.isType<float>() ) {
            std::copy( x.begin<float>(), x.end<float>(), y.begin<float>() );
        } else {
            AMP_ERROR( "Unable to discern data types" );
        }
    }
}

void MultiVectorOperations::scale( double alpha, VectorData &x )
{
    AMP_ASSERT( getMultiVectorData( x ) );
    if ( d_operations.empty() ) {
        return;
    }
    for ( size_t i = 0; i != d_operations.size(); i++ )
        d_operations[i]->scale( alpha, *getVectorDataComponent( x, i ) );
}

void MultiVectorOperations::scale( double alpha, const VectorData &x, VectorData &y )
{
    if ( d_operations.empty() ) {
        return;
    }
    auto x2 = getMultiVectorData( x );
    auto y2 = getMultiVectorData( y );
    if ( x2 && y2 ) {
        AMP_ASSERT( d_operations.size() == x2->numberOfComponents() );
        for ( size_t i = 0; i != d_operations.size(); i++ )
            d_operations[i]->scale(
                alpha, *getVectorDataComponent( x, i ), *getVectorDataComponent( y, i ) );

    } else {
        AMP_ERROR( "MultiVectorOperations::scale requires both x and y to be MultiVectorData" );
    }
}

void MultiVectorOperations::add( const VectorData &x, const VectorData &y, VectorData &z )
{
    if ( d_operations.empty() ) {
        return;
    }
    auto x2 = getMultiVectorData( x );
    auto y2 = getMultiVectorData( y );
    if ( x2 && y2 ) {
        auto z2 = getMultiVectorData( y );
        AMP_ASSERT( z2 );
        AMP_ASSERT( d_operations.size() == x2->numberOfComponents() );
        AMP_ASSERT( d_operations.size() == y2->numberOfComponents() );
        for ( size_t i = 0; i != d_operations.size(); i++ )
            d_operations[i]->add( *getVectorDataComponent( x, i ),
                                  *getVectorDataComponent( y, i ),
                                  *getVectorDataComponent( z, i ) );
    } else {
        AMP_ERROR( "MultiVectorOperations::add requires x, y, z to be MultiVectorData" );
    }
}

void MultiVectorOperations::subtract( const VectorData &x, const VectorData &y, VectorData &z )
{
    if ( d_operations.empty() ) {
        return;
    }
    auto x2 = getMultiVectorData( x );
    auto y2 = getMultiVectorData( y );
    if ( x2 && y2 ) {
        auto z2 = getMultiVectorData( y );
        AMP_ASSERT( z2 );
        AMP_ASSERT( d_operations.size() == x2->numberOfComponents() );
        AMP_ASSERT( d_operations.size() == y2->numberOfComponents() );
        for ( size_t i = 0; i != d_operations.size(); i++ )
            d_operations[i]->subtract( *getVectorDataComponent( x, i ),
                                       *getVectorDataComponent( y, i ),
                                       *getVectorDataComponent( z, i ) );
    } else {
        AMP_ERROR( "MultiVectorOperations::subtract requires x, y, z to be MultiVectorData" );
    }
}

void MultiVectorOperations::multiply( const VectorData &x, const VectorData &y, VectorData &z )
{
    if ( d_operations.empty() ) {
        return;
    }
    auto x2 = getMultiVectorData( x );
    auto y2 = getMultiVectorData( y );
    if ( x2 && y2 ) {
        auto z2 = getMultiVectorData( y );
        AMP_ASSERT( z2 );
        AMP_ASSERT( d_operations.size() == x2->numberOfComponents() );
        AMP_ASSERT( d_operations.size() == y2->numberOfComponents() );
        for ( size_t i = 0; i != d_operations.size(); i++ )
            d_operations[i]->multiply( *getVectorDataComponent( x, i ),
                                       *getVectorDataComponent( y, i ),
                                       *getVectorDataComponent( z, i ) );
    } else {
        AMP_ERROR( "MultiVectorOperations::multiply requires x, y, z to be MultiVectorData" );
    }
}

void MultiVectorOperations::divide( const VectorData &x, const VectorData &y, VectorData &z )
{
    if ( d_operations.empty() ) {
        return;
    }
    auto x2 = getMultiVectorData( x );
    auto y2 = getMultiVectorData( y );
    if ( x2 && y2 ) {
        auto z2 = getMultiVectorData( y );
        AMP_ASSERT( z2 );
        AMP_ASSERT( d_operations.size() == x2->numberOfComponents() );
        AMP_ASSERT( d_operations.size() == y2->numberOfComponents() );
        for ( size_t i = 0; i != d_operations.size(); i++ )
            d_operations[i]->divide( *getVectorDataComponent( x, i ),
                                     *getVectorDataComponent( y, i ),
                                     *getVectorDataComponent( z, i ) );
    } else {
        AMP_ERROR( "MultiVectorOperations::divide requires x, y, z to be MultiVectorData" );
    }
}

void MultiVectorOperations::reciprocal( const VectorData &x, VectorData &y )
{
    if ( d_operations.empty() ) {
        return;
    }
    auto x2 = getMultiVectorData( x );
    auto y2 = getMultiVectorData( y );
    if ( x2 && y2 ) {
        AMP_ASSERT( d_operations.size() == y2->numberOfComponents() );
        AMP_ASSERT( x2->numberOfComponents() == y2->numberOfComponents() );
        for ( size_t i = 0; i != d_operations.size(); i++ )
            d_operations[i]->reciprocal( *getVectorDataComponent( x, i ),
                                         *getVectorDataComponent( y, i ) );
    } else {
        AMP_ERROR(
            "MultiVectorOperations::reciprocal requires both x and y to be MultiVectorData" );
    }
}

void MultiVectorOperations::linearSum(
    double alpha_in, const VectorData &x, double beta_in, const VectorData &y, VectorData &z )
{
    if ( d_operations.empty() ) {
        return;
    }
    auto x2 = getMultiVectorData( x );
    auto y2 = getMultiVectorData( y );
    if ( x2 && y2 ) {
        auto z2 = getMultiVectorData( y );
        AMP_ASSERT( z2 );
        AMP_ASSERT( d_operations.size() == x2->numberOfComponents() );
        AMP_ASSERT( d_operations.size() == y2->numberOfComponents() );
        AMP_ASSERT( d_operations.size() == z2->numberOfComponents() );
        for ( size_t i = 0; i != d_operations.size(); i++ )
            d_operations[i]->linearSum( alpha_in,
                                        *getVectorDataComponent( x, i ),
                                        beta_in,
                                        *getVectorDataComponent( y, i ),
                                        *getVectorDataComponent( z, i ) );

    } else {
        AMP_ASSERT( x.getLocalSize() == y.getLocalSize() );
        AMP_ASSERT( x.getLocalSize() == z.getLocalSize() );
        if ( x.isType<double>() && y.isType<double>() ) {
            auto xit  = x.begin<double>();
            auto yit  = y.begin<double>();
            auto zit  = z.begin<double>();
            auto xend = x.end<double>();
            while ( xit != xend ) {
                *zit = alpha_in * ( *xit ) + beta_in * ( *yit );
                ++xit;
                ++yit;
                ++zit;
            }
        } else if ( x.isType<float>() && y.isType<float>() ) {
            auto xit  = x.begin<float>();
            auto yit  = y.begin<float>();
            auto zit  = z.begin<float>();
            auto xend = x.end<float>();
            while ( xit != xend ) {
                *zit = alpha_in * ( *xit ) + beta_in * ( *yit );
                ++xit;
                ++yit;
                ++zit;
            }
        } else {
            AMP_ERROR( "Unable to discern data types" );
        }
    }
}

void MultiVectorOperations::axpy( double alpha_in,
                                  const VectorData &x,
                                  const VectorData &y,
                                  VectorData &z )
{
    linearSum( alpha_in, x, 1.0, y, z );
}

void MultiVectorOperations::axpby( double alpha_in,
                                   double beta_in,
                                   const VectorData &x,
                                   VectorData &z )
{
    linearSum( alpha_in, x, beta_in, z, z );
}

void MultiVectorOperations::abs( const VectorData &x, VectorData &y )
{
    if ( d_operations.empty() ) {
        return;
    }
    auto x2 = getMultiVectorData( x );
    auto y2 = getMultiVectorData( y );
    if ( x2 && y2 ) {
        AMP_ASSERT( d_operations.size() == x2->numberOfComponents() );
        AMP_ASSERT( d_operations.size() == y2->numberOfComponents() );
        for ( size_t i = 0; i != d_operations.size(); i++ ) {
            d_operations[i]->abs( *getVectorDataComponent( x, i ),
                                  *getVectorDataComponent( y, i ) );
        }
    } else {
        AMP_ERROR( "MultiVectorOperations::abs requires x, y to be MultiVectorData" );
    }
}

void MultiVectorOperations::addScalar( const VectorData &x, double alpha_in, VectorData &y )
{
    if ( d_operations.empty() ) {
        return;
    }
    auto x2 = getMultiVectorData( x );
    auto y2 = getMultiVectorData( y );
    if ( x2 && y2 ) {
        AMP_ASSERT( d_operations.size() == x2->numberOfComponents() );
        AMP_ASSERT( d_operations.size() == y2->numberOfComponents() );
        for ( size_t i = 0; i != d_operations.size(); i++ )
            d_operations[i]->addScalar(
                *getVectorDataComponent( x, i ), alpha_in, *getVectorDataComponent( y, i ) );
    } else {
        AMP_ERROR( "MultiVectorOperations::addScalar requires x, y to be MultiVectorData" );
    }
}

double MultiVectorOperations::localMin( const VectorData &x ) const
{
    double ans = 1e300;
    AMP_ASSERT( getMultiVectorData( x ) );
    if ( d_operations.empty() ) {
        return 0;
    }
    for ( size_t i = 0; i != d_operations.size(); i++ )
        ans = std::min( ans, d_operations[i]->localMin( *getVectorDataComponent( x, i ) ) );
    return ans;
}

double MultiVectorOperations::localMax( const VectorData &x ) const
{
    double ans = -1e300;
    AMP_ASSERT( getMultiVectorData( x ) );
    if ( d_operations.empty() ) {
        return 0;
    }
    for ( size_t i = 0; i != d_operations.size(); i++ )
        ans = std::max( ans, d_operations[i]->localMax( *getVectorDataComponent( x, i ) ) );
    return ans;
}

double MultiVectorOperations::localL1Norm( const VectorData &x ) const
{
    double ans = 0.0;
    AMP_ASSERT( getMultiVectorData( x ) );
    if ( d_operations.empty() ) {
        return 0;
    }
    for ( size_t i = 0; i != d_operations.size(); i++ )
        ans += d_operations[i]->localL1Norm( *getVectorDataComponent( x, i ) );
    return ans;
}

double MultiVectorOperations::localL2Norm( const VectorData &x ) const
{
    double ans = 0.0;
    AMP_ASSERT( getMultiVectorData( x ) );
    if ( d_operations.empty() ) {
        return 0;
    }
    for ( size_t i = 0; i != d_operations.size(); i++ ) {
        const auto tmp = d_operations[i]->localL2Norm( *getVectorDataComponent( x, i ) );
        ans += tmp * tmp;
    }
    return sqrt( ans );
}

double MultiVectorOperations::localMaxNorm( const VectorData &x ) const
{
    double ans = 0.0;
    AMP_ASSERT( getMultiVectorData( x ) );
    if ( d_operations.empty() ) {
        return 0;
    }
    for ( size_t i = 0; i != d_operations.size(); i++ )
        ans = std::max( ans, d_operations[i]->localMaxNorm( *getVectorDataComponent( x, i ) ) );
    return ans;
}

double MultiVectorOperations::localDot( const VectorData &x, const VectorData &y ) const
{
    if ( d_operations.empty() ) {
        return 0;
    }

    double ans = 0.0;
    auto x2    = getMultiVectorData( x );
    auto y2    = getMultiVectorData( y );
    if ( x2 && y2 ) {
        AMP_ASSERT( d_operations.size() == x2->numberOfComponents() );
        AMP_ASSERT( d_operations.size() == y2->numberOfComponents() );
        for ( size_t i = 0; i != d_operations.size(); i++ )
            ans += d_operations[i]->localDot( *getVectorDataComponent( x, i ),
                                              *getVectorDataComponent( y, i ) );
    } else {
        AMP_ERROR( "MultiVectorOperations::localMinQuotient requires x, y to be MultiVectorData" );
    }
    return ans;
}

double MultiVectorOperations::localMinQuotient( const VectorData &x, const VectorData &y ) const
{
    if ( d_operations.empty() ) {
        return std::numeric_limits<double>::max();
    }

    double ans = std::numeric_limits<double>::max();
    auto x2    = getMultiVectorData( x );
    auto y2    = getMultiVectorData( y );
    if ( x2 && y2 ) {
        AMP_ASSERT( d_operations.size() == x2->numberOfComponents() );
        AMP_ASSERT( d_operations.size() == y2->numberOfComponents() );
        for ( size_t i = 0; i != d_operations.size(); i++ )
            ans = std::min( ans,
                            d_operations[i]->localMinQuotient( *getVectorDataComponent( x, i ),
                                                               *getVectorDataComponent( y, i ) ) );
    } else {
        AMP_ERROR( "MultiVectorOperations::localMinQuotient requires x, y to be MultiVectorData" );
    }
    return ans;
}

double MultiVectorOperations::localWrmsNorm( const VectorData &x, const VectorData &y ) const
{
    if ( d_operations.empty() ) {
        return 0;
    }
    double ans = 0.0;
    auto x2    = getMultiVectorData( x );
    auto y2    = getMultiVectorData( y );
    if ( x2 && y2 ) {
        AMP_ASSERT( d_operations.size() == x2->numberOfComponents() );
        AMP_ASSERT( d_operations.size() == y2->numberOfComponents() );
        for ( size_t i = 0; i < d_operations.size(); i++ ) {
            auto yi    = getVectorDataComponent( y, i );
            double tmp = d_operations[i]->localWrmsNorm( *getVectorDataComponent( x, i ), *yi );
            size_t N1  = yi->getLocalSize();
            ans += tmp * tmp * N1;
        }
        size_t N = y.getLocalSize();
        return sqrt( ans / N );
    } else {
        AMP_ERROR( "MultiVectorOperations::localWrmsNorm requires x, y to be MultiVectorData" );
    }
    return ans;
}

double MultiVectorOperations::localWrmsNormMask( const VectorData &x,
                                                 const VectorData &mask,
                                                 const VectorData &y ) const
{
    if ( d_operations.empty() ) {
        return 0;
    }
    double ans = 0.0;
    auto x2    = getMultiVectorData( x );
    auto m2    = getMultiVectorData( mask );
    auto y2    = getMultiVectorData( y );
    if ( x2 && m2 && y2 ) {
        AMP_ASSERT( d_operations.size() == x2->numberOfComponents() );
        AMP_ASSERT( d_operations.size() == m2->numberOfComponents() );
        AMP_ASSERT( d_operations.size() == y2->numberOfComponents() );
        for ( size_t i = 0; i < d_operations.size(); i++ ) {
            auto yi    = getVectorDataComponent( y, i );
            double tmp = d_operations[i]->localWrmsNormMask(
                *getVectorDataComponent( x, i ), *getVectorDataComponent( mask, i ), *yi );
            size_t N1 = yi->getLocalSize();
            ans += tmp * tmp * N1;
        }
        size_t N = y.getLocalSize();
        return sqrt( ans / N );
    } else {
        AMP_ERROR(
            "MultiVectorOperations::localWrmsNormMask requires x, mask, y to be MultiVectorData" );
    }
    return ans;
}

bool MultiVectorOperations::localEquals( const VectorData &x,
                                         const VectorData &y,
                                         double tol ) const
{
    if ( d_operations.empty() ) {
        return false;
    }
    bool ans = true;
    auto x2  = getMultiVectorData( x );
    auto y2  = getMultiVectorData( y );
    if ( x2 && y2 ) {
        AMP_ASSERT( d_operations.size() == x2->numberOfComponents() );
        AMP_ASSERT( d_operations.size() == y2->numberOfComponents() );
        for ( size_t i = 0; i < d_operations.size(); i++ ) {
            ans = ans && d_operations[i]->localEquals( *getVectorDataComponent( x, i ),
                                                       *getVectorDataComponent( y, i ),
                                                       tol );
        }
    } else {
        AMP_ERROR( "MultiVectorOperations::localEquals requires x, y to be MultiVectorData" );
    }
    return ans;
}

} // namespace LinearAlgebra
} // namespace AMP
