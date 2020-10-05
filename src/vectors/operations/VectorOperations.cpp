#include "AMP/vectors/operations/VectorOperations.h"
#include "AMP/vectors/Vector.h"
#include "AMP/vectors/data/VectorData.h"


namespace AMP {
namespace LinearAlgebra {


/****************************************************************
 * Constructors                                                  *
 ****************************************************************/
VectorOperations::VectorOperations() {}

/****************************************************************
 * equals                                                        *
 * Note: these routines require communication                    *
 ****************************************************************/
bool VectorOperations::equals( const VectorData &a, const VectorData &b, const Scalar &tol ) const
{
    bool equal = localEquals( a, b, tol );
    if ( b.hasComm() )
        equal = b.getComm().allReduce( equal );
    return equal;
}

Scalar VectorOperations::min( const VectorData &x ) const
{
    auto ans = localMin( x );
    if ( x.hasComm() )
        ans = minReduce( x.getComm(), ans );
    return ans;
}
Scalar VectorOperations::max( const VectorData &x ) const
{
    auto ans = localMax( x );
    if ( x.hasComm() )
        ans = maxReduce( x.getComm(), ans );
    return ans;
}
Scalar VectorOperations::dot( const VectorData &x, const VectorData &y ) const
{
    auto ans = localDot( x, y );
    if ( x.hasComm() )
        ans = sumReduce( x.getComm(), ans );
    return ans;
}
Scalar VectorOperations::L1Norm( const VectorData &x ) const
{
    Scalar ans = localL1Norm( x );
    if ( x.hasComm() )
        ans = sumReduce( x.getComm(), ans );
    return ans;
}
Scalar VectorOperations::maxNorm( const VectorData &x ) const
{
    Scalar ans = localMaxNorm( x );
    if ( x.hasComm() )
        ans = maxReduce( x.getComm(), ans );
    return ans;
}
Scalar VectorOperations::L2Norm( const VectorData &x ) const
{
    auto ans = localL2Norm( x );
    if ( x.hasComm() )
        ans = sumReduce( x.getComm(), ans * ans ).sqrt();
    return ans;
}
Scalar VectorOperations::minQuotient( const VectorData &x, const VectorData &y ) const
{
    auto ans = localMinQuotient( x, y );
    if ( y.getCommunicationList() )
        ans = minReduce( y.getComm(), ans );
    AMP_INSIST( ans < std::numeric_limits<double>::max(),
                "denominator is the zero vector on an entire process" );
    return ans;
}
Scalar VectorOperations::wrmsNorm( const VectorData &x, const VectorData &y ) const
{
    auto ans = localWrmsNorm( x, y );
    if ( y.getCommunicationList() ) {
        double N1 = y.getCommunicationList()->numLocalRows();
        double N2 = y.getCommunicationList()->getTotalSize();
        auto tmp  = ans * ans * ( N1 / N2 );
        ans       = sumReduce( x.getComm(), tmp );
        ans       = ans.sqrt();
    }
    return ans;
}
Scalar VectorOperations::wrmsNormMask( const VectorData &x,
                                       const VectorData &mask,
                                       const VectorData &y ) const
{
    auto ans = localWrmsNormMask( x, mask, y );
    if ( y.getCommunicationList() ) {
        double N1 = y.getCommunicationList()->numLocalRows();
        double N2 = y.getCommunicationList()->getTotalSize();
        auto tmp  = ans * ans * ( N1 / N2 );
        ans       = sumReduce( x.getComm(), tmp );
        ans       = ans.sqrt();
    }
    return ans;
}

} // namespace LinearAlgebra
} // namespace AMP
