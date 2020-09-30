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

double VectorOperations::min( const VectorData &x ) const
{
    double ans = localMin( x );
    if ( x.hasComm() )
        ans = x.getComm().minReduce( ans );
    return ans;
}
double VectorOperations::max( const VectorData &x ) const
{
    double ans = localMax( x );
    if ( x.hasComm() )
        ans = x.getComm().maxReduce( ans );
    return ans;
}
double VectorOperations::dot( const VectorData &x, const VectorData &y ) const
{
    double ans = localDot( x, y );
    if ( x.hasComm() )
        ans = x.getComm().sumReduce( ans );
    return ans;
}
double VectorOperations::L1Norm( const VectorData &x ) const
{
    double ans = localL1Norm( x );
    if ( x.hasComm() )
        ans = x.getComm().sumReduce( ans );
    return ans;
}
double VectorOperations::maxNorm( const VectorData &x ) const
{
    double ans = localMaxNorm( x );
    if ( x.hasComm() )
        ans = x.getComm().maxReduce( ans );
    return ans;
}
double VectorOperations::L2Norm( const VectorData &x ) const
{
    double ans = localL2Norm( x );
    if ( x.hasComm() )
        ans = sqrt( x.getComm().sumReduce( ans * ans ) );
    return ans;
}
double VectorOperations::minQuotient( const VectorData &x, const VectorData &y ) const
{
    double ans = localMinQuotient( x, y );
    if ( y.getCommunicationList() )
        ans = y.getComm().minReduce( ans );
    AMP_INSIST( ans < std::numeric_limits<double>::max(),
                "denominator is the zero vector on an entire process" );
    return ans;
}

double VectorOperations::wrmsNorm( const VectorData &x, const VectorData &y ) const
{
    double ans = localWrmsNorm( x, y );
    if ( y.getCommunicationList() ) {
        size_t N1 = y.getCommunicationList()->numLocalRows();
        size_t N2 = y.getCommunicationList()->getTotalSize();
        ans       = sqrt( y.getComm().sumReduce( ans * ans * N1 ) / N2 );
    }
    return ans;
}
double VectorOperations::wrmsNormMask( const VectorData &x,
                                       const VectorData &mask,
                                       const VectorData &y ) const
{
    double ans = localWrmsNormMask( x, mask, y );
    if ( y.getCommunicationList() ) {
        size_t N1 = y.getCommunicationList()->numLocalRows();
        size_t N2 = y.getCommunicationList()->getTotalSize();
        ans       = sqrt( y.getComm().sumReduce( ans * ans * N1 ) / N2 );
    }
    return ans;
}

} // namespace LinearAlgebra
} // namespace AMP
