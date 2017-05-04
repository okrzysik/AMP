#include "vectors/VectorOperations.h"
#include "vectors/VectorData.h"
#include "vectors/Vector.h"


namespace AMP {
namespace LinearAlgebra {


/****************************************************************
* Constructors                                                  *
****************************************************************/
VectorOperations::VectorOperations():
     d_VectorData(nullptr)
{
}


/****************************************************************
* min, max, norms, etc.                                         *
* Note: these routines require communication                    *
****************************************************************/
double VectorOperations::min() const
{
    double ans = localMin();
    if ( hasComm() )
        ans = getComm().minReduce( ans );
    return ans;
}
double VectorOperations::max() const
{
    double ans = localMax();
    if ( hasComm() )
        ans = getComm().maxReduce( ans );
    return ans;
}
double VectorOperations::dot( const VectorOperations &x ) const
{
    double ans = x.localDot( x );
    if ( hasComm() )
        ans = getComm().sumReduce( ans );
    return ans;
}
double VectorOperations::L1Norm( void ) const
{
    double ans = localL1Norm();
    if ( hasComm() )
        ans = getComm().sumReduce( ans );
    return ans;
}
double VectorOperations::maxNorm( void ) const
{
    double ans = localMaxNorm();
    if ( hasComm() )
        ans = getComm().maxReduce( ans );
    return ans;
}
double VectorOperations::L2Norm( void ) const
{
    double ans = localL2Norm();
    if ( hasComm() )
        ans = sqrt( getComm().sumReduce( ans * ans ) );
    return ans;
}


} // LinearAlgebra namespace
} // AMP namespace

