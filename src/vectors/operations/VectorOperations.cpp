#include "AMP/vectors/operations/VectorOperations.h"
#include "AMP/IO/RestartManager.h"
#include "AMP/vectors/Vector.h"
#include "AMP/vectors/data/VectorData.h"


namespace AMP::LinearAlgebra {


/****************************************************************
 * Constructors                                                  *
 ****************************************************************/
VectorOperations::VectorOperations() = default;


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
        ans = x.getComm().minReduce( ans );
    return ans;
}
Scalar VectorOperations::max( const VectorData &x ) const
{
    auto ans = localMax( x );
    if ( x.hasComm() )
        ans = x.getComm().maxReduce( ans );
    return ans;
}
Scalar VectorOperations::sum( const VectorData &x ) const
{
    auto ans = localSum( x );
    if ( x.hasComm() )
        ans = x.getComm().sumReduce( ans );
    return ans;
}
Scalar VectorOperations::mean( const VectorData &x ) const
{
    return sum( x ) / Scalar( x.getGlobalSize() );
}
Scalar VectorOperations::dot( const VectorData &x, const VectorData &y ) const
{
    auto ans = localDot( x, y );
    if ( x.hasComm() )
        ans = x.getComm().sumReduce( ans );
    return ans;
}
Scalar VectorOperations::L1Norm( const VectorData &x ) const
{
    Scalar ans = localL1Norm( x );
    if ( x.hasComm() )
        ans = x.getComm().sumReduce( ans );
    return ans;
}
Scalar VectorOperations::maxNorm( const VectorData &x ) const
{
    Scalar ans = localMaxNorm( x );
    if ( x.hasComm() )
        ans = x.getComm().maxReduce( ans );
    return ans;
}
Scalar VectorOperations::L2Norm( const VectorData &x ) const
{
    auto ans = localL2Norm( x );
    if ( x.hasComm() )
        ans = x.getComm().sumReduce( ans * ans ).sqrt();
    return ans;
}
Scalar VectorOperations::minQuotient( const VectorData &x, const VectorData &y ) const
{
    auto ans = localMinQuotient( x, y );
    if ( y.getCommunicationList() )
        ans = y.getComm().minReduce( ans );
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
        ans       = x.getComm().sumReduce( tmp );
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
        ans       = x.getComm().sumReduce( tmp );
        ans       = ans.sqrt();
    }
    return ans;
}


/****************************************************************
 * Get an id                                                     *
 ****************************************************************/
uint64_t VectorOperations::getID() const
{
    AMP_ERROR( "Not finished" );
    return 0;
}


} // namespace AMP::LinearAlgebra


/********************************************************
 *  Restart operations                                   *
 ********************************************************/
template<>
AMP::IO::RestartManager::DataStoreType<AMP::LinearAlgebra::VectorOperations>::DataStoreType(
    const std::string &name,
    std::shared_ptr<const AMP::LinearAlgebra::VectorOperations> data,
    RestartManager * )
    : d_data( data )
{
    d_name = name;
    d_hash = data->getID();
    AMP_ERROR( "Not finished" );
}
template<>
void AMP::IO::RestartManager::DataStoreType<AMP::LinearAlgebra::VectorOperations>::write(
    hid_t fid, const std::string &name ) const
{
    hid_t gid = createGroup( fid, name );
    AMP_ERROR( "Not finished" );
    closeGroup( gid );
}
template<>
std::shared_ptr<AMP::LinearAlgebra::VectorOperations>
AMP::IO::RestartManager::getData<AMP::LinearAlgebra::VectorOperations>( const std::string &name )
{
    hid_t gid = openGroup( d_fid, name );
    AMP_ERROR( "Not finished" );
    closeGroup( gid );
}
