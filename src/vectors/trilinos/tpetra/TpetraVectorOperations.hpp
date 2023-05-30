#ifndef included_AMP_NativeTpetraVectorOperations_HPP_
#define included_AMP_NativeTpetraVectorOperations_HPP_

#include "AMP/vectors/trilinos/tpetra/TpetraVectorData.h"
#include "AMP/vectors/trilinos/tpetra/TpetraVectorOperations.h"

#include <Kokkos_Core.hpp>

namespace AMP::LinearAlgebra {


template<typename ST, typename LO, typename GO, typename NT>
static inline const Tpetra::Vector<ST, LO, GO, NT> &getTpetraVector( const VectorData &vec )
{
    auto tpetraData = dynamic_cast<const TpetraVectorData<ST, LO, GO, NT> *>( &vec );
    AMP_INSIST( tpetraData, "Not TpetraVectorData" );
    return *( tpetraData->getTpetraVector() );
}

template<typename ST, typename LO, typename GO, typename NT>
static inline Tpetra::Vector<ST, LO, GO, NT> &getTpetraVector( VectorData &vec )
{
    auto data = dynamic_cast<TpetraVectorData<ST, LO, GO, NT> *>( &vec );
    AMP_INSIST( data, "Not TpetraVectorData" );
    return *( data->getTpetraVector() );
}

template<typename ST, typename LO, typename GO, typename NT>
void TpetraVectorOperations<ST, LO, GO, NT>::setToScalar( const AMP::Scalar &alpha, VectorData &x )
{
    getTpetraVector<ST, LO, GO, NT>( x ).putScalar( static_cast<ST>( alpha ) );
}

template<typename ST, typename LO, typename GO, typename NT>
void TpetraVectorOperations<ST, LO, GO, NT>::addScalar( const VectorData &,
                                                        const Scalar &,
                                                        VectorData & )
{
    AMP_ERROR( "TpetraVectorOperations::addScalar not implemented" );
}

template<typename ST, typename LO, typename GO, typename NT>
void TpetraVectorOperations<ST, LO, GO, NT>::setRandomValues( VectorData &x )
{
    getTpetraVector<ST, LO, GO, NT>( x ).randomize();
}

template<typename ST, typename LO, typename GO, typename NT>
void TpetraVectorOperations<ST, LO, GO, NT>::copy( const VectorData &x, VectorData &z )
{
    deep_copy( getTpetraVector<ST, LO, GO, NT>( z ), getTpetraVector<ST, LO, GO, NT>( x ) );
}

template<typename ST, typename LO, typename GO, typename NT>
void TpetraVectorOperations<ST, LO, GO, NT>::scale( const Scalar &alpha,
                                                    const VectorData &x,
                                                    VectorData &y )
{
    getTpetraVector<ST, LO, GO, NT>( y ).scale( static_cast<ST>( alpha ),
                                                getTpetraVector<ST, LO, GO, NT>( x ) );
}

template<typename ST, typename LO, typename GO, typename NT>
void TpetraVectorOperations<ST, LO, GO, NT>::scale( const Scalar &alpha, VectorData &x )
{
    getTpetraVector<ST, LO, GO, NT>( x ).scale( static_cast<ST>( alpha ) );
}

template<typename ST, typename LO, typename GO, typename NT>
void TpetraVectorOperations<ST, LO, GO, NT>::add( const VectorData &x,
                                                  const VectorData &y,
                                                  VectorData &z )
{
    getTpetraVector<ST, LO, GO, NT>( z ).update( static_cast<ST>( 1.0 ),
                                                 getTpetraVector<ST, LO, GO, NT>( x ),
                                                 static_cast<ST>( 1.0 ),
                                                 getTpetraVector<ST, LO, GO, NT>( y ),
                                                 static_cast<ST>( 0.0 ) );
}

template<typename ST, typename LO, typename GO, typename NT>
void TpetraVectorOperations<ST, LO, GO, NT>::subtract( const VectorData &x,
                                                       const VectorData &y,
                                                       VectorData &z )
{
    getTpetraVector<ST, LO, GO, NT>( z ).update( static_cast<ST>( 1.0 ),
                                                 getTpetraVector<ST, LO, GO, NT>( x ),
                                                 static_cast<ST>( -1.0 ),
                                                 getTpetraVector<ST, LO, GO, NT>( y ),
                                                 static_cast<ST>( 0.0 ) );
}

template<typename ST, typename LO, typename GO, typename NT>
void TpetraVectorOperations<ST, LO, GO, NT>::multiply( const VectorData &x,
                                                       const VectorData &y,
                                                       VectorData &z )
{
    getTpetraVector<ST, LO, GO, NT>( z ).elementWiseMultiply( static_cast<ST>( 1.0 ),
                                                              getTpetraVector<ST, LO, GO, NT>( x ),
                                                              getTpetraVector<ST, LO, GO, NT>( y ),
                                                              static_cast<ST>( 0.0 ) );
}

template<typename ST, typename LO, typename GO, typename NT>
void TpetraVectorOperations<ST, LO, GO, NT>::divide( const VectorData &x,
                                                     const VectorData &y,
                                                     VectorData &z )
{
    getTpetraVector<ST, LO, GO, NT>( z ).reciprocal( getTpetraVector<ST, LO, GO, NT>( y ) );
    getTpetraVector<ST, LO, GO, NT>( z ).elementWiseMultiply( static_cast<ST>( 1.0 ),
                                                              getTpetraVector<ST, LO, GO, NT>( x ),
                                                              getTpetraVector<ST, LO, GO, NT>( z ),
                                                              static_cast<ST>( 0.0 ) );
}

template<typename ST, typename LO, typename GO, typename NT>
void TpetraVectorOperations<ST, LO, GO, NT>::reciprocal( const VectorData &x, VectorData &y )
{
    getTpetraVector<ST, LO, GO, NT>( y ).reciprocal( getTpetraVector<ST, LO, GO, NT>( x ) );
}

template<typename ST, typename LO, typename GO, typename NT>
void TpetraVectorOperations<ST, LO, GO, NT>::linearSum( const Scalar &alpha,
                                                        const VectorData &x,
                                                        const Scalar &beta,
                                                        const VectorData &y,
                                                        VectorData &z )
{
    getTpetraVector<ST, LO, GO, NT>( z ).update( static_cast<ST>( alpha ),
                                                 getTpetraVector<ST, LO, GO, NT>( x ),
                                                 static_cast<ST>( beta ),
                                                 getTpetraVector<ST, LO, GO, NT>( y ),
                                                 static_cast<ST>( 0.0 ) );
}

template<typename ST, typename LO, typename GO, typename NT>
void TpetraVectorOperations<ST, LO, GO, NT>::axpy( const Scalar &alpha,
                                                   const VectorData &x,
                                                   const VectorData &y,
                                                   VectorData &z )
{
    getTpetraVector<ST, LO, GO, NT>( z ).update( static_cast<ST>( alpha ),
                                                 getTpetraVector<ST, LO, GO, NT>( x ),
                                                 static_cast<ST>( 1.0 ),
                                                 getTpetraVector<ST, LO, GO, NT>( y ),
                                                 static_cast<ST>( 0.0 ) );
}

template<typename ST, typename LO, typename GO, typename NT>
void TpetraVectorOperations<ST, LO, GO, NT>::axpby( const Scalar &alpha,
                                                    const Scalar &beta,
                                                    const VectorData &x,
                                                    VectorData &y )
{
    getTpetraVector<ST, LO, GO, NT>( y ).update(
        static_cast<ST>( alpha ), getTpetraVector<ST, LO, GO, NT>( x ), static_cast<ST>( beta ) );
}

template<typename ST, typename LO, typename GO, typename NT>
void TpetraVectorOperations<ST, LO, GO, NT>::abs( const VectorData &x, VectorData &z )
{
    getTpetraVector<ST, LO, GO, NT>( z ).abs( getTpetraVector<ST, LO, GO, NT>( x ) );
}

template<typename ST, typename LO, typename GO, typename NT>
Scalar TpetraVectorOperations<ST, LO, GO, NT>::min( const VectorData &x ) const
{
    AMP_ERROR( "Not implemented" );
    return 0;
}

template<typename ST, typename LO, typename GO, typename NT>
Scalar TpetraVectorOperations<ST, LO, GO, NT>::max( const VectorData &x ) const
{
    AMP_ERROR( "Not implemented" );
    return 0;
}

template<typename ST, typename LO, typename GO, typename NT>
Scalar TpetraVectorOperations<ST, LO, GO, NT>::L1Norm( const VectorData &x ) const
{
    return getTpetraVector<ST, LO, GO, NT>( x ).norm1();
}

template<typename ST, typename LO, typename GO, typename NT>
Scalar TpetraVectorOperations<ST, LO, GO, NT>::L2Norm( const VectorData &x ) const
{
    return getTpetraVector<ST, LO, GO, NT>( x ).norm2();
}

template<typename ST, typename LO, typename GO, typename NT>
Scalar TpetraVectorOperations<ST, LO, GO, NT>::maxNorm( const VectorData &x ) const
{
    return getTpetraVector<ST, LO, GO, NT>( x ).normInf();
}

template<typename ST, typename LO, typename GO, typename NT>
Scalar TpetraVectorOperations<ST, LO, GO, NT>::dot( const VectorData &x, const VectorData &y ) const
{
    return getTpetraVector<ST, LO, GO, NT>( x ).dot( getTpetraVector<ST, LO, GO, NT>( y ) );
}

template<typename ST, typename LO, typename GO, typename NT>
static Teuchos::RCP<Tpetra::Map<LO, GO, NT>>
getLocalMap( Teuchos::RCP<const Tpetra::Map<LO, GO, NT>> map )
{
    Teuchos::RCP<Tpetra::Map<LO, GO, NT>> lmap(
        new Tpetra::Map<LO, GO, NT>( map->getLocalNumElements(),
                                     map->getIndexBase(),
                                     map->getComm(),
                                     Tpetra::LocallyReplicated ) );
    return lmap;
}

template<typename ST, typename LO, typename GO, typename NT>
Scalar TpetraVectorOperations<ST, LO, GO, NT>::localL1Norm( const VectorData &x ) const
{
    const auto &xt   = getTpetraVector<ST, LO, GO, NT>( x );
    const auto &lmap = getLocalMap<ST, LO, GO, NT>( xt.getMap() );
    const auto &xl   = xt.offsetView( lmap, 0 );
    return xl->norm1();
}

template<typename ST, typename LO, typename GO, typename NT>
Scalar TpetraVectorOperations<ST, LO, GO, NT>::localL2Norm( const VectorData &x ) const
{
    const auto &xt   = getTpetraVector<ST, LO, GO, NT>( x );
    const auto &lmap = getLocalMap<ST, LO, GO, NT>( xt.getMap() );
    const auto &xl   = xt.offsetView( lmap, 0 );
    return xl->norm2();
}

template<typename ST, typename LO, typename GO, typename NT>
Scalar TpetraVectorOperations<ST, LO, GO, NT>::localMaxNorm( const VectorData &x ) const
{
    const auto &xt   = getTpetraVector<ST, LO, GO, NT>( x );
    const auto &lmap = getLocalMap<ST, LO, GO, NT>( xt.getMap() );
    const auto &xl   = xt.offsetView( lmap, 0 );
    return xl->normInf();
}

template<typename ST, typename LO, typename GO, typename NT>
Scalar TpetraVectorOperations<ST, LO, GO, NT>::localDot( const VectorData &x,
                                                         const VectorData &y ) const
{
    const auto &xt    = getTpetraVector<ST, LO, GO, NT>( x );
    const auto &xlmap = getLocalMap<ST, LO, GO, NT>( xt.getMap() );
    const auto &xl    = xt.offsetView( xlmap, 0 );

    const auto &yt    = getTpetraVector<ST, LO, GO, NT>( y );
    const auto &ylmap = getLocalMap<ST, LO, GO, NT>( yt.getMap() );
    const auto &yl    = yt.offsetView( ylmap, 0 );

    return xl->dot( *yl );
}

template<typename ST, typename LO, typename GO, typename NT>
Scalar
TpetraVectorOperations<ST, LO, GO, NT>::localMin( const AMP::LinearAlgebra::VectorData & ) const
{
    AMP_ERROR( "TpetraVectorOperations::localMax not implemented" );
    return 0;
}

template<typename ST, typename LO, typename GO, typename NT>
Scalar
TpetraVectorOperations<ST, LO, GO, NT>::localMax( const AMP::LinearAlgebra::VectorData & ) const
{
    AMP_ERROR( "TpetraVectorOperations::localMax not implemented" );
    return 0;
}

template<typename ST, typename LO, typename GO, typename NT>
Scalar
TpetraVectorOperations<ST, LO, GO, NT>::localSum( const AMP::LinearAlgebra::VectorData & ) const
{
    AMP_ERROR( "TpetraVectorOperations::localSum not implemented" );
    return 0;
}

template<typename ST, typename LO, typename GO, typename NT>
Scalar TpetraVectorOperations<ST, LO, GO, NT>::localMinQuotient(
    const AMP::LinearAlgebra::VectorData &, const AMP::LinearAlgebra::VectorData & ) const
{
    AMP_ERROR( "TpetraVectorOperations::localMinQuotient not implemented" );
    return 0;
}

template<typename ST, typename LO, typename GO, typename NT>
Scalar TpetraVectorOperations<ST, LO, GO, NT>::localWrmsNorm(
    const AMP::LinearAlgebra::VectorData &, const AMP::LinearAlgebra::VectorData & ) const
{
    AMP_ERROR( "TpetraVectorOperations::localWrmsNorm not implemented" );
}

template<typename ST, typename LO, typename GO, typename NT>
Scalar TpetraVectorOperations<ST, LO, GO, NT>::localWrmsNormMask(
    const AMP::LinearAlgebra::VectorData &,
    const AMP::LinearAlgebra::VectorData &,
    const AMP::LinearAlgebra::VectorData & ) const
{
    AMP_ERROR( "TpetraVectorOperations::localWrmsNormMask not implemented" );
    return 0;
}

template<typename ST, typename LO, typename GO, typename NT>
bool TpetraVectorOperations<ST, LO, GO, NT>::localEquals( const AMP::LinearAlgebra::VectorData &,
                                                          const AMP::LinearAlgebra::VectorData &,
                                                          const Scalar & ) const
{
    AMP_ERROR( "TpetraVectorOperations::localEquals not implemented" );
    return false;
}

} // namespace AMP::LinearAlgebra

#endif
