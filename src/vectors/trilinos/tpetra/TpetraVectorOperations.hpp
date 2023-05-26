#ifndef included_AMP_NativeTpetraVectorOperations_HPP_
#define included_AMP_NativeTpetraVectorOperations_HPP_

#include "AMP/vectors/trilinos/tpetra/TpetraVectorData.h"
#include "AMP/vectors/trilinos/tpetra/TpetraVectorOperations.h"

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

} // namespace AMP::LinearAlgebra

#endif
