#include "AMP/vectors/trilinos/epetra/EpetraVectorOperations.h"
#include "AMP/vectors/trilinos/epetra/EpetraVectorData.h"

namespace AMP {
namespace LinearAlgebra {

static inline const Epetra_Vector &getEpetraVector( const VectorData &vec )
{
    auto epetraData = dynamic_cast<const EpetraVectorData *>( &vec );
    AMP_INSIST( epetraData, "Not EpetraVectorData" );
    return epetraData->getEpetra_Vector();
}
static inline Epetra_Vector &getEpetraVector( VectorData &vec )
{
    auto epetraData = dynamic_cast<EpetraVectorData *>( &vec );
    AMP_INSIST( epetraData, "Not EpetraVectorData" );
    return epetraData->getEpetra_Vector();
}

//**********************************************************************
// Functions that operate on VectorData objects

void EpetraVectorOperations::setToScalar( const Scalar &alpha, VectorData &x )
{
    getEpetraVector( x ).PutScalar( alpha.get<double>() );
}

void EpetraVectorOperations::setRandomValues( VectorData &x )
{
    getEpetraVector( x ).Random();
    abs( x, x );
}

void EpetraVectorOperations::scale( const Scalar &alpha, const VectorData &x, VectorData &y )
{
    getEpetraVector( y ).Scale( alpha.get<double>(), getEpetraVector( x ) );
}

void EpetraVectorOperations::scale( const Scalar &alpha, VectorData &x )
{
    getEpetraVector( x ).Scale( alpha.get<double>() );
}

void EpetraVectorOperations::add( const VectorData &x, const VectorData &y, VectorData &z )
{
    getEpetraVector( z ).Update( 1., getEpetraVector( x ), 1., getEpetraVector( y ), 0. );
}

void EpetraVectorOperations::subtract( const VectorData &x, const VectorData &y, VectorData &z )
{
    getEpetraVector( z ).Update( 1., getEpetraVector( x ), -1., getEpetraVector( y ), 0. );
}

void EpetraVectorOperations::multiply( const VectorData &x, const VectorData &y, VectorData &z )
{
    getEpetraVector( z ).Multiply( 1., getEpetraVector( x ), getEpetraVector( y ), 0. );
}

void EpetraVectorOperations::divide( const VectorData &x, const VectorData &y, VectorData &z )
{
    getEpetraVector( z ).ReciprocalMultiply( 1., getEpetraVector( y ), getEpetraVector( x ), 0. );
}

void EpetraVectorOperations::reciprocal( const VectorData &x, VectorData &y )
{
    getEpetraVector( y ).Reciprocal( getEpetraVector( x ) );
}

void EpetraVectorOperations::linearSum( const Scalar &alpha,
                                        const VectorData &x,
                                        const Scalar &beta,
                                        const VectorData &y,
                                        VectorData &z )
{
    getEpetraVector( z ).Update(
        alpha.get<double>(), getEpetraVector( x ), beta.get<double>(), getEpetraVector( y ), 0. );
}

void EpetraVectorOperations::axpy( const Scalar &alpha,
                                   const VectorData &x,
                                   const VectorData &y,
                                   VectorData &z )
{
    linearSum( alpha.get<double>(), x, 1.0, y, z );
}

void EpetraVectorOperations::axpby( const Scalar &alpha,
                                    const Scalar &beta,
                                    const VectorData &x,
                                    VectorData &z )
{
    linearSum( alpha.get<double>(), x, beta.get<double>(), z, z );
}

void EpetraVectorOperations::abs( const VectorData &x, VectorData &y )
{
    getEpetraVector( y ).Abs( getEpetraVector( x ) );
}

Scalar EpetraVectorOperations::min( const VectorData &x ) const
{
    double retVal;
    getEpetraVector( x ).MinValue( &retVal );
    return retVal;
}

Scalar EpetraVectorOperations::max( const VectorData &x ) const
{
    double retVal;
    getEpetraVector( x ).MaxValue( &retVal );
    return retVal;
}

Scalar EpetraVectorOperations::dot( const VectorData &x, const VectorData &y ) const
{
    double retVal;
    getEpetraVector( y ).Dot( getEpetraVector( x ), &retVal );
    return retVal;
}

Scalar EpetraVectorOperations::L1Norm( const VectorData &x ) const
{
    double retVal;
    getEpetraVector( x ).Norm1( &retVal );
    return retVal;
}

Scalar EpetraVectorOperations::L2Norm( const VectorData &x ) const
{
    double retVal;
    getEpetraVector( x ).Norm2( &retVal );
    return retVal;
}

Scalar EpetraVectorOperations::maxNorm( const VectorData &x ) const
{
    double retVal;
    getEpetraVector( x ).NormInf( &retVal );
    return retVal;
}

} // namespace LinearAlgebra
} // namespace AMP
