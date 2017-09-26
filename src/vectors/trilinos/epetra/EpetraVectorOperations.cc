#include "vectors/trilinos/epetra/EpetraVectorOperations.h"
#include "vectors/trilinos/epetra/EpetraVectorEngine.h"


namespace AMP {
namespace LinearAlgebra {


static inline Epetra_Vector &getEpetraVector( VectorOperations &vec )
{
    auto epetra = dynamic_cast<EpetraVectorEngine *>( &vec );
    AMP_INSIST( epetra != nullptr, "Not an EpetraVectorEngine" );
    return epetra->getEpetra_Vector();
}
static inline const Epetra_Vector &getEpetraVector( const VectorOperations &vec )
{
    auto epetra = dynamic_cast<const EpetraVectorEngine *>( &vec );
    AMP_INSIST( epetra != nullptr, "Not an EpetraVectorEngine" );
    return epetra->getEpetra_Vector();
}

Epetra_Vector &EpetraVectorOperations::getEpetra_Vector()
{
    auto epetra = dynamic_cast<EpetraVectorEngine *>( this );
    AMP_INSIST( epetra != nullptr, "Not an EpetraVectorEngine" );
    return epetra->getEpetra_Vector();
}
const Epetra_Vector &EpetraVectorOperations::getEpetra_Vector() const
{
    auto epetra = dynamic_cast<const EpetraVectorEngine *>( this );
    AMP_INSIST( epetra != nullptr, "Not an EpetraVectorEngine" );
    return epetra->getEpetra_Vector();
}

void EpetraVectorOperations::setToScalar( const double alpha )
{
    getEpetra_Vector().PutScalar( alpha );
}

void EpetraVectorOperations::scale( double alpha, const VectorOperations &x )
{
    getEpetra_Vector().Scale( alpha, getEpetraVector( x ) );
}

void EpetraVectorOperations::add( const VectorOperations &x, const VectorOperations &y )
{
    getEpetra_Vector().Update( 1., getEpetraVector( x ), 1., getEpetraVector( y ), 0. );
}

void EpetraVectorOperations::subtract( const VectorOperations &x, const VectorOperations &y )
{
    getEpetra_Vector().Update( 1., getEpetraVector( x ), -1., getEpetraVector( y ), 0. );
}

void EpetraVectorOperations::multiply( const VectorOperations &x, const VectorOperations &y )
{
    getEpetra_Vector().Multiply( 1., getEpetraVector( x ), getEpetraVector( y ), 0. );
}

void EpetraVectorOperations::divide( const VectorOperations &x, const VectorOperations &y )
{
    getEpetra_Vector().ReciprocalMultiply( 1., getEpetraVector( y ), getEpetraVector( x ), 0. );
}

void EpetraVectorOperations::reciprocal( const VectorOperations &x )
{
    getEpetra_Vector().Reciprocal( getEpetraVector( x ) );
}

void EpetraVectorOperations::linearSum( double alpha,
                                        const VectorOperations &x,
                                        double beta,
                                        const VectorOperations &y )
{
    getEpetra_Vector().Update( alpha, getEpetraVector( x ), beta, getEpetraVector( y ), 0. );
}

void EpetraVectorOperations::axpy( double alpha,
                                   const VectorOperations &x,
                                   const VectorOperations &y )
{
    linearSum( alpha, x, 1., y );
}

void EpetraVectorOperations::axpby( double alpha, double beta, const VectorOperations &x )
{
    linearSum( alpha, x, beta, *this );
}

void EpetraVectorOperations::abs( const VectorOperations &x )
{
    getEpetra_Vector().Abs( getEpetraVector( x ) );
}

double EpetraVectorOperations::min() const
{
    double retVal;
    getEpetra_Vector().MinValue( &retVal );
    return retVal;
}

double EpetraVectorOperations::max() const
{
    double retVal;
    getEpetra_Vector().MaxValue( &retVal );
    return retVal;
}

void EpetraVectorOperations::setRandomValues()
{
    getEpetra_Vector().Random();
    this->abs( *this );
}

void EpetraVectorOperations::scale( double alpha )
{
    getEpetra_Vector().Scale( alpha );
}

double EpetraVectorOperations::L1Norm() const
{
    double retVal;
    getEpetra_Vector().Norm1( &retVal );
    return retVal;
}

double EpetraVectorOperations::L2Norm() const
{
    double retVal;
    getEpetra_Vector().Norm2( &retVal );
    return retVal;
}


double EpetraVectorOperations::maxNorm() const
{
    double retVal;
    getEpetra_Vector().NormInf( &retVal );
    return retVal;
}

double EpetraVectorOperations::dot( const VectorOperations &x ) const
{
    double retVal;
    getEpetra_Vector().Dot( getEpetraVector( x ), &retVal );
    return retVal;
}

double EpetraVectorOperations::localDot( const VectorOperations & ) const
{
    AMP_ERROR( "Not implimented" );
    return 0;
}

double EpetraVectorOperations::localMin() const
{
    AMP_ERROR( "Not implimented" );
    return 0;
}
double EpetraVectorOperations::localMax() const
{
    AMP_ERROR( "Not implimented" );
    return 0;
}
double EpetraVectorOperations::localL1Norm() const
{
    AMP_ERROR( "Not implimented" );
    return 0;
}
double EpetraVectorOperations::localL2Norm() const
{
    AMP_ERROR( "Not implimented" );
    return 0;
}
double EpetraVectorOperations::localMaxNorm() const
{
    AMP_ERROR( "Not implimented" );
    return 0;
}


} // namespace LinearAlgebra
} // namespace AMP
