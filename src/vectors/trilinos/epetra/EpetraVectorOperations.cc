#include "AMP/vectors/trilinos/epetra/EpetraVectorOperations.h"
#include "AMP/vectors/trilinos/epetra/EpetraVectorEngine.h"
#include "AMP/vectors/trilinos/epetra/NativeEpetraVector.h"


namespace AMP {
namespace LinearAlgebra {


static inline const Epetra_Vector &getEpetraVector( const VectorOperations &vec )
{
    auto epetraEngine = dynamic_cast<const EpetraVectorEngine *>( &vec );
    if(epetraEngine)
      return epetraEngine->getEpetra_Vector();
    else {
      auto epetraVec = dynamic_cast<const NativeEpetraVector *>( &vec );      
      AMP_INSIST( epetraVec != nullptr, "Not a NativeEpetraVector" );
      return epetraVec->getEpetra_Vector();
    }
}

Epetra_Vector &EpetraVectorOperations::getEpetra_Vector()
{
    auto epetraEngine = dynamic_cast<EpetraVectorEngine *>( this );
    if(epetraEngine)
      return epetraEngine->getEpetra_Vector();
    else {
      auto epetraVec = dynamic_cast<NativeEpetraVector *>( this );      
      AMP_INSIST( epetraVec != nullptr, "Not a NativeEpetraVector" );
      return epetraVec->getEpetra_Vector();
    }
}
const Epetra_Vector &EpetraVectorOperations::getEpetra_Vector() const
{
    auto epetraEngine = dynamic_cast<const EpetraVectorEngine *>( this );
    if(epetraEngine)
      return epetraEngine->getEpetra_Vector();
    else {
      auto epetraVec = dynamic_cast<const NativeEpetraVector *>( this );      
      AMP_INSIST( epetraVec != nullptr, "Not a NativeEpetraVector" );
      return epetraVec->getEpetra_Vector();
    }
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

void EpetraVectorOperations::scale( double alpha ) { getEpetra_Vector().Scale( alpha ); }

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

} // namespace LinearAlgebra
} // namespace AMP
