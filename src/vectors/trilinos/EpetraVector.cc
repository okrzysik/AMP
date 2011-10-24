#include "EpetraVector.h"
#include "ManagedEpetraVector.h"
#include "MultiVector.h"
#include "CommCollectVector.h"


namespace AMP {
namespace LinearAlgebra {

  const Vector::shared_ptr  EpetraVector::constView ( const Vector::shared_ptr inVector )
  {
    if ( inVector->isA<CommCollectVector>() )
    {
      return constView ( inVector->castTo<const CommCollectVector>().getSmallCommVector() );
    }
    if ( inVector->isA<EpetraVector> () )
    {
      return inVector;
    }
    else if ( inVector->isA<ManagedVector> () )
    {
      Vector::shared_ptr retVal;
      retVal = Vector::shared_ptr ( new ManagedEpetraVector ( inVector ) );
      return retVal;
    }
    else if ( inVector->isA<MultiVector> () )
    {
      if ( inVector->numberOfDataBlocks() == 1 )
      {
        return constView ( inVector->castTo<MultiVector>().getVector ( 0 ) );
      }
    }
    AMP_ERROR( "Cannot create view!" );
    return Vector::shared_ptr ();
  }

  Vector::shared_ptr  EpetraVector::view ( Vector::shared_ptr inVector )
  {
    Vector::shared_ptr  retVal;

    if ( inVector->isA<EpetraVector> () )
    {
      retVal = inVector;
    }
    else if ( inVector->isA<CommCollectVector>() )
    {
      retVal = view ( inVector->castTo<CommCollectVector>().getSmallCommVector() );
    }
    else if ( inVector->isA<MultiVector> () )
    {
      if ( inVector->numberOfDataBlocks() == 1 )
        retVal = view ( inVector->castTo<MultiVector>().getVector ( 0 ) );
    }
    else if ( inVector->isA<ManagedVector> () )
    {
      retVal = Vector::shared_ptr ( new ManagedEpetraVector ( inVector->castTo<ManagedVector>().getRootVector() ) );
    }

    if ( !retVal)
    {
      AMP_ERROR( "Cannot create view!" );
    }

    return retVal;
  }

  Vector::shared_ptr  EpetraVector::createView ( Vector::shared_ptr  inVector )
  {
    DEPRECATED("createView","view");
    return view ( inVector );
  }

  const Vector::shared_ptr  EpetraVector::createConstView ( const Vector::shared_ptr  inVector )
  {
    DEPRECATED("createConstView","constView");
    return constView ( inVector );
  }

}
}

