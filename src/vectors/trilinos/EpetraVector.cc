#include "EpetraVector.h"
#include "ManagedEpetraVector.h"
#include "MultiVector.h"


namespace AMP {
namespace LinearAlgebra {


  EpetraVector::EpetraVector()
  {
  }

  EpetraVector::~EpetraVector ()
  {
  }

  const Vector::shared_ptr  EpetraVector::constView ( const Vector::shared_ptr inVector )
  {
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

    if ( inVector->isA<EpetraVector> () ) {
        retVal = inVector;
    } else if ( inVector->isA<MultiVector> () ) {
        if ( inVector->numberOfDataBlocks() == 1 ) {
            Vector::shared_ptr localVector = inVector->castTo<MultiVector>().getVector( 0 );
            retVal = view ( localVector );
        } else {
            AMP_ERROR("View of multi-block MultiVector is not supported yet");
        }
    } else if ( inVector->isA<ManagedVector> () ) {
        boost::shared_ptr<ManagedEpetraVector> managed( new ManagedEpetraVector ( inVector->castTo<ManagedVector>().getRootVector() ) );
        retVal = managed;
    }

    if ( !retVal)
        AMP_ERROR( "Cannot create view!" );

    return retVal;
  }


}
}

