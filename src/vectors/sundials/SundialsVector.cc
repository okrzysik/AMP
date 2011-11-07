#include "vectors/sundials/ManagedSundialsVector.h"


namespace AMP {
namespace LinearAlgebra {

  const Vector::shared_ptr  SundialsVector::constView ( const Vector::shared_ptr inVector )
  {
    Vector::shared_ptr  retVal;

    if ( inVector->isA<SundialsVector> () )
    {
      return inVector;
    }

    if ( inVector->hasView<SundialsVector> () )
    {
      return inVector->getView<SundialsVector>();
    }


    if ( inVector->isA<ManagedVector> () )
    {
      retVal = Vector::shared_ptr ( new ManagedSundialsVector ( inVector ) );
      inVector->registerView ( retVal );
    }
    else if ( inVector->isA<VectorEngine> () )
    {
      ManagedSundialsVectorParameters *new_params = new ManagedSundialsVectorParameters;
      new_params->d_Engine = boost::dynamic_pointer_cast<VectorEngine> ( inVector );
      new_params->d_CloneEngine = false;
      new_params->d_CommList = inVector->getCommunicationList() ? inVector->getCommunicationList()
                                                                : CommunicationList::createEmpty ( inVector->getLocalSize() );
      ManagedSundialsVector *t = new ManagedSundialsVector ( VectorParameters::shared_ptr ( new_params ) );
      t->setVariable ( inVector->getVariable() );
      t->setUpdateStatus ( inVector->getUpdateStatus () );
      retVal = Vector::shared_ptr ( t );
      inVector->registerView ( retVal );
    }
    else
    {
      AMP_ERROR( "Cannot create view!" );
    }

    return retVal;
  }

  Vector::shared_ptr  SundialsVector::view ( Vector::shared_ptr inVector )
  {
    Vector::shared_ptr  retVal;

    if ( inVector->isA<SundialsVector> () )
    {
      return inVector;
    }

    if ( inVector->hasView<SundialsVector> () )
    {
      return inVector->getView<SundialsVector>();
    }

    if ( inVector->isA<ManagedVector> () )
    {
      retVal = Vector::shared_ptr ( new ManagedSundialsVector ( inVector ) );
      inVector->registerView ( retVal );
    }
    else if ( inVector->isA<VectorEngine> () )
    {
      ManagedSundialsVectorParameters *new_params = new ManagedSundialsVectorParameters;
      new_params->d_Engine = boost::dynamic_pointer_cast<VectorEngine> ( inVector );
      new_params->d_CloneEngine = false;
      new_params->d_CommList = inVector->getCommunicationList() ? inVector->getCommunicationList()
                                                                : CommunicationList::createEmpty ( inVector->getLocalSize() );
      ManagedSundialsVector *t = new ManagedSundialsVector ( VectorParameters::shared_ptr ( new_params ) );
      t->setVariable ( inVector->getVariable() );
      t->setUpdateStatus ( inVector->getUpdateStatus () );
      retVal = Vector::shared_ptr ( t );
      inVector->registerView ( retVal );
    }
    else
    {
      AMP_ERROR( "Cannot create view!" );
    }

    return retVal;
  }

}
}

