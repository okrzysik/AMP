
#include "vectors/MultiVector.h"
#include "vectors/SimpleVector.h"

#include "ManagedPetscVector.h"


namespace AMP {
namespace LinearAlgebra {


void PetscVector::dataChanged ()
{
    PetscObjectStateIncrease ( reinterpret_cast< ::PetscObject> ( getVec() ) );
}


const Vector::shared_ptr  PetscVector::constView ( const Vector::shared_ptr inVector )
{
    Vector::shared_ptr  retVal;

    if ( inVector->isA<PetscVector> () )
        return inVector;
    if ( inVector->hasView<PetscVector> () )
        return inVector->getView<PetscVector>();
  
    if ( inVector->isA<ManagedVector> () ) {
        retVal = Vector::shared_ptr ( new ManagedPetscVector ( inVector ) );
        retVal->setVariable ( inVector->getVariable() );
        inVector->registerView ( retVal );
    } else if ( inVector->isA<VectorEngine> () ) {
        ManagedPetscVectorParameters *newParams = new ManagedPetscVectorParameters;
        newParams->d_Engine = boost::dynamic_pointer_cast<VectorEngine> ( inVector );
        newParams->d_CloneEngine = false;
        AMP_INSIST(inVector->getCommunicationList().get()!=NULL,"All vectors must have a communication list");
        newParams->d_CommList = inVector->getCommunicationList();
        AMP_INSIST(inVector->getDOFManager().get()!=NULL,"All vectors must have a DOFManager list");
        newParams->d_DOFManager = inVector->getDOFManager();
        ManagedPetscVector *t = new ManagedPetscVector ( VectorParameters::shared_ptr ( newParams ) );
        inVector->castTo<DataChangeFirer>().registerListener( t );
        t->setVariable ( inVector->getVariable() );
        t->setUpdateStatusPtr ( inVector->getUpdateStatusPtr () );
        retVal = Vector::shared_ptr ( t );
        inVector->registerView ( retVal );
    } else {
        AMP_ERROR( "Nobody uses constView, anyway" );
    }

    return retVal;
}


Vector::shared_ptr  PetscVector::view ( Vector::shared_ptr inVector )
{
    Vector::shared_ptr  retVal;

    if ( inVector->isA<PetscVector> () ) {
        retVal = inVector;
    } else if ( inVector->hasView<PetscVector> () ) {
        retVal = inVector->getView<PetscVector>();
    } else if ( inVector->isA<ManagedVector> () ) {
        retVal = Vector::shared_ptr ( new ManagedPetscVector ( inVector ) );
        inVector->registerView ( retVal );
    } else if ( inVector->isA<VectorEngine> () ) {
        ManagedPetscVectorParameters *newParams = new ManagedPetscVectorParameters;
        newParams->d_Engine = boost::dynamic_pointer_cast<VectorEngine> ( inVector );
        newParams->d_CloneEngine = false;
        AMP_INSIST(inVector->getCommunicationList().get()!=NULL,"All vectors must have a communication list");
        newParams->d_CommList = inVector->getCommunicationList();
        AMP_INSIST(inVector->getDOFManager().get()!=NULL,"All vectors must have a DOFManager list");
        newParams->d_DOFManager = inVector->getDOFManager();
        ManagedPetscVector *newVector = new ManagedPetscVector ( VectorParameters::shared_ptr ( newParams ) );
        inVector->castTo<DataChangeFirer>().registerListener( newVector );
        newVector->setVariable ( inVector->getVariable() );
        newVector->setUpdateStatusPtr ( inVector->getUpdateStatusPtr () );
        retVal = Vector::shared_ptr ( newVector );
        inVector->registerView ( retVal );
    } else if ( inVector->isA<SimpleVector> () ) {
        retVal = view ( MultiVector::view ( inVector , AMP_MPI(AMP_COMM_SELF) ) );  // This is an extraordinary hack so for SimpleVectors
        inVector->registerView ( retVal );
    } else {
        AMP_ERROR( "Failed view" );
    }

    return retVal;
}


}
}

