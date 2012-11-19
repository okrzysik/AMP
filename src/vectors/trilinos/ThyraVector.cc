#include "vectors/trilinos/ThyraVector.h"
#include "vectors/trilinos/ManagedThyraVector.h"
#include "vectors/SimpleVector.h"


namespace AMP {
namespace LinearAlgebra {


/****************************************************************
* constView                                                     *
****************************************************************/
Vector::const_shared_ptr  ThyraVector::constView ( Vector::const_shared_ptr inVector )
{
    // Check if we have an exisiting view
    if ( boost::dynamic_pointer_cast<const ThyraVector>(inVector)!=NULL )
        return inVector;
    if ( inVector->hasView<ManagedThyraVector> () )
        return inVector->getView<ManagedThyraVector>();
    // Create a new view
    Vector::shared_ptr  retVal;
    if ( inVector->isA<ManagedVector> () ) {
        Vector::shared_ptr inVector2 = boost::const_pointer_cast<Vector>( inVector );
        retVal = Vector::shared_ptr ( new ManagedThyraVector( inVector2 ) );
        retVal->setVariable ( inVector->getVariable() );
        inVector->registerView ( retVal );
    } else if ( inVector->isA<VectorEngine> () ) {
        Vector::shared_ptr inVector2 = boost::const_pointer_cast<Vector>( inVector );
        ManagedThyraVectorParameters *newParams = new ManagedThyraVectorParameters;
        newParams->d_Engine = boost::dynamic_pointer_cast<VectorEngine>( inVector2 );
        newParams->d_CloneEngine = false;
        AMP_INSIST(inVector->getCommunicationList().get()!=NULL,"All vectors must have a communication list");
        newParams->d_CommList = inVector->getCommunicationList();
        AMP_INSIST(inVector->getDOFManager().get()!=NULL,"All vectors must have a DOFManager list");
        newParams->d_DOFManager = inVector->getDOFManager();
        ManagedThyraVector *t = new ManagedThyraVector ( VectorParameters::shared_ptr ( newParams ) );
        retVal = Vector::shared_ptr ( t );
        inVector->registerView ( retVal );
    } else {
        Vector::shared_ptr inVector2 = boost::const_pointer_cast<Vector>( inVector );
        retVal = view ( MultiVector::view ( inVector2, inVector->getComm() ) );
        inVector->registerView ( retVal );
    }
    return retVal;
}


/****************************************************************
* View                                                          *
****************************************************************/
Vector::shared_ptr  ThyraVector::view ( Vector::shared_ptr inVector )
{
    // Check if we have an exisiting view
    if ( boost::dynamic_pointer_cast<ThyraVector>(inVector)!=NULL )
        return inVector;
    if ( inVector->hasView<ManagedThyraVector> () )
        return inVector->getView<ManagedThyraVector>();
    // Create a new view
    Vector::shared_ptr  retVal;
    if ( inVector->isA<ManagedVector> () ) {
        retVal = Vector::shared_ptr ( new ManagedThyraVector ( inVector ) );
        inVector->registerView ( retVal );
    } else if ( inVector->isA<VectorEngine> () ) {
        ManagedThyraVectorParameters *newParams = new ManagedThyraVectorParameters;
        newParams->d_Engine = boost::dynamic_pointer_cast<VectorEngine> ( inVector );
        newParams->d_CloneEngine = false;
        AMP_INSIST(inVector->getCommunicationList().get()!=NULL,"All vectors must have a communication list");
        newParams->d_CommList = inVector->getCommunicationList();
        AMP_INSIST(inVector->getDOFManager().get()!=NULL,"All vectors must have a DOFManager list");
        newParams->d_DOFManager = inVector->getDOFManager();
        ManagedThyraVector *newVector = new ManagedThyraVector ( VectorParameters::shared_ptr ( newParams ) );
        newVector->setVariable ( inVector->getVariable() );
        newVector->setUpdateStatusPtr ( inVector->getUpdateStatusPtr () );
        retVal = Vector::shared_ptr ( newVector );
        inVector->registerView ( retVal );
    } else {
        retVal = view ( MultiVector::view ( inVector, inVector->getComm() ) );
        inVector->registerView ( retVal );
    }
    return retVal;
}


/****************************************************************
* Return the thyra vector                                       *
****************************************************************/
Teuchos::RCP<Thyra::VectorBase<double> > ThyraVector::getVec()
{
    return d_thyraVec;
}
Teuchos::RCP<const Thyra::VectorBase<double> >  ThyraVector::getVec() const
{
    return d_thyraVec;
}


}
}

