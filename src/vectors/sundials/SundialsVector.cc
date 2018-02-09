
#include "AMP/vectors/MultiVector.h"
#include "AMP/vectors/sundials/ManagedSundialsVector.h"


namespace AMP {
namespace LinearAlgebra {

Vector::const_shared_ptr SundialsVector::constView( Vector::const_shared_ptr inVector )
{
    Vector::shared_ptr retVal;
    if ( dynamic_pointer_cast<const SundialsVector>( inVector ) ) {
        return inVector;
    } else if ( inVector->hasView<SundialsVector>() ) {
        return inVector->getView<SundialsVector>();
    } else if ( dynamic_pointer_cast<const ManagedVector>( inVector ) ) {
        Vector::shared_ptr inVector2 = AMP::const_pointer_cast<Vector>( inVector );
        retVal                       = Vector::shared_ptr( new ManagedSundialsVector( inVector2 ) );
        inVector->registerView( retVal );
    } else if ( dynamic_pointer_cast<const VectorEngine>( inVector ) ) {
        Vector::shared_ptr inVector2 = AMP::const_pointer_cast<Vector>( inVector );
        auto new_params              = new ManagedSundialsVectorParameters;
        new_params->d_Engine         = AMP::dynamic_pointer_cast<VectorEngine>( inVector2 );
        new_params->d_CloneEngine    = false;
        if ( inVector->getCommunicationList().get() != nullptr )
            new_params->d_CommList = inVector->getCommunicationList();
        else
            new_params->d_CommList =
                CommunicationList::createEmpty( inVector->getLocalSize(), inVector->getComm() );
        if ( inVector->getDOFManager().get() != nullptr )
            new_params->d_DOFManager = inVector->getDOFManager();
        else
            new_params->d_DOFManager = AMP::make_shared<AMP::Discretization::DOFManager>(
                inVector->getLocalSize(), inVector->getComm() );
        ManagedSundialsVector *t =
            new ManagedSundialsVector( VectorParameters::shared_ptr( new_params ) );
        t->setVariable( inVector->getVariable() );
        t->setUpdateStatusPtr( inVector->getUpdateStatusPtr() );
        retVal = Vector::shared_ptr( t );
        inVector->registerView( retVal );
    } else {
        // Create a multivector to wrap the given vector and create a view
        Vector::shared_ptr inVector2 = AMP::const_pointer_cast<Vector>( inVector );
        retVal                       = view( MultiVector::view( inVector2, inVector->getComm() ) );
        inVector2->registerView( retVal );
    }
    return retVal;
}


Vector::shared_ptr SundialsVector::view( Vector::shared_ptr inVector )
{
    Vector::shared_ptr retVal;
    if ( dynamic_pointer_cast<SundialsVector>( inVector ) ) {
        retVal = inVector;
    } else if ( inVector->hasView<SundialsVector>() ) {
        retVal = inVector->getView<SundialsVector>();
    } else if ( dynamic_pointer_cast<ManagedVector>( inVector ) ) {
        retVal = Vector::shared_ptr( new ManagedSundialsVector( inVector ) );
        inVector->registerView( retVal );
    } else if ( dynamic_pointer_cast<VectorEngine>( inVector ) ) {
        auto new_params           = new ManagedSundialsVectorParameters;
        new_params->d_Engine      = AMP::dynamic_pointer_cast<VectorEngine>( inVector );
        new_params->d_CloneEngine = false;
        if ( inVector->getCommunicationList().get() != nullptr )
            new_params->d_CommList = inVector->getCommunicationList();
        else
            new_params->d_CommList =
                CommunicationList::createEmpty( inVector->getLocalSize(), inVector->getComm() );
        if ( inVector->getDOFManager().get() != nullptr )
            new_params->d_DOFManager = inVector->getDOFManager();
        else
            new_params->d_DOFManager = AMP::make_shared<AMP::Discretization::DOFManager>(
                inVector->getLocalSize(), inVector->getComm() );
        ManagedSundialsVector *t =
            new ManagedSundialsVector( VectorParameters::shared_ptr( new_params ) );
        t->setVariable( inVector->getVariable() );
        t->setUpdateStatusPtr( inVector->getUpdateStatusPtr() );
        retVal = Vector::shared_ptr( t );
        inVector->registerView( retVal );
    } else {
        // Create a multivector to wrap the given vector and create a view
        retVal = view( MultiVector::view( inVector, inVector->getComm() ) );
        inVector->registerView( retVal );
    }
    return retVal;
}


} // namespace LinearAlgebra
} // namespace AMP
