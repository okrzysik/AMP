
#include "AMP/vectors/MultiVector.h"
#include "AMP/vectors/sundials/ManagedSundialsVector.h"


namespace AMP {
namespace LinearAlgebra {


/****************************************************************
 * view                                                          *
 ****************************************************************/
Vector::const_shared_ptr SundialsVector::constView( Vector::const_shared_ptr inVector )
{
    return view( std::const_pointer_cast<Vector>( inVector ) );
}
Vector::shared_ptr SundialsVector::view( Vector::shared_ptr inVector )
{
    Vector::shared_ptr retVal;
    if ( std::dynamic_pointer_cast<SundialsVector>( inVector ) ) {
        retVal = inVector;
    } else if ( inVector->hasView<SundialsVector>() ) {
        retVal = inVector->getView<SundialsVector>();
    } else if ( std::dynamic_pointer_cast<ManagedVector>( inVector ) ) {
        retVal = std::make_shared<ManagedSundialsVector>( inVector );
        inVector->registerView( retVal );
    } else if ( std::dynamic_pointer_cast<MultiVector>( inVector ) ) {
        auto new_params      = std::make_shared<ManagedSundialsVectorParameters>();
        new_params->d_Engine = std::dynamic_pointer_cast<Vector>( inVector );
        if ( inVector->getCommunicationList() )
            new_params->d_CommList = inVector->getCommunicationList();
        else
            new_params->d_CommList =
                CommunicationList::createEmpty( inVector->getLocalSize(), inVector->getComm() );
        if ( inVector->getDOFManager() )
            new_params->d_DOFManager = inVector->getDOFManager();
        else
            new_params->d_DOFManager = std::make_shared<AMP::Discretization::DOFManager>(
                inVector->getLocalSize(), inVector->getComm() );
        auto t = std::make_shared<ManagedSundialsVector>( new_params );
        t->setVariable( inVector->getVariable() );
        t->getVectorData()->setUpdateStatusPtr( inVector->getVectorData()->getUpdateStatusPtr() );
        retVal = t;
        inVector->registerView( retVal );
    } else {
        // Create a multivector to wrap the given vector and create a view
        retVal = view( MultiVector::view( inVector, inVector->getComm() ) );
        inVector->registerView( retVal );
    }
    return retVal;
}

SundialsVector::SundialsVector() {}


N_Vector &SundialsVector::getNVector() { return d_n_vector; }

const N_Vector &SundialsVector::getNVector() const { return d_n_vector; }


} // namespace LinearAlgebra
} // namespace AMP
