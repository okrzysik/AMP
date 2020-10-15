#include "AMP/vectors/MultiVector.h"
#include "AMP/vectors/data/ManagedVectorData.h"
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
    // Check if we have an existing view
    if ( std::dynamic_pointer_cast<SundialsVector>( inVector ) )
        return inVector;
    if ( inVector->hasView<SundialsVector>() )
        return inVector->getView<SundialsVector>();
    // Check if we are dealing with a managed vector
    auto managedData = std::dynamic_pointer_cast<ManagedVectorData>( inVector->getVectorData() );
    if ( managedData ) {
        auto retVal = view( managedData->getVectorEngine() );
        retVal->setVariable( inVector->getVariable() );
        return retVal;
    }
    // Check if we are dealing with a multivector
    if ( std::dynamic_pointer_cast<MultiVector>( inVector ) ) {
        auto retVal = std::make_shared<ManagedSundialsVector>( inVector );
        retVal->setVariable( inVector->getVariable() );
        retVal->getVectorData()->setUpdateStatusPtr(
            inVector->getVectorData()->getUpdateStatusPtr() );
        inVector->registerView( retVal );
        return retVal;
    }
    // Create a multivector to wrap the given vector and create a view
    auto retVal = view( MultiVector::view( inVector, inVector->getComm() ) );
    retVal->setVariable( inVector->getVariable() );
    inVector->registerView( retVal );
    return retVal;
}

SundialsVector::SundialsVector() {}


N_Vector &SundialsVector::getNVector() { return d_n_vector; }

const N_Vector &SundialsVector::getNVector() const { return d_n_vector; }


} // namespace LinearAlgebra
} // namespace AMP
