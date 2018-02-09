
#include "AMP/vectors/MultiVector.h"
#include "AMP/vectors/SimpleVector.h"

#include "ManagedPetscVector.h"


namespace AMP {
namespace LinearAlgebra {


void PetscVector::dataChanged()
{
    PetscObjectStateIncrease( reinterpret_cast<::PetscObject>( getVec() ) );
}


Vector::const_shared_ptr PetscVector::constView( Vector::const_shared_ptr inVector )
{
    Vector::shared_ptr retVal;
    if ( dynamic_pointer_cast<const PetscVector>( inVector ) ) {
        return inVector;
    } else if ( inVector->hasView<PetscVector>() ) {
        return inVector->getView<PetscVector>();
    } else if ( dynamic_pointer_cast<const ManagedVector>( inVector ) ) {
        auto inVector2 = AMP::const_pointer_cast<Vector>( inVector );
        retVal         = AMP::make_shared<ManagedPetscVector>( inVector2 );
        retVal->setVariable( inVector->getVariable() );
        inVector->registerView( retVal );
    } else if ( dynamic_pointer_cast<const VectorEngine>( inVector ) ) {
        auto inVector2           = AMP::const_pointer_cast<Vector>( inVector );
        auto newParams           = AMP::make_shared<ManagedPetscVectorParameters>();
        newParams->d_Engine      = AMP::dynamic_pointer_cast<VectorEngine>( inVector2 );
        newParams->d_CloneEngine = false;
        AMP_INSIST( inVector->getCommunicationList(),
                    "All vectors must have a communication list" );
        newParams->d_CommList = inVector->getCommunicationList();
        AMP_INSIST( inVector->getDOFManager(), "All vectors must have a DOFManager list" );
        newParams->d_DOFManager = inVector->getDOFManager();
        auto newVector          = AMP::make_shared<ManagedPetscVector>( newParams );
        dynamic_pointer_cast<DataChangeFirer>( inVector2 )->registerListener( newVector.get() );
        newVector->setVariable( inVector->getVariable() );
        newVector->setUpdateStatusPtr( inVector->getUpdateStatusPtr() );
        inVector->registerView( retVal );
        retVal = newVector;
    } else {
        auto inVector2 = AMP::const_pointer_cast<Vector>( inVector );
        retVal         = view( MultiVector::view( inVector2, inVector->getComm() ) );
        inVector->registerView( retVal );
    }
    return retVal;
}


Vector::shared_ptr PetscVector::view( Vector::shared_ptr inVector )
{
    Vector::shared_ptr retVal;
    if ( dynamic_pointer_cast<PetscVector>( inVector ) ) {
        retVal = inVector;
    } else if ( inVector->hasView<PetscVector>() ) {
        retVal = inVector->getView<PetscVector>();
    } else if ( dynamic_pointer_cast<ManagedVector>( inVector ) ) {
        retVal = Vector::shared_ptr( new ManagedPetscVector( inVector ) );
        inVector->registerView( retVal );
    } else if ( dynamic_pointer_cast<VectorEngine>( inVector ) ) {
        auto newParams           = AMP::make_shared<ManagedPetscVectorParameters>();
        newParams->d_Engine      = AMP::dynamic_pointer_cast<VectorEngine>( inVector );
        newParams->d_CloneEngine = false;
        AMP_INSIST( inVector->getCommunicationList(),
                    "All vectors must have a communication list" );
        newParams->d_CommList = inVector->getCommunicationList();
        AMP_INSIST( inVector->getDOFManager(), "All vectors must have a DOFManager list" );
        newParams->d_DOFManager = inVector->getDOFManager();
        auto newVector          = AMP::make_shared<ManagedPetscVector>( newParams );
        dynamic_pointer_cast<DataChangeFirer>( inVector )->registerListener( newVector.get() );
        newVector->setVariable( inVector->getVariable() );
        newVector->setUpdateStatusPtr( inVector->getUpdateStatusPtr() );
        inVector->registerView( newVector );
        retVal = newVector;
    } else {
        // Create a multivector to wrap the given vector and create a view
        retVal = view( MultiVector::view( inVector, inVector->getComm() ) );
        inVector->registerView( retVal );
    }
    return retVal;
}


} // namespace LinearAlgebra
} // namespace AMP
