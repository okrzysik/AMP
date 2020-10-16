#include "AMP/vectors/MultiVector.h"
#include "AMP/vectors/data/ManagedVectorData.h"
#include "AMP/vectors/petsc/ManagedPetscVector.h"

#include "petsc/private/petscimpl.h"


namespace AMP {
namespace LinearAlgebra {


/****************************************************************
 * view                                                          *
 ****************************************************************/
Vector::const_shared_ptr PetscVector::constView( Vector::const_shared_ptr inVector )
{
    return view( std::const_pointer_cast<Vector>( inVector ) );
}
Vector::shared_ptr PetscVector::view( Vector::shared_ptr inVector )
{
    // Check if we have an existing view
    if ( std::dynamic_pointer_cast<PetscVector>( inVector ) )
        return inVector;
    if ( inVector->hasView<PetscVector>() )
        return inVector->getView<PetscVector>();
    // Check if we are dealing with a managed vector
    auto managedData = std::dynamic_pointer_cast<ManagedVectorData>( inVector->getVectorData() );
    if ( managedData ) {
        auto retVal = view( managedData->getVectorEngine() );
        retVal->setVariable( inVector->getVariable() );
        return retVal;
    }
    // Check if we are dealing with a multivector
    if ( std::dynamic_pointer_cast<MultiVector>( inVector ) ) {
        auto newVector = std::make_shared<ManagedPetscVector>( inVector );
        newVector->setVariable( inVector->getVariable() );
        newVector->getVectorData()->setUpdateStatusPtr(
            inVector->getVectorData()->getUpdateStatusPtr() );
        inVector->registerView( newVector );
        return newVector;
    }
    // Create a multivector to wrap the given vector and create a view
    // Note: this is required so that we call the native vector's operations
    auto newVector = view( MultiVector::view( inVector, inVector->getComm() ) );
    inVector->registerView( newVector );
    return newVector;
}

PetscVector::PetscVector() : d_petscVec( nullptr ) {}


Vec &PetscVector::getVec() { return d_petscVec; }


const Vec &PetscVector::getVec() const { return d_petscVec; }


PetscVector::~PetscVector() {}


} // namespace LinearAlgebra
} // namespace AMP
