#include "AMP/vectors/MultiVector.h"
#include "AMP/vectors/SimpleVector.h"

#include "ManagedPetscVector.h"

#include "petsc/private/petscimpl.h"

namespace AMP {
namespace LinearAlgebra {


void PetscVector::dataChanged()
{
    PetscObjectStateIncrease( reinterpret_cast<::PetscObject>( getVec() ) );
}


/****************************************************************
 * view                                                          *
 ****************************************************************/
Vector::const_shared_ptr PetscVector::constView( Vector::const_shared_ptr inVector )
{
    return view( std::const_pointer_cast<Vector>( inVector ) );
}
Vector::shared_ptr PetscVector::view( Vector::shared_ptr inVector )
{
    Vector::shared_ptr retVal;
    if ( std::dynamic_pointer_cast<PetscVector>( inVector ) ) {
        retVal = inVector;
    } else if ( inVector->hasView<PetscVector>() ) {
        retVal = inVector->getView<PetscVector>();
    } else if ( std::dynamic_pointer_cast<ManagedVector>( inVector ) ) {
        retVal = std::make_shared<ManagedPetscVector>( inVector );
        inVector->registerView( retVal );
    } else if ( std::dynamic_pointer_cast<MultiVector>( inVector ) ) {
        auto newParams      = std::make_shared<ManagedPetscVectorParameters>();
        newParams->d_Engine = std::dynamic_pointer_cast<Vector>( inVector );
        newParams->d_Buffer = std::dynamic_pointer_cast<VectorData>( inVector );
        AMP_INSIST( inVector->getCommunicationList(),
                    "All vectors must have a communication list" );
        newParams->d_CommList = inVector->getCommunicationList();
        AMP_INSIST( inVector->getDOFManager(), "All vectors must have a DOFManager list" );
        newParams->d_DOFManager = inVector->getDOFManager();
        auto newVector          = std::make_shared<ManagedPetscVector>( newParams );
        std::dynamic_pointer_cast<DataChangeFirer>( inVector )->registerListener( newVector.get() );
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

PetscVector::PetscVector() : d_PetscRandom( 0 ), d_petscVec( nullptr ) {}


PetscRandom &PetscVector::getPetscRandom( const AMP_MPI &comm )
{
    if ( d_PetscRandom == 0 ) {
        d_PetscRandom = new PetscRandom;
        PetscRandomCreate( comm.getCommunicator(), d_PetscRandom );
        PetscRandomSetType( *d_PetscRandom, PETSCRAND48 ); // This is a horrible RNG for
                                                           // stochastic simulation.  Do not
                                                           // use.
    }
    return *d_PetscRandom;
}


Vec &PetscVector::getVec() { return d_petscVec; }


const Vec &PetscVector::getVec() const { return d_petscVec; }


PetscVector::~PetscVector()
{
    if ( d_PetscRandom ) {
        PETSC::randomDestroy( d_PetscRandom );
        delete d_PetscRandom;
    }
}


} // namespace LinearAlgebra
} // namespace AMP
