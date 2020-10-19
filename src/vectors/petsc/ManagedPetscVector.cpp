#include "AMP/vectors/petsc/ManagedPetscVector.h"
#include "AMP/utils/Utilities.h"
#include "AMP/vectors/data/ManagedVectorData.h"
#include "AMP/vectors/data/VectorDataCPU.h"
#include "AMP/vectors/data/VectorDataIterator.h"
#include "AMP/vectors/operations/ManagedVectorOperations.h"
#include "AMP/vectors/petsc/PetscHelpers.h"
#include "AMP/vectors/petsc/PetscVector.h"

#include "petsc/private/vecimpl.h"
#include "petscsys.h"


namespace AMP {
namespace LinearAlgebra {


static inline auto getVectorEngine( const std::shared_ptr<VectorData> &data )
{
    auto managed = std::dynamic_pointer_cast<ManagedVectorData>( data );
    AMP_ASSERT( managed );
    return managed->getVectorEngine();
}
static inline auto getVectorEngine( const std::shared_ptr<const VectorData> &data )
{
    auto managed = std::dynamic_pointer_cast<const ManagedVectorData>( data );
    AMP_ASSERT( managed );
    return managed->getVectorEngine();
}


ManagedPetscVector::ManagedPetscVector( Vector::shared_ptr vec ) : PetscVector()
{
    AMP_ASSERT( !std::dynamic_pointer_cast<ManagedVectorData>( vec->getVectorData() ) );
    d_VectorOps  = std::make_shared<ManagedVectorOperations>();
    d_VectorData = std::make_shared<ManagedVectorData>( vec );
    d_DOFManager = vec->getDOFManager();
    setVariable( vec->getVariable() );
    d_wrapper     = std::make_shared<PETSC::PetscVectorWrapper>( this );
    auto listener = std::dynamic_pointer_cast<DataChangeListener>( shared_from_this() );
    d_VectorData->registerListener( listener );
}


ManagedPetscVector::~ManagedPetscVector() {}


bool ManagedPetscVector::petscHoldsView() const { return d_wrapper->petscHoldsView(); }


ManagedPetscVector *ManagedPetscVector::petscDuplicate()
{
    ManagedPetscVector *pAns = rawClone();
    pAns->setVariable( getVariable() );
    pAns->d_wrapper->setMadeWithPetscDuplicate( true );
    return pAns;
}


void ManagedPetscVector::swapVectors( Vector &other )
{
    auto tmp = dynamic_cast<ManagedPetscVector *>( &other );
    AMP_ASSERT( tmp != nullptr );
    d_VectorData->swapData( *other.getVectorData() );
}

bool ManagedPetscVector::constructedWithPetscDuplicate()
{
    return d_wrapper->constructedWithPetscDuplicate();
}

ManagedPetscVector *ManagedPetscVector::rawClone() const
{
    auto vec    = getVectorEngine( getVectorData() );
    auto vec2   = vec->cloneVector( vec->getVariable() );
    auto retVal = new ManagedPetscVector( vec2 );
    return retVal;
}

Vector::shared_ptr ManagedPetscVector::cloneVector( const Variable::shared_ptr p ) const
{
    Vector::shared_ptr retVal( rawClone() );
    retVal->setVariable( p );
    return retVal;
}


void ManagedPetscVector::receiveDataChanged()
{
    PetscObjectStateIncrease( reinterpret_cast<::PetscObject>( getVec() ) );
}

bool ManagedPetscVector::isAnAliasOf( const ManagedPetscVector &rhs ) const
{
    auto managedData = std::dynamic_pointer_cast<const ManagedVectorData>( getVectorData() );
    return managedData->isAnAliasOf( *( rhs.getVectorData() ) );
}


/********************************************************
 * Subset                                                *
 ********************************************************/
Vector::shared_ptr ManagedPetscVector::subsetVectorForVariable( Variable::const_shared_ptr name )
{
    Vector::shared_ptr retVal;
    if ( !retVal )
        retVal = Vector::subsetVectorForVariable( name );
    if ( !retVal ) {
        auto vec = getVectorEngine( getVectorData() );
        if ( vec )
            retVal = vec->subsetVectorForVariable( name );
    }
    return retVal;
}
Vector::const_shared_ptr
ManagedPetscVector::constSubsetVectorForVariable( Variable::const_shared_ptr name ) const
{
    Vector::const_shared_ptr retVal;
    if ( !retVal )
        retVal = Vector::constSubsetVectorForVariable( name );
    if ( !retVal ) {
        auto const vec = getVectorEngine( getVectorData() );
        if ( vec )
            retVal = vec->constSubsetVectorForVariable( name );
    }
    if ( !retVal ) {
        auto const vec = getVectorEngine( getVectorData() );
        printf( "Unable to subset for %s in %s:%s\n",
                name->getName().data(),
                getVariable()->getName().data(),
                vec->getVariable()->getName().data() );
    }
    return retVal;
}


} // namespace LinearAlgebra
} // namespace AMP
