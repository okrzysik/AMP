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


void ManagedPetscVector::initPetsc()
{
    AMP_MPI comm = getVectorEngine( getVectorData() )->getComm();
    VecCreate( comm.getCommunicator(), &d_petscVec );

    d_petscVec->data        = this;
    d_petscVec->petscnative = PETSC_FALSE;

    PETSC::reset_vec_ops( d_petscVec );

#if ( PETSC_VERSION_MAJOR == 3 && PETSC_VERSION_MINOR == 0 )
    PetscMapInitialize( comm.getCommunicator(), d_petscVec->map );
    PetscMapSetBlockSize( d_petscVec->map, 1 );
    PetscMapSetSize( d_petscVec->map, this->getGlobalSize() );
    PetscMapSetLocalSize( d_petscVec->map, this->getLocalSize() );
    d_petscVec->map->rstart = static_cast<PetscInt>( this->getDOFManager()->beginDOF() );
    d_petscVec->map->rend   = static_cast<PetscInt>( this->getDOFManager()->endDOF() );
#elif PETSC_VERSION_GE( 3, 2, 0 )
    PetscLayoutSetBlockSize( d_petscVec->map, 1 );
    PetscLayoutSetSize( d_petscVec->map, this->getGlobalSize() );
    PetscLayoutSetLocalSize( d_petscVec->map, this->getLocalSize() );
    PetscLayoutSetUp( d_petscVec->map );
#else
#error Not programmed for this version yet
#endif

    d_bMadeWithPetscDuplicate = false;

    const std::string my_name = "AMPManagedPetscVectorReal";
    int ierr                  = 0;

    if ( ( (PetscObject) d_petscVec )->type_name ) {
        ierr = PetscFree( ( (PetscObject) d_petscVec )->type_name );
    }

    ierr =
        PetscObjectChangeTypeName( reinterpret_cast<PetscObject>( d_petscVec ), my_name.c_str() );
    AMP_INSIST( ierr == 0, "PetscObjectChangeTypeName returned non-zero error code" );

    VecSetFromOptions( d_petscVec );
}


ManagedPetscVector::ManagedPetscVector( Vector::shared_ptr vec ) : PetscVector()
{
    AMP_ASSERT( !std::dynamic_pointer_cast<ManagedVectorData>( vec->getVectorData() ) );
    d_VectorOps  = std::make_shared<ManagedVectorOperations>();
    d_VectorData = std::make_shared<ManagedVectorData>( vec );
    d_DOFManager = vec->getDOFManager();
    setVariable( vec->getVariable() );
    initPetsc();
    auto listener = std::dynamic_pointer_cast<DataChangeListener>( shared_from_this() );
    d_VectorData->registerListener( listener );
}


ManagedPetscVector::~ManagedPetscVector()
{
    int refct = ( ( (PetscObject) d_petscVec )->refct );
    if ( !d_bMadeWithPetscDuplicate ) {
        if ( refct > 1 )
            AMP_ERROR( "Deleting a vector still held by PETSc" );
        PETSC::vecDestroy( &d_petscVec );
    }
}


bool ManagedPetscVector::petscHoldsView() const
{
    int refct = ( ( (PetscObject) d_petscVec )->refct );
    if ( !d_bMadeWithPetscDuplicate && refct > 1 )
        return true;
    return false;
}


ManagedPetscVector *ManagedPetscVector::petscDuplicate()
{
    ManagedPetscVector *pAns = rawClone();
    pAns->setVariable( getVariable() );
    pAns->d_bMadeWithPetscDuplicate = true;
    return pAns;
}


void ManagedPetscVector::swapVectors( Vector &other )
{
    auto tmp = dynamic_cast<ManagedPetscVector *>( &other );
    AMP_ASSERT( tmp != nullptr );
    d_VectorData->swapData( *other.getVectorData() );
}

bool ManagedPetscVector::constructedWithPetscDuplicate() { return d_bMadeWithPetscDuplicate; }

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
