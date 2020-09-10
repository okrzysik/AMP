#include "AMP/vectors/petsc/ManagedPetscVector.h"
#include "AMP/utils/Utilities.h"
#include "AMP/vectors/data/VectorDataCPU.h"
#include "AMP/vectors/data/VectorDataIterator.h"
#include "AMP/vectors/petsc/PetscHelpers.h"
#include "AMP/vectors/petsc/PetscVector.h"
#ifdef USE_EXT_TRILINOS
#include "AMP/vectors/trilinos/epetra/EpetraVectorEngine.h"
#endif


#include "petsc/private/vecimpl.h"
#include "petscsys.h"


namespace AMP {
namespace LinearAlgebra {


void ManagedPetscVector::initPetsc()
{
    AMP_MPI comm = getVectorEngine()->getComm();
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


ManagedPetscVector::ManagedPetscVector( VectorParameters::shared_ptr params )
    : ManagedVector( params ), PetscVector()
{
    initPetsc();
    auto listener = std::dynamic_pointer_cast<DataChangeListener>( shared_from_this() );
    d_VectorData->registerListener( listener );
}


ManagedPetscVector::ManagedPetscVector( Vector::shared_ptr alias )
    : ManagedVector( alias ), PetscVector()
{
    initPetsc();
    auto listener = std::dynamic_pointer_cast<DataChangeListener>( shared_from_this() );
    //    alias->registerListener( listener );
    alias->getVectorData()->registerListener( listener );
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


void ManagedPetscVector::copyFromPetscVec( Vector &dest, Vec source )
{
    if ( sizeof( PetscInt ) < 8 )
        AMP_INSIST( dest.getGlobalSize() < 0x80000000,
                    "PETsc is compiled with 32-bit integers and "
                    "we are trying to use a vector with more "
                    "than 2^31 elements" );

    auto ids       = new PetscInt[dest.getLocalSize()];
    PetscInt begin = dest.getLocalStartID();
    PetscInt end   = begin + dest.getLocalSize() - 1;

    for ( PetscInt i = begin; i < end; i++ )
        ids[i - begin] = i;
    VecGetValues( source, dest.getLocalSize(), ids, dest.getRawDataBlock<double>() );
    delete[] ids;
}


std::shared_ptr<AMP::LinearAlgebra::Vector> ManagedPetscVector::createFromPetscVec( Vec source,
                                                                                    AMP_MPI &comm )
{
#ifdef USE_EXT_TRILINOS
    PetscInt local_size, global_size, local_start, local_end;
    VecGetLocalSize( source, &local_size );
    VecGetSize( source, &global_size );
    VecGetOwnershipRange( source, &local_start, &local_end );
    auto buffer = std::make_shared<VectorDataCPU<double>>( local_start, local_size, global_size );
    auto t      = std::make_shared<ManagedPetscVectorParameters>();
    auto ve_params =
        std::make_shared<EpetraVectorEngineParameters>( local_size, global_size, comm );
    t->d_Engine = std::make_shared<EpetraVectorEngine>( ve_params, buffer );
    auto pRetVal =
        std::make_shared<ManagedPetscVector>( std::dynamic_pointer_cast<VectorParameters>( t ) );
    return pRetVal;
#else
    AMP_ERROR( "General case not programmed yet" );
    NULL_USE( source );
    NULL_USE( comm );
    return std::shared_ptr<AMP::LinearAlgebra::Vector>();
#endif
}


void ManagedPetscVector::swapVectors( Vector &other )
{
    auto tmp = dynamic_cast<ManagedPetscVector *>( &other );
    AMP_ASSERT( tmp != nullptr );
    ParentVector::swapVectors( *tmp );
}

ManagedVector *ManagedPetscVector::getNewRawPtr() const
{
    return new ManagedPetscVector( std::dynamic_pointer_cast<VectorParameters>( d_pParameters ) );
}

bool ManagedPetscVector::constructedWithPetscDuplicate() { return d_bMadeWithPetscDuplicate; }

ManagedPetscVector *ManagedPetscVector::rawClone() const
{
    auto p   = std::make_shared<ManagedPetscVectorParameters>();
    auto vec = getVectorEngine();
    if ( vec ) {
        auto vec2   = vec->cloneVector( "ManagedPetscVectorClone" );
        p->d_Buffer = std::dynamic_pointer_cast<VectorData>( vec2 );
        p->d_Engine = std::dynamic_pointer_cast<Vector>( vec2 );
    } else {
        AMP_ERROR( "ManagedPetscVector::rawClone() should not have reached here!" );
    }
    p->d_CommList   = getCommunicationList();
    p->d_DOFManager = getDOFManager();
    auto retVal     = new ManagedPetscVector( p );
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

} // namespace LinearAlgebra
} // namespace AMP
