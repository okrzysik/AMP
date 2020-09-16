#include "AMP/vectors/petsc/NativePetscVectorData.h"
#include "AMP/vectors/petsc/NativePetscVectorOperations.h"

#include "AMP/vectors/petsc/ManagedPetscVector.h"

#include "petsc.h"
#include "petsc/private/vecimpl.h"
#include "petscvec.h"

namespace AMP {
namespace LinearAlgebra {

NativePetscVectorData::NativePetscVectorData( Vec v, bool deleteable, AMP_MPI comm )
    : VectorData(), PetscVector()
{
    // Set the vector
    d_petscVec  = v;
    d_pArray    = nullptr;
    d_bDeleteMe = deleteable;
    // Get the correct communicator if it is not set
    MPI_Comm comm2 = comm.getCommunicator(); // Get a MPI_comm object from AMP_MPI to pass to PETSc
    PetscObjectGetComm( reinterpret_cast<PetscObject>( v ), &comm2 );
    if ( comm2 != comm.getCommunicator() )
        comm = AMP_MPI( comm2 );
    // Create the communication list
    auto params         = std::make_shared<CommunicationListParameters>();
    params->d_comm      = comm;
    params->d_localsize = getLocalSize();
    setCommunicationList( std::make_shared<CommunicationList>( params ) );
}


NativePetscVectorData::~NativePetscVectorData()
{
    resetArray();
    if ( d_bDeleteMe )
        PETSC::vecDestroy( &d_petscVec );
}

void NativePetscVectorData::putRawData( const double *in )
{
    int a, b;
    VecGetOwnershipRange( d_petscVec, &a, &b );
    AMP_ASSERT( b - a == (int) getLocalSize() );
    std::vector<int> offs( b - a );
    for ( size_t j = 0; j != offs.size(); j++ )
        offs[j] = a + j;
    VecSetValues( d_petscVec, offs.size(), offs.data(), in, INSERT_VALUES );
}


void NativePetscVectorData::copyOutRawData( double *out ) const
{
    std::copy( getRawDataBlock<double>( 0 ), getRawDataBlock<double>( 0 ) + getLocalSize(), out );
}


void NativePetscVectorData::swapData( VectorData &other )
{
    resetArray();
    auto otherData = dynamic_cast<NativePetscVectorData *>( &other );
    otherData->resetArray();
    VecSwap( d_petscVec, otherData->getVec() );
}

std::shared_ptr<VectorData> NativePetscVectorData::cloneData() const
{
    resetArray();
    Vec new_petscVec;
    VecDuplicate( d_petscVec, &new_petscVec );
    return std::make_shared<NativePetscVectorData>( new_petscVec, true, getComm() );
}

size_t NativePetscVectorData::numberOfDataBlocks() const { return 1; }


size_t NativePetscVectorData::sizeOfDataBlock( size_t i ) const
{
    if ( i != 0 )
        return 0;
    return getLocalSize();
}

void NativePetscVectorData::resetArray()
{
    if ( d_pArray ) {
        VecRestoreArray( d_petscVec, &d_pArray );
        d_pArray = 0;
    }
}

void NativePetscVectorData::resetArray() const
{
    if ( d_pArray ) {
        VecRestoreArray( d_petscVec, &d_pArray );
        d_pArray = 0;
    }
}

void NativePetscVectorData::setValuesByLocalID( int num, size_t *indices, const double *vals )
{
    for ( int i = 0; i != num; i++ )
        getRawDataBlock<double>()[indices[i]] = vals[i];
}


void NativePetscVectorData::setLocalValuesByGlobalID( int num, size_t *indices, const double *vals )
{
    resetArray();
    if ( sizeof( size_t ) == sizeof( PetscInt ) ) {
        VecSetValues( d_petscVec, num, (PetscInt *) indices, vals, INSERT_VALUES );
    } else {
        AMP_ASSERT( getGlobalSize() < 0x80000000 );
        std::vector<PetscInt> indices2( num, 0 );
        for ( int i = 0; i < num; i++ )
            indices2[i] = (PetscInt) indices[i];
        VecSetValues( d_petscVec, num, &indices2[0], vals, INSERT_VALUES );
    }
}


void NativePetscVectorData::addValuesByLocalID( int num, size_t *indices, const double *vals )
{
    for ( int i = 0; i != num; i++ )
        getRawDataBlock<double>()[indices[i]] += vals[i];
}


void NativePetscVectorData::addLocalValuesByGlobalID( int num, size_t *indices, const double *vals )
{
    resetArray();
    if ( sizeof( size_t ) == sizeof( PetscInt ) ) {
        VecSetValues( d_petscVec, num, (PetscInt *) indices, vals, ::ADD_VALUES );
    } else {
        AMP_ASSERT( getGlobalSize() < 0x80000000 );
        std::vector<PetscInt> indices2( num, 0 );
        for ( int i = 0; i < num; i++ )
            indices2[i] = (PetscInt) indices[i];
        VecSetValues( d_petscVec, num, &indices2[0], vals, ::ADD_VALUES );
    }
}

size_t NativePetscVectorData::getLocalSize() const
{
    int lsize;
    VecGetLocalSize( d_petscVec, &lsize );
    return static_cast<size_t>( lsize );
}

size_t NativePetscVectorData::getGlobalSize() const
{
    int gsize;
    VecGetSize( d_petscVec, &gsize );
    return static_cast<size_t>( gsize );
}


void *NativePetscVectorData::getRawDataBlockAsVoid( size_t i )
{
    if ( i > 0 )
        return 0;
    if ( d_pArray == 0 ) {
        VecGetArray( d_petscVec, &d_pArray );
    }
    return d_pArray;
}

const void *NativePetscVectorData::getRawDataBlockAsVoid( size_t i ) const
{
    if ( i > 0 )
        return 0;
    if ( d_pArray == 0 ) {
        VecGetArray( d_petscVec, &d_pArray );
    }
    return d_pArray;
}


void NativePetscVectorData::getLocalValuesByGlobalID( int numVals, size_t *ndx, double *vals ) const
{
    if ( numVals == 0 )
        return;
    if ( sizeof( size_t ) == sizeof( PetscInt ) ) {
        VecGetValues( d_petscVec, numVals, (PetscInt *) ndx, vals );
    } else {
        AMP_ASSERT( getGlobalSize() < 0x80000000 );
        std::vector<PetscInt> ndx2( numVals, 0 );
        for ( int i = 0; i < numVals; i++ )
            ndx2[i] = (PetscInt) ndx[i];
        VecGetValues( d_petscVec, numVals, &ndx2[0], vals );
    }
}

void NativePetscVectorData::assemble()
{
    VecAssemblyBegin( d_petscVec );
    VecAssemblyEnd( d_petscVec );
}


} // namespace LinearAlgebra
} // namespace AMP
