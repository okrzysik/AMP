#include "vectors/petsc/NativePetscVector.h"
#include "vectors/petsc/ManagedPetscVector.h"

#include "petsc.h"
#include "petscvec.h"


namespace AMP {
namespace LinearAlgebra {

NativePetscVector::NativePetscVector( VectorParameters::shared_ptr in_params )
    : NativeVector(), PetscVector(), VectorEngine()
{
    auto npvParams = dynamic_pointer_cast<NativePetscVectorParameters>( in_params );
    d_petscVec     = npvParams->d_InVec;
    d_pArray       = nullptr;
    CommunicationListParameters::shared_ptr params( new CommunicationListParameters() );
    params->d_comm      = npvParams->d_Comm;
    params->d_localsize = npvParams->d_localsize;
    setCommunicationList( CommunicationList::shared_ptr( new CommunicationList( params ) ) );
    d_bDeleteMe  = npvParams->d_Deleteable;
    d_DOFManager = AMP::Discretization::DOFManager::shared_ptr(
        new AMP::Discretization::DOFManager( npvParams->d_localsize, npvParams->d_Comm ) );
}


NativePetscVector::~NativePetscVector()
{
    resetArray();
    if ( d_bDeleteMe )
        PETSC::vecDestroy( &d_petscVec );
}


Vector::shared_ptr NativePetscVector::cloneVector( const Variable::shared_ptr var ) const
{
    resetArray();
    Vec new_petscVec;
    VecDuplicate( d_petscVec, &new_petscVec );
    auto npvParams            = make_shared<NativePetscVectorParameters>( new_petscVec, true );
    npvParams->d_Comm         = getComm();
    Vector::shared_ptr retVal = Vector::shared_ptr( new NativePetscVector( npvParams ) );
    retVal->setVariable( var );
    return retVal;
}

void NativePetscVector::putRawData( const double *in )
{
    int a, b;
    VecGetOwnershipRange( d_petscVec, &a, &b );
    AMP_ASSERT( b - a == (int) getLocalSize() );
    std::vector<int> offs( b - a );
    for ( size_t j = 0; j != offs.size(); j++ )
        offs[j]    = a + j;
    VecSetValues( d_petscVec, offs.size(), offs.data(), in, INSERT_VALUES );
}


void NativePetscVector::copyOutRawData( double *out ) const
{
    std::copy( getRawDataBlock<double>( 0 ), getRawDataBlock<double>( 0 ) + getLocalSize(), out );
}
}
}
