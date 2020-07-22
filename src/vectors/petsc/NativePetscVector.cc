#include "AMP/vectors/petsc/NativePetscVector.h"

#include "AMP/vectors/petsc/ManagedPetscVector.h"

#include "petsc.h"
#include "petsc/private/vecimpl.h"
#include "petscvec.h"

namespace AMP {
namespace LinearAlgebra {

NativePetscVector::NativePetscVector( VectorParameters::shared_ptr in_params )
  : NativeVector(), PetscVector()
{
    auto npvParams = std::dynamic_pointer_cast<NativePetscVectorParameters>( in_params );
    d_petscVec     = npvParams->d_InVec;
    d_pArray       = nullptr;
    CommunicationListParameters::shared_ptr params( new CommunicationListParameters() );
    params->d_comm      = npvParams->d_Comm;
    params->d_localsize = npvParams->d_localsize;
    setCommunicationList( std::make_shared<CommunicationList>( params ) );
    d_bDeleteMe  = npvParams->d_Deleteable;
    d_DOFManager = std::make_shared<AMP::Discretization::DOFManager>( npvParams->d_localsize,
                                                                      npvParams->d_Comm );
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
    auto npvParams            = std::make_shared<NativePetscVectorParameters>( new_petscVec, true );
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
        offs[j] = a + j;
    VecSetValues( d_petscVec, offs.size(), offs.data(), in, INSERT_VALUES );
}


void NativePetscVector::copyOutRawData( double *out ) const
{
    std::copy( getRawDataBlock<double>( 0 ), getRawDataBlock<double>( 0 ) + getLocalSize(), out );
}


void NativePetscVector::swapData( VectorData & ) { AMP_ERROR( "Not finished" ); }


double NativePetscVector::localL1Norm( void ) const
{
    resetArray();
    double ans;
    PetscErrorCode ierr;
    ierr = ( *d_petscVec->ops->norm_local )( d_petscVec, NORM_1, &ans );
    CHKERRQ( ierr );
    return ans;
}


double NativePetscVector::localL2Norm( void ) const
{
    resetArray();
    double ans;
    PetscErrorCode ierr;

    ierr = ( *d_petscVec->ops->norm_local )( d_petscVec, NORM_2, &ans );
    CHKERRQ( ierr );
    return ans;
}


double NativePetscVector::localMaxNorm( void ) const
{
    resetArray();
    double ans;
    PetscErrorCode ierr;

    ierr = ( *d_petscVec->ops->norm_local )( d_petscVec, NORM_INFINITY, &ans );
    CHKERRQ( ierr );
    return ans;
}

} // namespace LinearAlgebra
} // namespace AMP
