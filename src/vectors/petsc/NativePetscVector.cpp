#include "AMP/vectors/petsc/NativePetscVector.h"

#include "AMP/vectors/petsc/ManagedPetscVector.h"

#include "petsc.h"
#include "petsc/private/vecimpl.h"
#include "petscvec.h"

namespace AMP {
namespace LinearAlgebra {

NativePetscVector::NativePetscVector( VectorParameters::shared_ptr in_params )
    : Vector(), PetscVector()
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

NativePetscVectorParameters::NativePetscVectorParameters( Vec v, bool deleteable )
{
    // Get the communicator from the PETSc vector
    d_InVec       = v;
    MPI_Comm comm = d_Comm.getCommunicator(); // Get a MPI_comm object from AMP_MPI to pass to PETSc
    PetscObjectGetComm( reinterpret_cast<PetscObject>( v ), &comm );
    if ( comm != d_Comm.getCommunicator() )
        d_Comm = AMP_MPI( comm );
    d_Deleteable = deleteable;
    int lsize;
    VecGetLocalSize( v, &lsize );
    d_localsize = (size_t) lsize;
}


size_t NativePetscVector::numberOfDataBlocks() const { return 1; }


size_t NativePetscVector::sizeOfDataBlock( size_t i ) const
{
    if ( i != 0 )
        return 0;
    return getLocalSize();
}

void NativePetscVector::getValuesByLocalID( int numVals, size_t *ndx, double *vals ) const
{
    Vector::getValuesByLocalID( numVals, ndx, vals );
}

std::shared_ptr<ParameterBase> NativePetscVector::getParameters() { return d_pParameters; }


void NativePetscVector::aliasVector( Vector & ) { AMP_ERROR( "not implemented" ); }


void NativePetscVector::resetArray()
{
    if ( d_pArray ) {
        VecRestoreArray( d_petscVec, &d_pArray );
        d_pArray = 0;
    }
}


void NativePetscVector::resetArray() const
{
    if ( d_pArray ) {
        VecRestoreArray( d_petscVec, &d_pArray );
        d_pArray = 0;
    }
}


void NativePetscVector::swapVectors( Vector &other )
{
    resetArray();
    VecSwap( d_petscVec, dynamic_cast<NativePetscVector *>( &other )->getVec() );
}

void NativePetscVector::setValuesByLocalID( int num, size_t *indices, const double *vals )
{
    for ( int i = 0; i != num; i++ )
        getRawDataBlock<double>()[indices[i]] = vals[i];
}


void NativePetscVector::setLocalValuesByGlobalID( int num, size_t *indices, const double *vals )
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


void NativePetscVector::addValuesByLocalID( int num, size_t *indices, const double *vals )
{
    for ( int i = 0; i != num; i++ )
        getRawDataBlock<double>()[indices[i]] += vals[i];
}


void NativePetscVector::addLocalValuesByGlobalID( int num, size_t *indices, const double *vals )
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


void NativePetscVector::assemble()
{
    VecAssemblyBegin( d_petscVec );
    VecAssemblyEnd( d_petscVec );
}


size_t NativePetscVector::getLocalSize() const
{
    int lsize;
    VecGetLocalSize( d_petscVec, &lsize );
    return static_cast<size_t>( lsize );
}

size_t NativePetscVector::getGlobalSize() const
{
    int gsize;
    VecGetSize( d_petscVec, &gsize );
    return static_cast<size_t>( gsize );
}


void *NativePetscVector::getRawDataBlockAsVoid( size_t i )
{
    if ( i > 0 )
        return 0;
    if ( d_pArray == 0 ) {
        VecGetArray( d_petscVec, &d_pArray );
    }
    return d_pArray;
}

const void *NativePetscVector::getRawDataBlockAsVoid( size_t i ) const
{
    if ( i > 0 )
        return 0;
    if ( d_pArray == 0 ) {
        VecGetArray( d_petscVec, &d_pArray );
    }
    return d_pArray;
}


void NativePetscVector::getLocalValuesByGlobalID( int numVals, size_t *ndx, double *vals ) const
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
#if 0

double NativePetscVector::localL1Norm( void ) const
{
  return localL1Norm( *getVectorData());
}


double NativePetscVector::localL2Norm( void ) const
{
  return localL2Norm( *getVectorData());
}


double NativePetscVector::localMaxNorm( void ) const
{
  return localMaxNorm( *getVectorData());
}

void NativePetscVector::copy( const VectorOperations &src )
{
    resetArray();
    VecCopy( dynamic_cast<const NativePetscVector *>( &src )->getVec(), d_petscVec );
    copyGhostValues( *dynamic_cast<const VectorData *>( &src ) );
}

void NativePetscVector::setToScalar( double alpha )
{
  setToScalar(alpha, *getVectorData() );
}

void NativePetscVector::setRandomValues( void )
{
   setRandomValues( *getVectorData() );
}

void NativePetscVector::scale( double alpha, const VectorOperations &x )
{
  scale(alpha, *(x.getVectorData()), *getVectorData());
}


void NativePetscVector::scale( double alpha )
{
  scale(alpha, *getVectorData());
}


// Function to perform  this = alpha x + beta y + gamma this
void NativePetscVector::axpbypcz( double alpha,
                                  const VectorOperations &vx,
                                  double beta,
                                  const VectorOperations &vy,
                                  double gamma )
{
  axpbypcz( alpha,
	    *(vx.getVectorData()),
	    beta,
	    *(vy.getVectorData()),
	    gamma,
	    *getVectorData() );
}


void NativePetscVector::add( const VectorOperations &x, const VectorOperations &y )
{
  add( *(x.getVectorData()), *(y.getVectorData()), *getVectorData() );
}


void NativePetscVector::subtract( const VectorOperations &x, const VectorOperations &y )
{
  subtract( *(x.getVectorData()), *(y.getVectorData()), *getVectorData() );
}


void NativePetscVector::multiply( const VectorOperations &x, const VectorOperations &y )
{
  multiply( *(x.getVectorData()), *(y.getVectorData()), *getVectorData() );
}


void NativePetscVector::divide( const VectorOperations &x, const VectorOperations &y )
{
  divide( *(x.getVectorData()), *(y.getVectorData()), *getVectorData() );
}


void NativePetscVector::reciprocal( const VectorOperations &x )
{
  reciprocal( *(x.getVectorData()), *getVectorData() );
}


void NativePetscVector::linearSum( double alpha,
                                   const VectorOperations &x,
                                   double beta,
                                   const VectorOperations &y )
{
  linearSum( alpha,
	     *(x.getVectorData()),
	     beta,
	     *(y.getVectorData()),
	     *getVectorData() );
}


void NativePetscVector::axpy( double alpha, const VectorOperations &x, const VectorOperations &y )
{
  axpy( alpha,
	*(x.getVectorData()),
	*(y.getVectorData()),
	*getVectorData() );
}


void NativePetscVector::axpby( double alpha, double beta, const VectorOperations &x )
{
  axpby( alpha,
	 beta,
	 *(x.getVectorData()),
	 *getVectorData() );
}


void NativePetscVector::abs( const VectorOperations &x )
{
    abs( *(x.getVectorData()), *getVectorData() );
}

double NativePetscVector::min( void ) const
{
  return min( *getVectorData() );
}

double NativePetscVector::max( void ) const
{
  return max( *getVectorData() );
}

double NativePetscVector::L1Norm( void ) const
{
  return L1Norm( *getVectorData() );
}


double NativePetscVector::L2Norm( void ) const
{
  return L2Norm( *getVectorData() );
}


double NativePetscVector::maxNorm( void ) const
{
  return maxNorm( *getVectorData() );
}


double NativePetscVector::dot( const VectorOperations &x ) const
{
    return dot( *(x.getVectorData()), *getVectorData() );
}
#endif

//**********************************************************************
// Static functions that operate on VectorData objects

// Helper function
Vec NativePetscVector::getPetscVec( VectorData &vx )
{
    auto nx = dynamic_cast<NativePetscVector *>( &vx );
    nx->resetArray();
    return nx->getVec();
}

Vec NativePetscVector::getPetscVec( const VectorData &vx )
{
    auto nx = dynamic_cast<const NativePetscVector *>( &vx );
    nx->resetArray();
    return nx->getVec();
}

Vec NativePetscVector::getConstPetscVec( const VectorData &vx )
{
    return dynamic_cast<const NativePetscVector *>( &vx )->getVec();
}

NativePetscVector *
NativePetscVector::getNativeVec( VectorData &vx  )
{
  return dynamic_cast<NativePetscVector *>( &vx );
}

const NativePetscVector *
NativePetscVector::getNativeVec( const VectorData &vx  )
{
  return dynamic_cast<const NativePetscVector *>( &vx );
}

// Function to perform  this = alpha x + beta y + gamma z
void NativePetscVector::axpbypcz( double alpha,
                                  const VectorData &vx,
                                  double beta,
                                  const VectorData &vy,
                                  double gamma,
				  VectorData &vz)
{
    Vec x = getConstPetscVec(vx);
    Vec y = getConstPetscVec(vy);
    Vec z = getPetscVec(vz);

    if ( x != y && x != z && y != z ) {
        // We can safely perform  z = alpha x + beta y + gamma z
        VecAXPBYPCZ( z, alpha, beta, gamma, x, y );
    } else if ( x != y && x == z ) {
        // x==z:  z = (alpha+gamma)*z + beta*y
        double scale = alpha + gamma;
        VecAXPBY( z, beta, scale, y );
    } else if ( x != y && y == z ) {
        // y==z:  z = (beta+gamma)*z + alpha*x
        double scale = beta + gamma;
        VecAXPBY( z, alpha, scale, x );
    } else if ( x == y && x == z ) {
        // x==y==z:  z = (alpha+beta+gamma)*z
        double scale = alpha + beta + gamma;
        VecScale( z, scale );
    } else {
        AMP_ERROR( "Internal error\n" );
    }
}
void NativePetscVector::copy( const VectorData &x, VectorData &y )
{
    auto nx = getNativeVec(x);
    auto ny = getNativeVec(y);
    ny->resetArray();    
    VecCopy( nx->getVec(), ny->getVec() );
    y.copyGhostValues(x);
}

void NativePetscVector::setToScalar( double alpha, VectorData &x )
{
    auto vec = getPetscVec(x);
    VecSet( vec, alpha );
}

void NativePetscVector::setRandomValues( VectorData &x )
{
    auto nx = getNativeVec(x);
    nx->resetArray();
    VecSetRandom( nx->getVec(), nx->getPetscRandom( nx->getComm() ) );
}

void NativePetscVector::scale( double alpha, const VectorData &x, VectorData &y )
{
    auto nx = getNativeVec(x);
    auto ny = getNativeVec(y);
    ny->resetArray();
    ny->copyVector( nx->shared_from_this() );
    VecScale( ny->getVec(), alpha );
}


void NativePetscVector::scale( double alpha, VectorData &x )
{
    VecScale( getPetscVec(x), alpha );
}

void NativePetscVector::add( const VectorData &x, const VectorData &y, VectorData &z )
{
  axpbypcz( 1.0, x, 1.0, y, 0.0, z );
}


void NativePetscVector::subtract( const VectorData &x, const VectorData &y, VectorData &z  )
{
  axpbypcz( 1.0, x, -1.0, y, 0.0, z );
}


void NativePetscVector::multiply( const VectorData &x, const VectorData &y, VectorData &z )
{
    VecPointwiseMult( getPetscVec(z),
                      getConstPetscVec(x),
                      getConstPetscVec(y) );
}


void NativePetscVector::divide( const VectorData &x, const VectorData &y, VectorData &z )
{
    VecPointwiseDivide( getPetscVec(z),
			getConstPetscVec(x),
			getConstPetscVec(y) );
}


void NativePetscVector::reciprocal( const VectorData &x, VectorData &y )
{
    auto nx = getNativeVec(x);
    auto ny = getNativeVec(y);
    ny->resetArray();
    ny->copyVector( nx->shared_from_this() );
    VecReciprocal( ny->getVec() );
}


void NativePetscVector::linearSum( double alpha,
                                   const VectorData &x,
                                   double beta,
                                   const VectorData &y,
				   VectorData &z)
{
  axpbypcz( alpha, x, beta, y, 0.0, z );
}


void NativePetscVector::axpy( double alpha, const VectorData &x, const VectorData &y, VectorData &z )
{
  axpbypcz( alpha, x, 1.0, y, 0.0, z );
}


void NativePetscVector::axpby( double alpha, double beta, const VectorData &x, VectorData &vz )
{
  auto &z = *getNativeVec(vz);
  axpbypcz( alpha, x, beta, z, 0.0, z );
}


void NativePetscVector::abs( const VectorData &x, VectorData &y )
{
    auto nx = getNativeVec(x);
    auto ny = getNativeVec(y);
    ny->resetArray();
    ny->copyVector( nx->shared_from_this() );
    VecAbs( ny->getVec() );
}

double NativePetscVector::min( const VectorData &x )  const
{
    double val;
    VecMin( getConstPetscVec(x), PETSC_NULL, &val );
    return val;
}

double NativePetscVector::max( const VectorData &x )  const
{
    double val;
    VecMax( getConstPetscVec(x), PETSC_NULL, &val );
    return val;
}

double NativePetscVector::dot( const VectorData &x, const VectorData &y ) const
{
    double ans;
    VecDot( getPetscVec(y), getConstPetscVec(x), &ans );
    return ans;
}

double NativePetscVector::L1Norm( const VectorData &x )  const
{
    double ans;
    VecNorm( getConstPetscVec(x), NORM_1, &ans );
    return ans;
}

double NativePetscVector::L2Norm( const VectorData &x ) const 
{
    double ans;
    VecNorm( getConstPetscVec(x), NORM_2, &ans );
    return ans;
}

double NativePetscVector::maxNorm( const VectorData &x )  const
{
    double ans;
    VecNorm( getConstPetscVec(x), NORM_INFINITY, &ans );
    return ans;
}

double NativePetscVector::localL1Norm( const VectorData &vx ) const
{
    Vec x = getPetscVec(vx);

    double ans;
    PetscErrorCode ierr;
    ierr = ( *x->ops->norm_local )( x, NORM_1, &ans );
    CHKERRQ( ierr );
    return ans;
}

double NativePetscVector::localL2Norm( const VectorData &vx  ) const
{
    Vec x = getPetscVec(vx);

    double ans;
    PetscErrorCode ierr;

    ierr = ( *x->ops->norm_local )( x, NORM_2, &ans );
    CHKERRQ( ierr );
    return ans;
}

double NativePetscVector::localMaxNorm( const VectorData &vx ) const
{
    Vec x = getPetscVec(vx);

    double ans;
    PetscErrorCode ierr;

    ierr = ( *x->ops->norm_local )( x, NORM_INFINITY, &ans );
    CHKERRQ( ierr );
    return ans;
}

} // namespace LinearAlgebra
} // namespace AMP
