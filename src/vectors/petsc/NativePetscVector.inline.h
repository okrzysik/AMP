#include "vectors/petsc/ManagedPetscVector.h"
#include "vectors/petsc/NativePetscVector.h"

#include "petscvec.h"


namespace AMP {
namespace LinearAlgebra {


inline NativePetscVectorParameters::NativePetscVectorParameters( Vec v, bool deleteable )
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


inline size_t NativePetscVector::numberOfDataBlocks() const { return 1; }


inline size_t NativePetscVector::sizeOfDataBlock( size_t i ) const
{
    if ( i != 0 )
        return 0;
    return getLocalSize();
}


inline VectorEngine::BufferPtr NativePetscVector::getNewBuffer()
{
    return VectorEngine::BufferPtr();
}


inline bool NativePetscVector::sameEngine( VectorEngine &e ) const
{
    return dynamic_cast<const NativePetscVector *>( &e ) != nullptr;
}


inline VectorEngine::shared_ptr NativePetscVector::cloneEngine( BufferPtr ) const
{
    return AMP::dynamic_pointer_cast<VectorEngine>( Vector::cloneVector( "engine_clone" ) );
}


inline void NativePetscVector::swapEngines( VectorEngine::shared_ptr p )
{
    Vector::shared_ptr p2 = AMP::dynamic_pointer_cast<Vector>( p );
    Vector::swapVectors( p2 );
}


inline void *NativePetscVector::getDataBlock( size_t i )
{
    return static_cast<void *>( getRawDataBlock<double>( i ) );
}


inline const void *NativePetscVector::getDataBlock( size_t i ) const
{
    return static_cast<const void *>( getRawDataBlock<double>( i ) );
}


inline void NativePetscVector::getValuesByLocalID( int numVals, size_t *ndx, double *vals ) const
{
    Vector::getValuesByLocalID( numVals, ndx, vals );
}


inline Vector::shared_ptr NativePetscVector::getManagedVectorCopy( AMP_MPI comm )
{
    Vector::shared_ptr pRet = ManagedPetscVector::createFromPetscVec( d_petscVec, comm );
    ManagedPetscVector::copyFromPetscVec( *AMP::dynamic_pointer_cast<ManagedPetscVector>( pRet ),
                                          d_petscVec );
    return pRet;
}


inline AMP_MPI NativePetscVector::getComm() const { return Vector::getComm(); }


inline Vector::shared_ptr NativePetscVector::getManagedVectorDuplicate( AMP_MPI comm )
{
    Vector::shared_ptr pRet = ManagedPetscVector::createFromPetscVec( d_petscVec, comm );
    return pRet;
}


inline AMP::shared_ptr<ParameterBase> NativePetscVector::getParameters() { return d_pParameters; }


inline void NativePetscVector::aliasVector( Vector & ) { AMP_ERROR( "not implemented" ); }


inline void NativePetscVector::resetArray()
{
    if ( d_pArray ) {
        VecRestoreArray( d_petscVec, &d_pArray );
        d_pArray = 0;
    }
}


inline void NativePetscVector::resetArray() const
{
    if ( d_pArray ) {
        VecRestoreArray( d_petscVec, &d_pArray );
        d_pArray = 0;
    }
}


inline void NativePetscVector::swapVectors( Vector &other )
{
    resetArray();
    VecSwap( d_petscVec, dynamic_cast<NativePetscVector *>( &other )->getVec() );
}


inline void NativePetscVector::copy( const VectorOperations &src )
{
    resetArray();
    VecCopy( dynamic_cast<const NativePetscVector *>( &src )->getVec(), d_petscVec );
    copyGhostValues( *dynamic_cast<const VectorData *>( &src ) );
}


inline void NativePetscVector::setToScalar( double alpha )
{
    resetArray();
    VecSet( d_petscVec, alpha );
}


inline void NativePetscVector::scale( double alpha, const VectorOperations &x )
{
    resetArray();
    copyVector( dynamic_cast<const NativePetscVector *>( &x )->shared_from_this() );
    VecScale( d_petscVec, alpha );
}


inline void NativePetscVector::scale( double alpha )
{
    resetArray();
    VecScale( d_petscVec, alpha );
}


// Function to perform  this = alpha x + beta y + gamma this
inline void NativePetscVector::axpbypcz( double alpha,
                                         const VectorOperations &vx,
                                         double beta,
                                         const VectorOperations &vy,
                                         double gamma )
{
    resetArray();
    Vec x = dynamic_cast<const NativePetscVector *>( &vx )->getVec();
    Vec y = dynamic_cast<const NativePetscVector *>( &vy )->getVec();
    Vec z = d_petscVec;
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


inline void NativePetscVector::add( const VectorOperations &x, const VectorOperations &y )
{
    axpbypcz( 1.0, x, 1.0, y, 0.0 );
}


inline void NativePetscVector::subtract( const VectorOperations &x, const VectorOperations &y )
{
    axpbypcz( 1.0, x, -1.0, y, 0.0 );
}


inline void NativePetscVector::multiply( const VectorOperations &x, const VectorOperations &y )
{
    resetArray();
    VecPointwiseMult( d_petscVec,
                      dynamic_cast<const NativePetscVector *>( &x )->getVec(),
                      dynamic_cast<const NativePetscVector *>( &y )->getVec() );
}


inline void NativePetscVector::divide( const VectorOperations &x, const VectorOperations &y )
{
    resetArray();
    VecPointwiseDivide( d_petscVec,
                        dynamic_cast<const NativePetscVector *>( &x )->getVec(),
                        dynamic_cast<const NativePetscVector *>( &y )->getVec() );
}


inline void NativePetscVector::reciprocal( const VectorOperations &x )
{
    resetArray();
    copyVector( dynamic_cast<const Vector *>( &x )->shared_from_this() );
    VecReciprocal( d_petscVec );
}


inline void NativePetscVector::linearSum( double alpha,
                                          const VectorOperations &x,
                                          double beta,
                                          const VectorOperations &y )
{
    axpbypcz( alpha, x, beta, y, 0.0 );
}


inline void
NativePetscVector::axpy( double alpha, const VectorOperations &x, const VectorOperations &y )
{
    axpbypcz( alpha, x, 1.0, y, 0.0 );
}


inline void NativePetscVector::axpby( double alpha, double beta, const VectorOperations &x )
{
    axpbypcz( alpha, x, beta, *this, 0.0 );
}


inline void NativePetscVector::abs( const VectorOperations &x )
{
    resetArray();
    copyVector( dynamic_cast<const Vector *>( &x )->shared_from_this() );
    VecAbs( d_petscVec );
}

inline double NativePetscVector::min( void ) const
{
    resetArray();
    double val;
    VecMin( d_petscVec, PETSC_NULL, &val );
    return val;
}

inline double NativePetscVector::max( void ) const
{
    resetArray();
    double val;
    VecMax( d_petscVec, PETSC_NULL, &val );
    return val;
}


inline void NativePetscVector::setRandomValues( void )
{
    resetArray();
    VecSetRandom( d_petscVec, getPetscRandom( getComm() ) );
}


inline double NativePetscVector::L1Norm( void ) const
{
    resetArray();
    double ans;
    VecNorm( d_petscVec, NORM_1, &ans );
    return ans;
}


inline double NativePetscVector::L2Norm( void ) const
{
    resetArray();
    double ans;
    VecNorm( d_petscVec, NORM_2, &ans );
    return ans;
}


inline double NativePetscVector::maxNorm( void ) const
{
    resetArray();
    double ans;
    VecNorm( d_petscVec, NORM_INFINITY, &ans );
    return ans;
}


inline double NativePetscVector::dot( const VectorOperations &x ) const
{
    resetArray();
    double ans;
    VecDot( d_petscVec, dynamic_cast<const NativePetscVector *>( &x )->getVec(), &ans );
    return ans;
}


inline double NativePetscVector::localL1Norm( void ) const
{
    resetArray();
    double ans;
    PetscErrorCode ierr;
    ierr = ( *d_petscVec->ops->norm_local )( d_petscVec, NORM_1, &ans );
    CHKERRQ( ierr );
    return ans;
}


inline double NativePetscVector::localL2Norm( void ) const
{
    resetArray();
    double ans;
    PetscErrorCode ierr;
    ierr = ( *d_petscVec->ops->norm_local )( d_petscVec, NORM_2, &ans );
    CHKERRQ( ierr );
    return ans;
}


inline double NativePetscVector::localMaxNorm( void ) const
{
    resetArray();
    double ans;
    PetscErrorCode ierr;
    ierr = ( *d_petscVec->ops->norm_local )( d_petscVec, NORM_INFINITY, &ans );
    CHKERRQ( ierr );
    return ans;
}


inline void NativePetscVector::setValuesByLocalID( int num, size_t *indices, const double *vals )
{
    for ( int i = 0; i != num; i++ )
        getRawDataBlock<double>()[indices[i]] = vals[i];
}


inline void
NativePetscVector::setLocalValuesByGlobalID( int num, size_t *indices, const double *vals )
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


inline void NativePetscVector::addValuesByLocalID( int num, size_t *indices, const double *vals )
{
    for ( int i = 0; i != num; i++ )
        getRawDataBlock<double>()[indices[i]] += vals[i];
}


inline void
NativePetscVector::addLocalValuesByGlobalID( int num, size_t *indices, const double *vals )
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


inline void NativePetscVector::assemble()
{
    VecAssemblyBegin( d_petscVec );
    VecAssemblyEnd( d_petscVec );
}


inline size_t NativePetscVector::getLocalSize() const
{
    int lsize;
    VecGetLocalSize( d_petscVec, &lsize );
    return static_cast<size_t>( lsize );
}

inline size_t NativePetscVector::getGlobalSize() const
{
    int gsize;
    VecGetSize( d_petscVec, &gsize );
    return static_cast<size_t>( gsize );
}


inline void *NativePetscVector::getRawDataBlockAsVoid( size_t i )
{
    if ( i > 0 )
        return 0;
    if ( d_pArray == 0 ) {
        VecGetArray( d_petscVec, &d_pArray );
    }
    return d_pArray;
}

inline const void *NativePetscVector::getRawDataBlockAsVoid( size_t i ) const
{
    if ( i > 0 )
        return 0;
    if ( d_pArray == 0 ) {
        VecGetArray( d_petscVec, &d_pArray );
    }
    return d_pArray;
}


inline void
NativePetscVector::getLocalValuesByGlobalID( int numVals, size_t *ndx, double *vals ) const
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
} // namespace LinearAlgebra
} // namespace AMP
