#include "AMP/vectors/petsc/NativePetscVectorOperations.h"
#include "AMP/vectors/petsc/NativePetscVectorData.h"

#include "AMP/vectors/petsc/ManagedPetscVector.h"

namespace AMP {
namespace LinearAlgebra {

//**********************************************************************
// Static functions that operate on VectorData objects
// Helper function
Vec NativePetscVectorOperations::getPetscVec( VectorData &vx )
{
    auto nx = dynamic_cast<NativePetscVectorData *>( &vx );
    nx->resetArray();
    return nx->getVec();
}

Vec NativePetscVectorOperations::getPetscVec( const VectorData &vx )
{
    auto nx = dynamic_cast<const NativePetscVectorData *>( &vx );
    nx->resetArray();
    return nx->getVec();
}

Vec NativePetscVectorOperations::getConstPetscVec( const VectorData &vx )
{
    auto nx = dynamic_cast<const NativePetscVectorData *>( &vx );
    AMP_ASSERT( nx );
    return nx->getVec();
}

NativePetscVectorData *NativePetscVectorOperations::getNativeVec( VectorData &vx )
{
    return dynamic_cast<NativePetscVectorData *>( &vx );
}

const NativePetscVectorData *NativePetscVectorOperations::getNativeVec( const VectorData &vx )
{
    return dynamic_cast<const NativePetscVectorData *>( &vx );
}

// Function to perform  this = alpha x + beta y + gamma z
void NativePetscVectorOperations::axpbypcz( double alpha,
                                            const VectorData &vx,
                                            double beta,
                                            const VectorData &vy,
                                            double gamma,
                                            VectorData &vz )
{
    Vec x = getConstPetscVec( vx );
    Vec y = getConstPetscVec( vy );
    Vec z = getPetscVec( vz );

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
void NativePetscVectorOperations::copy( const VectorData &x, VectorData &y )
{
    auto nx = getNativeVec( x );
    auto ny = getNativeVec( y );
    ny->resetArray();
    VecCopy( nx->getVec(), ny->getVec() );
    y.copyGhostValues( x );
}

void NativePetscVectorOperations::zero( VectorData &x ) { VecZeroEntries( getPetscVec( x ) ); }

void NativePetscVectorOperations::setToScalar( double alpha, VectorData &x )
{
    auto vec = getPetscVec( x );
    VecSet( vec, alpha );
}

void NativePetscVectorOperations::setRandomValues( VectorData &x )
{
    auto nx = getNativeVec( x );
    nx->resetArray();
    VecSetRandom( nx->getVec(), nx->getPetscRandom( nx->getComm() ) );
}

void NativePetscVectorOperations::scale( double alpha, const VectorData &x, VectorData &y )
{
    VecCopy( getConstPetscVec( x ), getPetscVec( y ) );
    VecScale( getPetscVec( y ), alpha );
}


void NativePetscVectorOperations::scale( double alpha, VectorData &x )
{
    VecScale( getPetscVec( x ), alpha );
}

void NativePetscVectorOperations::add( const VectorData &x, const VectorData &y, VectorData &z )
{
    axpbypcz( 1.0, x, 1.0, y, 0.0, z );
}


void NativePetscVectorOperations::subtract( const VectorData &x,
                                            const VectorData &y,
                                            VectorData &z )
{
    axpbypcz( 1.0, x, -1.0, y, 0.0, z );
}


void NativePetscVectorOperations::multiply( const VectorData &x,
                                            const VectorData &y,
                                            VectorData &z )
{
    VecPointwiseMult( getPetscVec( z ), getConstPetscVec( x ), getConstPetscVec( y ) );
}


void NativePetscVectorOperations::divide( const VectorData &x, const VectorData &y, VectorData &z )
{
    VecPointwiseDivide( getPetscVec( z ), getConstPetscVec( x ), getConstPetscVec( y ) );
}


void NativePetscVectorOperations::reciprocal( const VectorData &x, VectorData &y )
{
    VecCopy( getConstPetscVec( x ), getPetscVec( y ) );
    VecReciprocal( getPetscVec( y ) );
}


void NativePetscVectorOperations::linearSum(
    double alpha, const VectorData &x, double beta, const VectorData &y, VectorData &z )
{
    axpbypcz( alpha, x, beta, y, 0.0, z );
}


void NativePetscVectorOperations::axpy( double alpha,
                                        const VectorData &x,
                                        const VectorData &y,
                                        VectorData &z )
{
    axpbypcz( alpha, x, 1.0, y, 0.0, z );
}


void NativePetscVectorOperations::axpby( double alpha,
                                         double beta,
                                         const VectorData &x,
                                         VectorData &vz )
{
    auto &z = *getNativeVec( vz );
    axpbypcz( alpha, x, beta, z, 0.0, z );
}


void NativePetscVectorOperations::abs( const VectorData &x, VectorData &y )
{
    VecCopy( getConstPetscVec( x ), getPetscVec( y ) );
    VecAbs( getPetscVec( y ) );
}

void NativePetscVectorOperations::addScalar( const VectorData &x, double alpha, VectorData &y )
{
    auto py = getPetscVec( y );
    VecCopy( getConstPetscVec( x ), py );
    VecShift( py, alpha );
}

double NativePetscVectorOperations::min( const VectorData &x ) const
{
    double val;
    VecMin( getConstPetscVec( x ), PETSC_NULL, &val );
    return val;
}

double NativePetscVectorOperations::max( const VectorData &x ) const
{
    double val;
    VecMax( getConstPetscVec( x ), PETSC_NULL, &val );
    return val;
}

double NativePetscVectorOperations::L1Norm( const VectorData &x ) const
{
    double ans;
    VecNorm( getConstPetscVec( x ), NORM_1, &ans );
    return ans;
}

double NativePetscVectorOperations::L2Norm( const VectorData &x ) const
{
    double ans;
    VecNorm( getConstPetscVec( x ), NORM_2, &ans );
    return ans;
}

double NativePetscVectorOperations::maxNorm( const VectorData &x ) const
{
    double ans;
    VecNorm( getConstPetscVec( x ), NORM_INFINITY, &ans );
    return ans;
}

double NativePetscVectorOperations::dot( const VectorData &x, const VectorData &y ) const
{
    double ans;
    VecDot( getConstPetscVec( x ), getConstPetscVec( y ), &ans );
    return ans;
}

double NativePetscVectorOperations::localL1Norm( const VectorData &vx ) const
{
    Vec x = getPetscVec( vx );

    double ans;
    PetscErrorCode ierr;
    ierr = ( *x->ops->norm_local )( x, NORM_1, &ans );
    CHKERRQ( ierr );
    return ans;
}

double NativePetscVectorOperations::localL2Norm( const VectorData &vx ) const
{
    Vec x = getPetscVec( vx );

    double ans;
    PetscErrorCode ierr;

    ierr = ( *x->ops->norm_local )( x, NORM_2, &ans );
    CHKERRQ( ierr );
    return ans;
}

double NativePetscVectorOperations::localMaxNorm( const VectorData &vx ) const
{
    Vec x = getPetscVec( vx );

    double ans;
    PetscErrorCode ierr;

    ierr = ( *x->ops->norm_local )( x, NORM_INFINITY, &ans );
    CHKERRQ( ierr );
    return ans;
}

double NativePetscVectorOperations::localDot( const VectorData &vx, const VectorData &vy ) const
{
    Vec x = getPetscVec( vx );

    PetscScalar ans;
    PetscErrorCode ierr;

    ierr = ( *x->ops->dot_local )( getConstPetscVec( vx ), getConstPetscVec( vy ), &ans );
    CHKERRQ( ierr );
    return static_cast<double>( ans );
}

} // namespace LinearAlgebra
} // namespace AMP
