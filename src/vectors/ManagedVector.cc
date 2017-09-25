#include "vectors/ManagedVector.h"
#include "utils/Utilities.h"
#include <iostream>
#include <stdexcept>
#include <string>


namespace AMP {
namespace LinearAlgebra {


// Helper functions
static inline ManagedVector *getManaged( Vector *x )
{
    auto y = dynamic_cast<ManagedVector *>( x );
    AMP_INSIST( y != nullptr, "x is not a ManagedVector" );
    return y;
}
static inline AMP::shared_ptr<ManagedVector> getManaged( AMP::shared_ptr<Vector> x )
{
    auto y = dynamic_pointer_cast<ManagedVector>( x );
    AMP_INSIST( y != nullptr, "x is not a ManagedVector" );
    return y;
}


/********************************************************
* Constructors                                          *
********************************************************/
ManagedVector::ManagedVector( VectorParameters::shared_ptr params_in ) : Vector( params_in )
{
    d_pParameters = AMP::dynamic_pointer_cast<ManagedVectorParameters>( params_in );
    if ( d_pParameters->d_Buffer.get() != nullptr )
        d_vBuffer = d_pParameters->d_Buffer;
    else
        d_vBuffer = d_pParameters->d_Engine->getNewBuffer();
    if ( d_pParameters->d_CloneEngine )
        d_Engine = d_pParameters->d_Engine->cloneEngine( d_vBuffer );
    else
        d_Engine                 = d_pParameters->d_Engine;
    d_pParameters->d_CloneEngine = true;
}
ManagedVector::ManagedVector( shared_ptr alias )
    : Vector( AMP::dynamic_pointer_cast<VectorParameters>( getManaged( alias )->getParameters() ) )
{
    auto vec      = getManaged( alias );
    d_vBuffer     = vec->d_vBuffer;
    d_Engine      = vec->d_Engine;
    d_pParameters = vec->d_pParameters;
    setVariable( vec->getVariable() );
    aliasGhostBuffer( vec );
}


/********************************************************
* Subset                                                *
********************************************************/
Vector::shared_ptr ManagedVector::subsetVectorForVariable( Variable::const_shared_ptr name )
{
    Vector::shared_ptr retVal;
    if ( !d_vBuffer )
        retVal = dynamic_pointer_cast<Vector>( d_Engine )->subsetVectorForVariable( name );
    if ( !retVal )
        retVal = Vector::subsetVectorForVariable( name );
    return retVal;
}
Vector::const_shared_ptr
ManagedVector::constSubsetVectorForVariable( Variable::const_shared_ptr name ) const
{
    Vector::const_shared_ptr retVal;
    if ( !d_vBuffer )
        retVal = dynamic_pointer_cast<Vector>( d_Engine )->constSubsetVectorForVariable( name );
    if ( !retVal )
        retVal = Vector::constSubsetVectorForVariable( name );
    return retVal;
}


bool ManagedVector::isAnAliasOf( Vector &rhs )
{
    bool retVal = false;
    auto other  = getManaged( &rhs );
    if ( other != nullptr ) {
        if ( d_vBuffer && ( other->d_vBuffer == d_vBuffer ) ) {
            retVal = true;
        }
    }
    return retVal;
}


void ManagedVector::copy( const VectorOperations &other )
{
    auto rhs_managed = dynamic_cast<const ManagedVector *>( &other );
    AMP::shared_ptr<Vector> vec1;
    AMP::shared_ptr<const Vector> vec2;
    if ( rhs_managed != nullptr ) {
        // We are dealing with two managed vectors, check if they both have data engines
        if ( d_Engine.get() != nullptr )
            vec1 = AMP::dynamic_pointer_cast<Vector>( d_Engine );
        if ( rhs_managed->d_Engine.get() != nullptr )
            vec2 = AMP::dynamic_pointer_cast<const Vector>( rhs_managed->d_Engine );
    }
    // Perform the copy
    if ( vec1 != nullptr && vec2 != nullptr ) {
        // We have two data engines, perform the copy between them
        vec1->copy( vec2 );
        fireDataChange();
        *d_UpdateState = *( other.getVectorData()->getUpdateStatusPtr() );
    } else {
        // Default, general case
        VectorOperationsDefault::copy( other );
    }
}


Vector::UpdateState ManagedVector::getUpdateStatus() const
{
    Vector::UpdateState state = *d_UpdateState;
    AMP::shared_ptr<const Vector> vec;
    if ( d_Engine.get() != nullptr ) {
        vec = AMP::dynamic_pointer_cast<const Vector>( d_Engine );
    }
    if ( vec.get() != nullptr ) {
        Vector::UpdateState sub_state = vec->getUpdateStatus();
        if ( sub_state == UpdateState::UNCHANGED ) {
            // No change in state
        } else if ( sub_state == UpdateState::LOCAL_CHANGED && state == UpdateState::UNCHANGED ) {
            state = UpdateState::LOCAL_CHANGED;
        } else if ( sub_state == UpdateState::LOCAL_CHANGED ) {
            // No change in state
        } else if ( sub_state == UpdateState::ADDING &&
                    ( state == UpdateState::UNCHANGED || state == UpdateState::LOCAL_CHANGED ||
                      state == UpdateState::ADDING ) ) {
            state = UpdateState::ADDING;
        } else if ( sub_state == UpdateState::SETTING &&
                    ( state == UpdateState::UNCHANGED || state == UpdateState::LOCAL_CHANGED ||
                      state == UpdateState::SETTING ) ) {
            state = UpdateState::SETTING;
        } else {
            state = UpdateState::MIXED;
        }
    }
    return state;
}


void ManagedVector::setUpdateStatus( UpdateState state )
{
    *d_UpdateState = state;
    AMP::shared_ptr<Vector> vec;
    if ( d_Engine.get() != nullptr )
        vec = AMP::dynamic_pointer_cast<Vector>( d_Engine );
    if ( vec.get() != nullptr )
        vec->setUpdateStatus( state );
}


void ManagedVector::swapVectors( Vector &other )
{
    auto in = getManaged( &other );
    d_vBuffer.swap( in->d_vBuffer );
    std::swap( d_pParameters, in->d_pParameters );
    d_Engine->swapEngines( in->d_Engine );
}


void ManagedVector::aliasVector( Vector &other )
{
    auto in       = getManaged( &other );
    d_pParameters = in->d_pParameters;
    d_vBuffer     = in->d_vBuffer;
}

void ManagedVector::getValuesByGlobalID( int numVals, size_t *ndx, double *vals ) const
{
    Vector::shared_ptr vec = AMP::dynamic_pointer_cast<Vector>( d_Engine );
    if ( vec.get() == nullptr ) {
        Vector::getValuesByGlobalID( numVals, ndx, vals );
    } else {
        vec->getValuesByGlobalID( numVals, ndx, vals );
    }
}

void ManagedVector::getLocalValuesByGlobalID( int numVals, size_t *ndx, double *vals ) const
{
    if ( d_vBuffer ) {
        for ( int i = 0; i != numVals; i++ )
            vals[i] = ( *d_vBuffer )[ndx[i] - d_CommList->getStartGID()];
    } else {
        Vector::shared_ptr vec = AMP::dynamic_pointer_cast<Vector>( d_Engine );
        vec->getLocalValuesByGlobalID( numVals, ndx, vals );
    }
}

void ManagedVector::getGhostValuesByGlobalID( int numVals, size_t *ndx, double *vals ) const
{
    Vector::shared_ptr vec = AMP::dynamic_pointer_cast<Vector>( d_Engine );
    if ( vec.get() == nullptr ) {
        Vector::getGhostValuesByGlobalID( numVals, ndx, vals );
    } else {
        vec->getGhostValuesByGlobalID( numVals, ndx, vals );
    }
}

void ManagedVector::setValuesByLocalID( int i, size_t *id, const double *val )
{
    AMP_ASSERT( *d_UpdateState != UpdateState::ADDING );
    if ( *d_UpdateState == UpdateState::UNCHANGED )
        *d_UpdateState = UpdateState::LOCAL_CHANGED;
    d_Engine->setValuesByLocalID( i, id, val );
    fireDataChange();
}

void ManagedVector::setLocalValuesByGlobalID( int numVals, size_t *ndx, const double *vals )
{
    AMP_ASSERT( *d_UpdateState != UpdateState::ADDING );
    if ( *d_UpdateState == UpdateState::UNCHANGED )
        *d_UpdateState = UpdateState::LOCAL_CHANGED;
    d_Engine->setLocalValuesByGlobalID( numVals, ndx, vals );
    fireDataChange();
}

void ManagedVector::setGhostValuesByGlobalID( int numVals, size_t *ndx, const double *vals )
{
    Vector::shared_ptr vec = AMP::dynamic_pointer_cast<Vector>( d_Engine );
    if ( vec.get() == nullptr ) {
        Vector::setGhostValuesByGlobalID( numVals, ndx, vals );
    } else {
        vec->setGhostValuesByGlobalID( numVals, ndx, vals );
    }
}

void ManagedVector::setValuesByGlobalID( int numVals, size_t *ndx, const double *vals )
{
    Vector::shared_ptr vec = AMP::dynamic_pointer_cast<Vector>( d_Engine );
    if ( vec.get() != nullptr ) {
        AMP_ASSERT( *d_UpdateState != UpdateState::ADDING );
        *d_UpdateState         = UpdateState::SETTING;
        Vector::shared_ptr vec = AMP::dynamic_pointer_cast<Vector>( d_Engine );
        vec->setValuesByGlobalID( numVals, ndx, vals );
        fireDataChange();
    } else {
        std::vector<size_t> local_ndx;
        local_ndx.reserve( numVals );
        std::vector<double> local_val;
        local_val.reserve( numVals );
        std::vector<size_t> ghost_ndx;
        ghost_ndx.reserve( numVals );
        std::vector<double> ghost_val;
        ghost_val.reserve( numVals );
        for ( int i = 0; i < numVals; i++ ) {
            if ( ( ndx[i] < getLocalStartID() ) ||
                 ( ndx[i] >= ( getLocalStartID() + getLocalMaxID() ) ) ) {
                ghost_ndx.push_back( ndx[i] );
                ghost_val.push_back( vals[i] );
            } else {
                local_ndx.push_back( ndx[i] );
                local_val.push_back( vals[i] );
            }
        }
        if ( !ghost_ndx.empty() )
            setGhostValuesByGlobalID( ghost_ndx.size(), &ghost_ndx[0], &ghost_val[0] );
        if ( !local_ndx.empty() )
            setLocalValuesByGlobalID( local_ndx.size(), &local_ndx[0], &local_val[0] );
    }
}

void ManagedVector::addValuesByLocalID( int i, size_t *id, const double *val )
{
    AMP_ASSERT( *d_UpdateState != UpdateState::SETTING );
    if ( *d_UpdateState == UpdateState::UNCHANGED )
        *d_UpdateState = UpdateState::LOCAL_CHANGED;
    d_Engine->addValuesByLocalID( i, id, val );
    fireDataChange();
}

void ManagedVector::addLocalValuesByGlobalID( int i, size_t *id, const double *val )
{
    AMP_ASSERT( *d_UpdateState != UpdateState::SETTING );
    if ( *d_UpdateState == UpdateState::UNCHANGED )
        *d_UpdateState = UpdateState::LOCAL_CHANGED;
    d_Engine->addLocalValuesByGlobalID( i, id, val );
    fireDataChange();
}

void ManagedVector::putRawData( const double *in ) { d_Engine->putRawData( in ); }

void ManagedVector::copyOutRawData( double *in ) const { d_Engine->copyOutRawData( in ); }


void ManagedVector::scale( double alpha, const VectorOperations &x )
{
    auto x2 = dynamic_cast<const ManagedVector *>( &x );
    if ( x2 != nullptr ) {
        d_Engine->scale( alpha, x2->d_Engine );
    } else {
        VectorOperationsDefault::scale( alpha, x );
    }
    dataChanged();
}


void ManagedVector::scale( double alpha )
{
    d_Engine->scale( alpha );
    dataChanged();
}


void ManagedVector::add( const VectorOperations &x, const VectorOperations &y )
{
    auto x2 = dynamic_cast<const ManagedVector *>( &x );
    auto y2 = dynamic_cast<const ManagedVector *>( &y );
    if ( x2 != nullptr && y2 != nullptr ) {
        d_Engine->add( x2->d_Engine, y2->d_Engine );
    } else {
        VectorOperationsDefault::add( x, y );
    }
    dataChanged();
}


void ManagedVector::subtract( const VectorOperations &x, const VectorOperations &y )
{
    auto x2 = dynamic_cast<const ManagedVector *>( &x );
    auto y2 = dynamic_cast<const ManagedVector *>( &y );
    if ( x2 != nullptr && y2 != nullptr ) {
        d_Engine->subtract( x2->d_Engine, y2->d_Engine );
    } else {
        VectorOperationsDefault::subtract( x, y );
    }
    dataChanged();
}


void ManagedVector::multiply( const VectorOperations &x, const VectorOperations &y )
{
    auto x2 = dynamic_cast<const ManagedVector *>( &x );
    auto y2 = dynamic_cast<const ManagedVector *>( &y );
    if ( x2 != nullptr && y2 != nullptr ) {
        d_Engine->multiply( x2->d_Engine, y2->d_Engine );
    } else {
        VectorOperationsDefault::multiply( x, y );
    }
    dataChanged();
}


void ManagedVector::divide( const VectorOperations &x, const VectorOperations &y )
{
    auto x2 = dynamic_cast<const ManagedVector *>( &x );
    auto y2 = dynamic_cast<const ManagedVector *>( &y );
    if ( x2 != nullptr && y2 != nullptr ) {
        d_Engine->divide( x2->d_Engine, y2->d_Engine );
    } else {
        VectorOperationsDefault::divide( x, y );
    }
    dataChanged();
}


void ManagedVector::reciprocal( const VectorOperations &x )
{
    auto x2 = dynamic_cast<const ManagedVector *>( &x );
    if ( x2 != nullptr ) {
        d_Engine->reciprocal( x2->d_Engine );
    } else {
        VectorOperationsDefault::reciprocal( x );
    }
    dataChanged();
}


void ManagedVector::linearSum( double alpha,
                               const VectorOperations &x,
                               double beta,
                               const VectorOperations &y )
{
    auto x2 = dynamic_cast<const ManagedVector *>( &x );
    auto y2 = dynamic_cast<const ManagedVector *>( &y );
    if ( x2 != nullptr && y2 != nullptr ) {
        d_Engine->linearSum( alpha, x2->d_Engine, beta, y2->d_Engine );
    } else {
        VectorOperationsDefault::linearSum( alpha, x, beta, y );
    }
    dataChanged();
}


void ManagedVector::axpy( double alpha, const VectorOperations &x, const VectorOperations &y )
{
    auto x2 = dynamic_cast<const ManagedVector *>( &x );
    auto y2 = dynamic_cast<const ManagedVector *>( &y );
    if ( x2 != nullptr && y2 != nullptr ) {
        d_Engine->axpy( alpha, x2->d_Engine, y2->d_Engine );
    } else {
        VectorOperationsDefault::axpy( alpha, x, y );
    }
    dataChanged();
}


void ManagedVector::axpby( double alpha, double beta, const VectorOperations &x )
{
    auto x2 = dynamic_cast<const ManagedVector *>( &x );
    if ( x2 != nullptr ) {
        d_Engine->axpby( alpha, beta, x2->d_Engine );
    } else {
        VectorOperationsDefault::axpby( alpha, beta, x );
    }
    dataChanged();
}


void ManagedVector::abs( const VectorOperations &x )
{
    auto x2 = dynamic_cast<const ManagedVector *>( &x );
    if ( x2 != nullptr ) {
        d_Engine->abs( x2->d_Engine );
    } else {
        VectorOperationsDefault::abs( x );
    }
    dataChanged();
}


double ManagedVector::min( void ) const { return d_Engine->min(); }


double ManagedVector::max( void ) const { return d_Engine->max(); }


void ManagedVector::setRandomValues( void )
{
    d_Engine->setRandomValues();
    dataChanged();
    this->makeConsistent( ScatterType::CONSISTENT_SET );
}


void ManagedVector::setToScalar( double alpha )
{
    d_Engine->setToScalar( alpha );
    dataChanged();
    this->makeConsistent( ScatterType::CONSISTENT_SET );
}


double ManagedVector::L1Norm( void ) const { return d_Engine->L1Norm(); }


double ManagedVector::L2Norm( void ) const { return d_Engine->L2Norm(); }


double ManagedVector::maxNorm( void ) const { return d_Engine->maxNorm(); }


double ManagedVector::dot( const VectorOperations &x ) const
{
    auto x2 = dynamic_cast<const ManagedVector *>( &x );
    if ( x2 != nullptr )
        return d_Engine->dot( x2->d_Engine );
    return VectorOperationsDefault::dot( *this );
}


AMP::shared_ptr<Vector> ManagedVector::cloneVector( const Variable::shared_ptr name ) const
{
    AMP::shared_ptr<Vector> retVal( getNewRawPtr() );
    if ( !d_vBuffer ) {
        getManaged( retVal )->d_Engine = d_Engine->cloneEngine( VectorEngine::BufferPtr() );
    }
    retVal->setVariable( name );
    return retVal;
}
}
}
