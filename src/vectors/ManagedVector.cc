#include "utils/Utilities.h"
#include "ManagedVector.h"
#include <iostream>
#include <stdexcept>
#include <string>


namespace AMP {
namespace LinearAlgebra {


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
    : Vector( AMP::dynamic_pointer_cast<VectorParameters>(
          alias->castTo<ManagedVector>().getParameters() ) )
{
    d_vBuffer = alias->castTo<ManagedVector>().d_vBuffer;
    d_Engine  = alias->castTo<ManagedVector>().d_Engine;
    setVariable( alias->getVariable() );
    d_pParameters = alias->castTo<ManagedVector>().d_pParameters;
    aliasGhostBuffer( alias );
}


/********************************************************
* Subset                                                *
********************************************************/
Vector::shared_ptr ManagedVector::subsetVectorForVariable( const Variable::shared_ptr &name )
{
    Vector::shared_ptr retVal;
    if ( !d_vBuffer )
        retVal = d_Engine->castTo<Vector>().subsetVectorForVariable( name );
    if ( !retVal )
        retVal = Vector::subsetVectorForVariable( name );
    return retVal;
}
Vector::const_shared_ptr
ManagedVector::constSubsetVectorForVariable( const Variable::shared_ptr &name ) const
{
    Vector::const_shared_ptr retVal;
    if ( !d_vBuffer )
        retVal = d_Engine->castTo<const Vector>().constSubsetVectorForVariable( name );
    if ( !retVal )
        retVal = Vector::constSubsetVectorForVariable( name );
    return retVal;
}


bool ManagedVector::isAnAliasOf( Vector &rhs )
{
    bool retVal = false;
    if ( rhs.isA<ManagedVector>() ) {
        ManagedVector &other = rhs.castTo<ManagedVector>();
        if ( d_vBuffer && ( other.d_vBuffer == d_vBuffer ) ) {
            retVal = true;
        }
    }
    return retVal;
}


void ManagedVector::copyVector( Vector::const_shared_ptr other )
{
    AMP::shared_ptr<const ManagedVector> rhs_managed =
        AMP::dynamic_pointer_cast<const ManagedVector>( other );
    AMP::shared_ptr<Vector> vec1;
    AMP::shared_ptr<const Vector> vec2;
    if ( rhs_managed.get() != nullptr ) {
        // We are dealing with two managed vectors, check if they both have data engines
        if ( d_Engine.get() != nullptr )
            vec1 = AMP::dynamic_pointer_cast<Vector>( d_Engine );
        if ( rhs_managed->d_Engine.get() != nullptr )
            vec2 = AMP::dynamic_pointer_cast<const Vector>( rhs_managed->d_Engine );
    }
    if ( vec1.get() != nullptr && vec2.get() != nullptr ) {
        // We have two data engines, perform the copy between them
        vec1->copyVector( vec2 );
        fireDataChange();
        *d_UpdateState = *( other->getUpdateStatusPtr() );
        return;
    }
    // Default, general case
    if ( other->getLocalSize() != getLocalSize() ) { // Another error condition
        AMP_ERROR( "Destination vector and source vector not the same size" );
    }
    fireDataChange();
    VectorDataIterator cur1      = begin();
    VectorDataIterator end1      = end();
    ConstVectorDataIterator cur2 = other->begin();
    while ( cur1 != end1 ) {
        *cur1 = *cur2;
        ++cur1;
        ++cur2;
    }
    copyGhostValues( other );
    // Copy the consistency state from other
    *d_UpdateState = *( other->getUpdateStatusPtr() );
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
        if ( sub_state == UNCHANGED ) {
            // No change in state
        } else if ( sub_state == LOCAL_CHANGED && state == UNCHANGED ) {
            state = LOCAL_CHANGED;
        } else if ( sub_state == LOCAL_CHANGED ) {
            // No change in state
        } else if ( sub_state == ADDING &&
                    ( state == UNCHANGED || state == LOCAL_CHANGED || state == ADDING ) ) {
            state = ADDING;
        } else if ( sub_state == SETTING &&
                    ( state == UNCHANGED || state == LOCAL_CHANGED || state == SETTING ) ) {
            state = SETTING;
        } else {
            state = MIXED;
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
    ManagedVector &in = other.castTo<ManagedVector>();
    d_vBuffer.swap( in.d_vBuffer );
    std::swap( d_pParameters, in.d_pParameters );
    d_Engine->swapEngines( in.d_Engine );
}


void ManagedVector::aliasVector( Vector &other )
{
    ManagedVector &in = other.castTo<ManagedVector>();
    d_pParameters     = in.d_pParameters;
    d_vBuffer         = in.d_vBuffer;
}

void ManagedVector::getValuesByGlobalID( int numVals, size_t *ndx, double *vals ) const
{
    INCREMENT_COUNT( "Virtual" );
    Vector::shared_ptr vec = AMP::dynamic_pointer_cast<Vector>( d_Engine );
    if ( vec.get() == nullptr ) {
        Vector::getValuesByGlobalID( numVals, ndx, vals );
    } else {
        vec->getValuesByGlobalID( numVals, ndx, vals );
    }
}

void ManagedVector::getLocalValuesByGlobalID( int numVals, size_t *ndx, double *vals ) const
{
    INCREMENT_COUNT( "Virtual" );
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
    INCREMENT_COUNT( "Virtual" );
    Vector::shared_ptr vec = AMP::dynamic_pointer_cast<Vector>( d_Engine );
    if ( vec.get() == nullptr ) {
        Vector::getGhostValuesByGlobalID( numVals, ndx, vals );
    } else {
        vec->getGhostValuesByGlobalID( numVals, ndx, vals );
    }
}

void ManagedVector::setValuesByLocalID( int i, size_t *id, const double *val )
{
    INCREMENT_COUNT( "Virtual" );
    AMP_ASSERT( *d_UpdateState != ADDING );
    if ( *d_UpdateState == UNCHANGED )
        *d_UpdateState = LOCAL_CHANGED;
    d_Engine->setValuesByLocalID( i, id, val );
    fireDataChange();
}

void ManagedVector::setLocalValuesByGlobalID( int numVals, size_t *ndx, const double *vals )
{
    INCREMENT_COUNT( "Virtual" );
    AMP_ASSERT( *d_UpdateState != ADDING );
    if ( *d_UpdateState == UNCHANGED )
        *d_UpdateState = LOCAL_CHANGED;
    d_Engine->setLocalValuesByGlobalID( numVals, ndx, vals );
    fireDataChange();
}

void ManagedVector::setGhostValuesByGlobalID( int numVals, size_t *ndx, const double *vals )
{
    INCREMENT_COUNT( "Virtual" );
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
        INCREMENT_COUNT( "Virtual" );
        AMP_ASSERT( *d_UpdateState != ADDING );
        *d_UpdateState         = SETTING;
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
    INCREMENT_COUNT( "Virtual" );
    AMP_ASSERT( *d_UpdateState != SETTING );
    if ( *d_UpdateState == UNCHANGED )
        *d_UpdateState = LOCAL_CHANGED;
    d_Engine->addValuesByLocalID( i, id, val );
    fireDataChange();
}

void ManagedVector::addLocalValuesByGlobalID( int i, size_t *id, const double *val )
{
    INCREMENT_COUNT( "Virtual" );
    AMP_ASSERT( *d_UpdateState != SETTING );
    if ( *d_UpdateState == UNCHANGED )
        *d_UpdateState = LOCAL_CHANGED;
    d_Engine->addLocalValuesByGlobalID( i, id, val );
    fireDataChange();
}

void ManagedVector::putRawData( const double *in ) { d_Engine->putRawData( in ); }

void ManagedVector::copyOutRawData( double *in ) const { d_Engine->copyOutRawData( in ); }
}
}
