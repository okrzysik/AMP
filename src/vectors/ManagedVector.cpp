#include "AMP/vectors/ManagedVector.h"
#include "AMP/utils/Utilities.h"
#include "AMP/vectors/MultiVector.h"
#include <iostream>
#include <stdexcept>
#include <string>
#include <typeinfo>

namespace AMP {
namespace LinearAlgebra {


// Helper functions
static inline ManagedVector *getManaged( Vector *x )
{
    auto y = dynamic_cast<ManagedVector *>( x );
    AMP_INSIST( y != nullptr, "x is not a ManagedVector" );
    return y;
}
static inline std::shared_ptr<ManagedVector> getManaged( std::shared_ptr<Vector> x )
{
    auto y = std::dynamic_pointer_cast<ManagedVector>( x );
    AMP_INSIST( y != nullptr, "x is not a ManagedVector" );
    return y;
}
static inline const ManagedVector *getManagedVector( const VectorData &x )
{
    auto y = dynamic_cast<const ManagedVector *>( &x );
    return y;
}
static inline ManagedVector *getManagedVector( VectorData &x )
{
    auto y = dynamic_cast<ManagedVector *>( &x );
    AMP_INSIST( y != nullptr, "x is not a ManagedVector" );
    return y;
}
static inline VectorData *getEngineData( VectorData &x )
{
    auto y = dynamic_cast<ManagedVector *>( &x );
    AMP_INSIST( y != nullptr, "x is not a ManagedVector" );
    auto engine = y->getVectorEngine();
    AMP_INSIST( engine, "ManagedVector Engine is Null" );
    auto vecEngine = std::dynamic_pointer_cast<Vector>( engine );
    if ( vecEngine )
        return vecEngine->getVectorData();
    else {
        AMP_ERROR( "Not programmed for as yet" );
    }
    return nullptr;
}
static inline const VectorData *getEngineData( const VectorData &x )
{
    auto y = dynamic_cast<const ManagedVector *>( &x );
    AMP_INSIST( y != nullptr, "x is not a ManagedVector" );
    auto engine = y->getVectorEngine();
    AMP_INSIST( engine, "ManagedVector Engine is Null" );
    auto vecEngine = std::dynamic_pointer_cast<const Vector>( engine );
    if ( vecEngine )
        return vecEngine->getVectorData();
    else {
        AMP_ERROR( "Not programmed for as yet" );
    }
    return nullptr;
}


/********************************************************
 * Constructors                                          *
 ********************************************************/
ManagedVector::ManagedVector( VectorParameters::shared_ptr params_in )
    : Vector( params_in ),
      d_pParameters( std::dynamic_pointer_cast<ManagedVectorParameters>( params_in ) )
{
    d_vBuffer = d_pParameters->d_Buffer;
    d_Engine  = d_pParameters->d_Engine;
    AMP_ASSERT( d_vBuffer );
    AMP_ASSERT( d_Engine );
    d_vBuffer->setUpdateStatusPtr( getUpdateStatusPtr() );
    auto vec = std::dynamic_pointer_cast<Vector>( d_Engine );
    if ( vec )
        vec->setUpdateStatusPtr( getUpdateStatusPtr() );
}
ManagedVector::ManagedVector( shared_ptr alias )
    : Vector( std::dynamic_pointer_cast<VectorParameters>( getManaged( alias )->getParameters() ) )
{
    auto vec      = getManaged( alias );
    d_vBuffer     = vec->d_vBuffer;
    d_Engine      = vec->d_Engine;
    d_pParameters = vec->d_pParameters;
    setVariable( vec->getVariable() );
    aliasGhostBuffer( vec );
    if ( d_vBuffer )
        d_vBuffer->setUpdateStatusPtr( getUpdateStatusPtr() );
    auto vec2 = std::dynamic_pointer_cast<Vector>( d_Engine );
    if ( vec2 )
        vec2->setUpdateStatusPtr( getUpdateStatusPtr() );
}

void ManagedVector::initializeVectorOperations( void )
{
  std::cout << "Entering ManagedVector::initializeVectorOperations" << std::endl;
    d_VectorOps = new ManagedVectorOperations();
  std::cout << "Exiting ManagedVector::initializeVectorOperations" << std::endl;
}
  
ManagedVector::~ManagedVector()
{
}
  
/********************************************************
 * Subset                                                *
 ********************************************************/
Vector::shared_ptr ManagedVector::subsetVectorForVariable( Variable::const_shared_ptr name )
{
    Vector::shared_ptr retVal;
    if ( !retVal )
        retVal = Vector::subsetVectorForVariable( name );
    if ( !retVal ) {
        auto vec = std::dynamic_pointer_cast<Vector>( d_Engine );
        if ( vec )
            retVal = vec->subsetVectorForVariable( name );
    }
    return retVal;
}
Vector::const_shared_ptr
ManagedVector::constSubsetVectorForVariable( Variable::const_shared_ptr name ) const
{
    Vector::const_shared_ptr retVal;
    if ( !retVal )
        retVal = Vector::constSubsetVectorForVariable( name );
    if ( !retVal ) {
        auto vec = std::dynamic_pointer_cast<const Vector>( d_Engine );
        if ( vec )
            retVal = vec->constSubsetVectorForVariable( name );
    }
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

Vector::UpdateState ManagedVector::getUpdateStatus() const
{
    Vector::UpdateState state = *d_UpdateState;
    std::shared_ptr<const Vector> vec;
    if ( d_Engine.get() != nullptr ) {
        vec = std::dynamic_pointer_cast<const Vector>( d_Engine );
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
    std::shared_ptr<Vector> vec;
    if ( d_Engine.get() != nullptr )
        vec = std::dynamic_pointer_cast<Vector>( d_Engine );
    if ( vec.get() != nullptr )
        vec->setUpdateStatus( state );
}


void ManagedVector::swapVectors( Vector &other )
{
    auto in = getManaged( &other );
    std::swap( d_vBuffer, in->d_vBuffer );
    std::swap( d_Engine, in->d_Engine );
    std::swap( d_pParameters, in->d_pParameters );

    d_vBuffer->setUpdateStatusPtr( getUpdateStatusPtr() );
    auto vec = std::dynamic_pointer_cast<Vector>( d_Engine );
    if ( vec )
        vec->setUpdateStatusPtr( getUpdateStatusPtr() );

    in->d_vBuffer->setUpdateStatusPtr( in->getUpdateStatusPtr() );
    vec = std::dynamic_pointer_cast<Vector>( in->d_Engine );
    if ( vec )
        vec->setUpdateStatusPtr( in->getUpdateStatusPtr() );
}


void ManagedVector::aliasVector( Vector &other )
{
    auto in       = getManaged( &other );
    d_pParameters = in->d_pParameters;
    d_vBuffer     = in->d_vBuffer;
}

void ManagedVector::getValuesByGlobalID( int numVals, size_t *ndx, double *vals ) const
{
    Vector::shared_ptr vec = std::dynamic_pointer_cast<Vector>( d_Engine );
    if ( vec.get() == nullptr ) {
        Vector::getValuesByGlobalID( numVals, ndx, vals );
    } else {
        vec->getValuesByGlobalID( numVals, ndx, vals );
    }
}

void ManagedVector::getLocalValuesByGlobalID( int numVals, size_t *ndx, double *vals ) const
{
    d_vBuffer->getLocalValuesByGlobalID( numVals, ndx, vals );
}

void ManagedVector::getGhostValuesByGlobalID( int numVals, size_t *ndx, double *vals ) const
{
    Vector::shared_ptr vec = std::dynamic_pointer_cast<Vector>( d_Engine );
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
    getEngineData( *this )->setValuesByLocalID( i, id, val );
    fireDataChange();
}

void ManagedVector::setLocalValuesByGlobalID( int numVals, size_t *ndx, const double *vals )
{
    AMP_ASSERT( *d_UpdateState != UpdateState::ADDING );
    if ( *d_UpdateState == UpdateState::UNCHANGED )
        *d_UpdateState = UpdateState::LOCAL_CHANGED;
    getEngineData( *this )->setLocalValuesByGlobalID( numVals, ndx, vals );
    fireDataChange();
}

void ManagedVector::setGhostValuesByGlobalID( int numVals, size_t *ndx, const double *vals )
{
    Vector::shared_ptr vec = std::dynamic_pointer_cast<Vector>( d_Engine );
    if ( vec.get() == nullptr ) {
        Vector::setGhostValuesByGlobalID( numVals, ndx, vals );
    } else {
        vec->setGhostValuesByGlobalID( numVals, ndx, vals );
    }
}

void ManagedVector::setValuesByGlobalID( int numVals, size_t *ndx, const double *vals )
{
    Vector::shared_ptr vec = std::dynamic_pointer_cast<Vector>( d_Engine );
    if ( vec.get() != nullptr ) {
        AMP_ASSERT( *d_UpdateState != UpdateState::ADDING );
        *d_UpdateState = UpdateState::SETTING;
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
    getEngineData( *this )->addValuesByLocalID( i, id, val );
    fireDataChange();
}

void ManagedVector::addLocalValuesByGlobalID( int i, size_t *id, const double *val )
{
    AMP_ASSERT( *d_UpdateState != UpdateState::SETTING );
    if ( *d_UpdateState == UpdateState::UNCHANGED )
        *d_UpdateState = UpdateState::LOCAL_CHANGED;

    getEngineData( *this )->addLocalValuesByGlobalID( i, id, val );
    fireDataChange();
}

void ManagedVector::putRawData( const double *in ) { getEngineData( *this )->putRawData( in ); }

void ManagedVector::copyOutRawData( double *in ) const
{
    getEngineData( *this )->copyOutRawData( in );
}

std::shared_ptr<Vector> ManagedVector::cloneVector( const Variable::shared_ptr name ) const
{
    std::shared_ptr<Vector> retVal( getNewRawPtr() );
    auto vec = std::dynamic_pointer_cast<Vector>( d_Engine );
    if ( vec ) {
        auto vec2                       = vec->cloneVector( "ManagedVectorClone" );
        getManaged( retVal )->d_vBuffer = std::dynamic_pointer_cast<VectorData>( vec2 );
        getManaged( retVal )->d_Engine  = std::dynamic_pointer_cast<Vector>( vec2 );
    } else {
        AMP_ERROR( "ManagedVector::cloneVector() should not have reached here!" );
    }
    retVal->setVariable( name );
    return retVal;
}

std::string ManagedVector::type() const
{
    if ( d_vBuffer )
        return " ( managed data )";
    std::string retVal = " ( managed view of ";
    auto vec           = std::dynamic_pointer_cast<Vector>( d_Engine );
    retVal += vec->type();
    retVal += " )";
    return retVal;
}


Vector::shared_ptr ManagedVector::getRootVector()
{
    if ( d_vBuffer )
        return shared_from_this();
    auto vec = std::dynamic_pointer_cast<ManagedVector>( d_Engine );
    if ( vec != nullptr )
        return vec->getRootVector();
    return std::dynamic_pointer_cast<Vector>( d_Engine )->shared_from_this();
}


void ManagedVector::dataChanged()
{
    if ( *d_UpdateState == UpdateState::UNCHANGED )
        *d_UpdateState = UpdateState::LOCAL_CHANGED;
}


Vector::shared_ptr ManagedVector::selectInto( const VectorSelector &s )
{
    Vector::shared_ptr result;
    if ( d_vBuffer ) {
        result = Vector::selectInto( s );
    } else {
        result = std::dynamic_pointer_cast<Vector>( d_Engine )->selectInto( s );
    }
    return result;
}


Vector::const_shared_ptr ManagedVector::selectInto( const VectorSelector &s ) const
{
    Vector::const_shared_ptr result;
    if ( d_vBuffer ) {
        result = Vector::selectInto( s );
    } else {
        result = std::dynamic_pointer_cast<Vector>( d_Engine )->selectInto( s );
    }
    return result;
}


ManagedVectorParameters::ManagedVectorParameters() : d_Buffer( nullptr ) {}


void *ManagedVector::getRawDataBlockAsVoid( size_t i )
{
    return getEngineData( *this )->getRawDataBlockAsVoid( i );
}


const void *ManagedVector::getRawDataBlockAsVoid( size_t i ) const
{
    return getEngineData( *this )->getRawDataBlockAsVoid( i );
}


void ManagedVector::addCommunicationListToParameters( CommunicationList::shared_ptr comm )
{
    d_pParameters->d_CommList = comm;
}


size_t ManagedVector::numberOfDataBlocks() const
{
    return getEngineData( *this )->numberOfDataBlocks();
}


size_t ManagedVector::sizeOfDataBlock( size_t i ) const
{
    return getEngineData( *this )->sizeOfDataBlock( i );
}


bool ManagedVector::isAnAliasOf( Vector::shared_ptr rhs ) { return isAnAliasOf( *rhs ); }


std::shared_ptr<ParameterBase> ManagedVector::getParameters()
{
    return std::dynamic_pointer_cast<ParameterBase>( d_pParameters );
}


std::shared_ptr<ManagedVectorParameters> ManagedVector::getManagedVectorParameters()
{
    return d_pParameters;
}


size_t ManagedVector::getLocalSize() const { return getEngineData( *this )->getLocalSize(); }


size_t ManagedVector::getGlobalSize() const { return getEngineData( *this )->getGlobalSize(); }

} // namespace LinearAlgebra
} // namespace AMP
