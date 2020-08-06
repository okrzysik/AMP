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
    return engine->getVectorData();
}
static inline const VectorData *getEngineData( const VectorData &x )
{
    auto y = dynamic_cast<const ManagedVector *>( &x );    
    AMP_INSIST( y != nullptr, "x is not a ManagedVector" );
    auto engine = y->getVectorEngine();
    AMP_INSIST( engine, "ManagedVector Engine is Null" );
    return engine->getVectorData();
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
    d_Engine->getVectorData()->setValuesByLocalID( i, id, val );
    fireDataChange();
}

void ManagedVector::setLocalValuesByGlobalID( int numVals, size_t *ndx, const double *vals )
{
    AMP_ASSERT( *d_UpdateState != UpdateState::ADDING );
    if ( *d_UpdateState == UpdateState::UNCHANGED )
        *d_UpdateState = UpdateState::LOCAL_CHANGED;
    d_Engine->getVectorData()->setLocalValuesByGlobalID( numVals, ndx, vals );
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
    d_Engine->getVectorData()->addValuesByLocalID( i, id, val );
    fireDataChange();
}

void ManagedVector::addLocalValuesByGlobalID( int i, size_t *id, const double *val )
{
    AMP_ASSERT( *d_UpdateState != UpdateState::SETTING );
    if ( *d_UpdateState == UpdateState::UNCHANGED )
        *d_UpdateState = UpdateState::LOCAL_CHANGED;

    d_Engine->getVectorData()->addLocalValuesByGlobalID( i, id, val );
    fireDataChange();
}

void ManagedVector::putRawData( const double *in ) { d_Engine->getVectorData()->putRawData( in ); }

void ManagedVector::copyOutRawData( double *in ) const
{
    d_Engine->getVectorData()->copyOutRawData( in );
}

std::shared_ptr<Vector> ManagedVector::cloneVector( const Variable::shared_ptr name ) const
{
    std::shared_ptr<Vector> retVal( getNewRawPtr() );
    auto vec = std::dynamic_pointer_cast<Vector>( d_Engine );
    if ( vec ) {
        auto vec2                       = vec->cloneVector( "ManagedVectorClone" );
        getManaged( retVal )->d_vBuffer = std::dynamic_pointer_cast<VectorData>( vec2 );
        getManaged( retVal )->d_Engine  = std::dynamic_pointer_cast<VectorOperations>( vec2 );
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
    return d_Engine->getVectorData()->getRawDataBlockAsVoid( i );
}


const void *ManagedVector::getRawDataBlockAsVoid( size_t i ) const
{
    return d_Engine->getVectorData()->getRawDataBlockAsVoid( i );
}


void ManagedVector::addCommunicationListToParameters( CommunicationList::shared_ptr comm )
{
    d_pParameters->d_CommList = comm;
}


size_t ManagedVector::numberOfDataBlocks() const
{
    return d_Engine->getVectorData()->numberOfDataBlocks();
}


size_t ManagedVector::sizeOfDataBlock( size_t i ) const
{
    return d_Engine->getVectorData()->sizeOfDataBlock( i );
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


size_t ManagedVector::getLocalSize() const { return d_Engine->getVectorData()->getLocalSize(); }


size_t ManagedVector::getGlobalSize() const { return d_Engine->getVectorData()->getGlobalSize(); }


void ManagedVector::copy( const VectorOperations &x )
{
  copy(*(x.getVectorData()), *getVectorData());
}

void ManagedVector::setToScalar( double alpha )
{
  setToScalar(alpha, *getVectorData() );
}

void ManagedVector::setRandomValues( void )
{
   setRandomValues( *getVectorData() );
}

void ManagedVector::scale( double alpha, const VectorOperations &x )
{
  scale(alpha, *(x.getVectorData()), *getVectorData());
}


void ManagedVector::scale( double alpha )
{
  scale(alpha, *getVectorData());
}

void ManagedVector::add( const VectorOperations &x, const VectorOperations &y )
{
  add( *(x.getVectorData()), *(y.getVectorData()), *getVectorData() );
}


void ManagedVector::subtract( const VectorOperations &x, const VectorOperations &y )
{
  subtract( *(x.getVectorData()), *(y.getVectorData()), *getVectorData() );
}


void ManagedVector::multiply( const VectorOperations &x, const VectorOperations &y )
{
  multiply( *(x.getVectorData()), *(y.getVectorData()), *getVectorData() );
}


void ManagedVector::divide( const VectorOperations &x, const VectorOperations &y )
{
  divide( *(x.getVectorData()), *(y.getVectorData()), *getVectorData() );
}


void ManagedVector::reciprocal( const VectorOperations &x )
{
  reciprocal( *(x.getVectorData()), *getVectorData() );
}


void ManagedVector::linearSum( double alpha,
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


void ManagedVector::axpy( double alpha, const VectorOperations &x, const VectorOperations &y )
{
  axpy( alpha,
	*(x.getVectorData()),
	*(y.getVectorData()),
	*getVectorData() );
}


void ManagedVector::axpby( double alpha, double beta, const VectorOperations &x )
{
  axpby( alpha,
	 beta,
	 *(x.getVectorData()),
	 *getVectorData() );
}


void ManagedVector::abs( const VectorOperations &x )
{
    abs( *(x.getVectorData()), *getVectorData() );
}

double ManagedVector::min( void ) const
{
  return min( *getVectorData() );
}

double ManagedVector::max( void ) const
{
  return max( *getVectorData() );
}

double ManagedVector::L1Norm( void ) const
{
  return L1Norm( *getVectorData() );
}


double ManagedVector::L2Norm( void ) const
{
  return L2Norm( *getVectorData() );
}


double ManagedVector::maxNorm( void ) const
{
  return maxNorm( *getVectorData() );
}


double ManagedVector::dot( const VectorOperations &x ) const
{
    return dot( *(x.getVectorData()), *getVectorData() );
}

//**********************************************************************
// Functions that operate on VectorData objects

void ManagedVector::copy( const VectorData &src, VectorData &dst )
{
    auto src_managed  = getManagedVector(src);
    auto dst_managed  = getManagedVector(dst);
    std::shared_ptr<Vector> vec1;
    std::shared_ptr<const Vector> vec2;
    if ( src_managed && dst_managed ) {
        // We are dealing with two managed vectors, check if they both have data engines
        if ( dst_managed->d_Engine.get() )
            vec1 = std::dynamic_pointer_cast<Vector>( dst_managed->d_Engine );
        if ( src_managed->d_Engine.get() != nullptr )
            vec2 = std::dynamic_pointer_cast<const Vector>( src_managed->d_Engine );
    }
    // Perform the copy
    if ( vec1 != nullptr && vec2 != nullptr ) {
        // We have two data engines, perform the copy between them
        vec1->copy( vec2, vec1 );
        fireDataChange();
        *(dst_managed->d_UpdateState) = *( src.getUpdateStatusPtr() );
    } else {
        // Default, general case
      VectorOperationsDefault::copy( src, dst );
    }
}

void ManagedVector::setToScalar( double alpha, VectorData &x )
{
  auto xm    = getManagedVector(x);
  auto xdata = getEngineData(x);
  xm->d_Engine->setToScalar( alpha, *xdata );
  dataChanged();
  xm->makeConsistent( ScatterType::CONSISTENT_SET );

}

void ManagedVector::setRandomValues( VectorData &x )
{
    auto xm    = getManagedVector(x);
    auto xdata = getEngineData(x);
    xm->d_Engine->setRandomValues(*xdata);
    dataChanged();
    xm->makeConsistent( ScatterType::CONSISTENT_SET );
}

void ManagedVector::scale( double alpha, const VectorData &x, VectorData &y )
{
    auto x2 = getManagedVector(x);
    if ( x2 != nullptr ) {
        auto y2 = getManagedVector(y);
        y2->d_Engine->scale( alpha, *getEngineData(x), *getEngineData(y) );
    } else {
      VectorOperationsDefault::scale( alpha, x, y );
    }
    dataChanged();
}

void ManagedVector::scale( double alpha, VectorData &x )
{
    auto y = getManagedVector(x);
    y->d_Engine->scale( alpha, *getEngineData(x) );
    dataChanged();
}

void ManagedVector::add( const VectorData &x, const VectorData &y, VectorData &z )
{
    auto x2 = getManagedVector(x);
    auto y2 = getManagedVector(y);
    if ( x2 != nullptr && y2 != nullptr ) {
        auto z2 = getManagedVector(z);
        z2->d_Engine->add( *getEngineData(x), *getEngineData(y), *getEngineData(z) );
    } else {
      VectorOperationsDefault::add( x, y, z );
    }
    dataChanged();
}

void ManagedVector::subtract( const VectorData &x, const VectorData &y, VectorData &z  )
{
    auto x2 = getManagedVector(x);
    auto y2 = getManagedVector(y);
    if ( x2 != nullptr && y2 != nullptr ) {
        auto z2 = getManagedVector(z);
        z2->d_Engine->subtract( *getEngineData(x), *getEngineData(y), *getEngineData(z) );
    } else {
      VectorOperationsDefault::subtract( x, y, z );
    }
    dataChanged();
}

void ManagedVector::multiply( const VectorData &x, const VectorData &y, VectorData &z )
{
    auto x2 = getManagedVector(x);
    auto y2 = getManagedVector(y);
    if ( x2 != nullptr && y2 != nullptr ) {
        auto z2 = getManagedVector(z);
        z2->d_Engine->multiply( *getEngineData(x), *getEngineData(y), *getEngineData(z) );
    } else {
      VectorOperationsDefault::multiply( x, y, z );
    }
    dataChanged();
}

void ManagedVector::divide( const VectorData &x, const VectorData &y, VectorData &z )
{
    auto x2 = getManagedVector(x);
    auto y2 = getManagedVector(y);
    if ( x2 != nullptr && y2 != nullptr ) {
        auto z2 = getManagedVector(z);
        z2->d_Engine->divide( *getEngineData(x), *getEngineData(y), *getEngineData(z) );
    } else {
      VectorOperationsDefault::divide( x, y, z );
    }
    dataChanged();
}

void ManagedVector::reciprocal( const VectorData &x, VectorData &y )
{
    auto x2 = getManagedVector(x);
    if ( x2 != nullptr ) {
        auto y2 = getManagedVector(y);
        y2->d_Engine->reciprocal( *getEngineData(x), *getEngineData(y) );
    } else {
      VectorOperationsDefault::reciprocal( x, y );
    }
    dataChanged();
}

void ManagedVector::linearSum( double alpha,
			       const VectorData &x,
			       double beta,
			       const VectorData &y,
			       VectorData &z)
{
    auto x2 = getManagedVector(x);
    auto y2 = getManagedVector(y);
    if ( x2 != nullptr && y2 != nullptr ) {
        auto z2 = getManagedVector(z);
        z2->d_Engine->linearSum( alpha, *getEngineData(x), beta, *getEngineData(y), *getEngineData(z) );
    } else {
      VectorOperationsDefault::linearSum( alpha, x, beta, y, z );
    }
    dataChanged();
}

void ManagedVector::axpy( double alpha, const VectorData &x, const VectorData &y, VectorData &z )
{
  linearSum(alpha, x, 1.0, y, z);
}

void ManagedVector::axpby( double alpha, double beta, const VectorData &x, VectorData &z )
{
  linearSum(alpha, x, beta, z, z);
}

void ManagedVector::abs( const VectorData &x, VectorData &y )
{
    auto x2 = getManagedVector(x);
    if ( x2 != nullptr ) {
        auto y2 = getManagedVector(y);      
        y2->d_Engine->abs( *getEngineData(x), *getEngineData(y) );
    } else {
      VectorOperationsDefault::abs( x, y );
    }
    dataChanged();
}

double ManagedVector::min( const VectorData &x )  const
{
    auto x2 = getManagedVector(x);
    return x2->d_Engine->min(*getEngineData(x));
}

double ManagedVector::max( const VectorData &x )  const
{
    auto x2 = getManagedVector(x);
    return x2->d_Engine->max(*getEngineData(x));
}

double ManagedVector::dot( const VectorData &x, const VectorData &y ) const
{
    auto x2 = getManagedVector(x);
    if ( x2 != nullptr ) {
      auto y2 = getManagedVector(y);      
      return y2->d_Engine->dot( *getEngineData(x), *getEngineData(y) );
    }
    return VectorOperationsDefault::dot( x, y );
}

double ManagedVector::L1Norm( const VectorData &x )  const
{
    auto x2 = getManagedVector(x);
    return x2->d_Engine->L1Norm(*getEngineData(x));
}

double ManagedVector::L2Norm( const VectorData &x ) const 
{
    auto x2 = getManagedVector(x);
    return x2->d_Engine->L2Norm(*getEngineData(x));
}

double ManagedVector::maxNorm( const VectorData &x )  const
{
    auto x2 = getManagedVector(x);
    return x2->d_Engine->maxNorm(*getEngineData(x));
}

} // namespace LinearAlgebra
} // namespace AMP
