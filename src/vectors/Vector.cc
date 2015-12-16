#include "Vector.h"


#include <float.h>
#include <math.h>
#include <string.h>
#include <typeinfo>

#include "utils/AMP_MPI.h"
#include "utils/Utilities.h"

#include "DataChangeFirer.h"

#include "MultiVector.h"
#include "VectorSelector.h"

namespace AMP {
namespace LinearAlgebra {


RNG::shared_ptr Vector::d_DefaultRNG;
#define DESCRIPTOR_ID_ARRAY_SCRATCH_SPACE ( 10 )


/****************************************************************
* Constructors                                                  *
****************************************************************/
Vector::Vector()
{
    d_Ghosts    = AMP::shared_ptr<std::vector<double>>( new std::vector<double> );
    d_AddBuffer = AMP::shared_ptr<std::vector<double>>( new std::vector<double> );
    d_UpdateState.reset( new UpdateState );
    *d_UpdateState = UNCHANGED;
    d_Views        = AMP::shared_ptr<std::vector<AMP::weak_ptr<Vector>>>(
        new std::vector<AMP::weak_ptr<Vector>>() );
}
Vector::Vector( VectorParameters::shared_ptr parameters )
{
    // Set default output stream
    d_output_stream = &AMP::plog;
    // Copy the relavent parameters
    AMP_INSIST( parameters->d_CommList.get() != nullptr,
                "d_CommList must be set in VectorParameters" );
    AMP_INSIST( parameters->d_DOFManager.get() != nullptr,
                "d_DOFManager must be set in VectorParameters" );
    setCommunicationList( parameters->d_CommList );
    d_UpdateState.reset( new UpdateState );
    *d_UpdateState = UNCHANGED;
    d_DOFManager   = parameters->d_DOFManager;
    d_Views        = AMP::shared_ptr<std::vector<AMP::weak_ptr<Vector>>>(
        new std::vector<AMP::weak_ptr<Vector>>() );
}


/****************************************************************
* De-Constructors                                               *
****************************************************************/
Vector::~Vector() {}


/****************************************************************
* Subset, View, and Select                                      *
****************************************************************/
Vector::shared_ptr Vector::selectInto( const VectorSelector &s )
{
    Vector::shared_ptr subvector;
    if ( s.isSelected( shared_from_this() ) ) {
        // Subset the vector
        subvector = s.subset( shared_from_this() );
        if ( subvector != nullptr ) {
            // Check the global size of the new vector to make sure it is <= the current size
            size_t N1 = this->getGlobalSize();
            size_t N2 = subvector->getGlobalSize();
            AMP_ASSERT( N2 <= N1 );
        }
    }
    return subvector;
}
Vector::const_shared_ptr Vector::selectInto( const VectorSelector &s ) const
{
    Vector::const_shared_ptr subvector;
    if ( s.isSelected( shared_from_this() ) ) {
        // Subset the vector
        subvector = s.subset( shared_from_this() );
        if ( subvector != nullptr ) {
            // Check the global size of the new vector to make sure it is <= the current size
            size_t N1 = this->getGlobalSize();
            size_t N2 = subvector->getGlobalSize();
            AMP_ASSERT( N2 <= N1 );
        }
    }
    return subvector;
}
Vector::shared_ptr Vector::select( const VectorSelector &s, const std::string &variable_name )
{
    if ( dynamic_cast<const VS_ByVariableName *>( &s ) ) {
        std::string name = dynamic_cast<const VS_ByVariableName *>( &s )->getName();
        if ( name == this->getVariable()->getName() )
            return shared_from_this();
    }
    Vector::shared_ptr retVal = this->selectInto( s );
    if ( retVal != nullptr ) {
        if ( AMP::dynamic_pointer_cast<MultiVector>( retVal ) == nullptr )
            retVal = MultiVector::view( retVal, retVal->getComm() );
        Variable::shared_ptr var( new Variable( variable_name ) );
        retVal->setVariable( var );
    }
    return retVal;
}
Vector::const_shared_ptr Vector::constSelect( const VectorSelector &s,
                                              const std::string &variable_name ) const
{
    if ( dynamic_cast<const VS_ByVariableName *>( &s ) ) {
        std::string name = dynamic_cast<const VS_ByVariableName *>( &s )->getName();
        if ( name == this->getVariable()->getName() )
            return shared_from_this();
    }
    Vector::const_shared_ptr retVal = this->selectInto( s );
    if ( retVal != nullptr ) {
        if ( AMP::dynamic_pointer_cast<const MultiVector>( retVal ) == nullptr )
            retVal = MultiVector::view( retVal, retVal->getComm() );
        Variable::shared_ptr var( new Variable( variable_name ) );
        AMP::const_pointer_cast<Vector>( retVal )->setVariable( var );
    }
    return retVal;
}
void Vector::registerView( Vector::shared_ptr v ) const
{
    for ( size_t i = 0; i != d_Views->size(); i++ )
        if ( ( *d_Views )[i].lock() == v )
            return;
    ( *d_Views ).push_back( v );
}
Vector::shared_ptr Vector::subsetVectorForVariable( const Variable::shared_ptr &name )
{
    Vector::shared_ptr retVal;
    if ( d_pVariable ) { // If there is a variable...
        if ( *d_pVariable == *name )
            retVal = shared_from_this();
    }
    return retVal;
}
Vector::const_shared_ptr
Vector::constSubsetVectorForVariable( const Variable::shared_ptr &name ) const
{
    Vector::const_shared_ptr retVal;
    if ( d_pVariable ) { // If there is a variable...
        if ( *d_pVariable == *name )
            retVal = shared_from_this();
    }
    return retVal;
}
Vector::shared_ptr Vector::cloneVector( const std::string &name ) const
{
    Vector::shared_ptr retVal;
    if ( getVariable() ) {
        retVal = cloneVector( getVariable()->cloneVariable( name ) );
    } else {
        retVal = cloneVector( Variable::shared_ptr( new Variable( name ) ) );
    }
    return retVal;
}


/****************************************************************
* Set/Get individual values                                     *
****************************************************************/
void Vector::setValuesByGlobalID( int numVals, size_t *ndx, const double *vals )
{
    INCREMENT_COUNT( "Virtual" );
    AMP_ASSERT( *d_UpdateState != ADDING );
    *d_UpdateState = SETTING;
    for ( int i = 0; i < numVals; i++ ) {
        if ( ( ndx[i] < getLocalStartID() ) ||
             ( ndx[i] >= ( getLocalStartID() + getLocalMaxID() ) ) ) {
            ( *d_Ghosts )[d_CommList->getLocalGhostID( ndx[i] )] = vals[i];
        } else {
            setLocalValuesByGlobalID( 1, ndx + i, vals + i );
        }
    }
}
void Vector::setGhostValuesByGlobalID( int numVals, size_t *ndx, const double *vals )
{
    INCREMENT_COUNT( "Virtual" );
    AMP_ASSERT( *d_UpdateState != ADDING );
    *d_UpdateState = SETTING;
    for ( int i = 0; i < numVals; i++ ) {
        if ( ( ndx[i] < getLocalStartID() ) ||
             ( ndx[i] >= ( getLocalStartID() + getLocalMaxID() ) ) ) {
            ( *d_Ghosts )[d_CommList->getLocalGhostID( ndx[i] )] = vals[i];
        } else {
            AMP_ERROR( "Non ghost index" );
        }
    }
}
void Vector::addValuesByGlobalID( int numVals, size_t *ndx, const double *vals )
{
    INCREMENT_COUNT( "Virtual" );
    AMP_ASSERT( *d_UpdateState != SETTING );
    *d_UpdateState = ADDING;
    for ( int i = 0; i < numVals; i++ ) {
        if ( ( ndx[i] < getLocalStartID() ) ||
             ( ndx[i] >= ( getLocalStartID() + getLocalMaxID() ) ) ) {
            ( *d_AddBuffer )[d_CommList->getLocalGhostID( ndx[i] )] += vals[i];
        } else {
            addLocalValuesByGlobalID( 1, ndx + i, vals + i );
        }
    }
}
void Vector::getValuesByLocalID( int num, size_t *ndx, double *vals ) const
{
    INCREMENT_COUNT( "Virtual" );
    for ( int i = 0; i != num; i++ ) {
        size_t block_number = 0;
        size_t offset       = ndx[i];
        while ( offset >= sizeOfDataBlock( block_number ) ) {
            offset -= sizeOfDataBlock( block_number );
            block_number++;
            if ( block_number >= numberOfDataBlocks() ) {
                AMP_ERROR( "Bad local id!" );
            }
        }
        vals[i] = getRawDataBlock<double>( block_number )[offset];
    }
}
void Vector::getValuesByGlobalID( int numVals, size_t *ndx, double *vals ) const
{
    INCREMENT_COUNT( "Virtual" );
    for ( int i = 0; i < numVals; i++ ) {
        if ( ( ndx[i] < getLocalStartID() ) ||
             ( ndx[i] >= ( getLocalStartID() + getLocalMaxID() ) ) ) {
            getGhostValuesByGlobalID( 1, ndx + i, vals + i );
        } else {
            getLocalValuesByGlobalID( 1, ndx + i, vals + i );
        }
    }
}
void Vector::getGhostValuesByGlobalID( int numVals, size_t *ndx, double *vals ) const
{
    INCREMENT_COUNT( "Virtual" );
    for ( int i = 0; i < numVals; i++ ) {
        if ( ( ndx[i] < getLocalStartID() ) ||
             ( ndx[i] >= ( getLocalStartID() + getLocalMaxID() ) ) ) {
            vals[i] = ( *d_Ghosts )[d_CommList->getLocalGhostID( ndx[i] )] +
                      ( *d_AddBuffer )[d_CommList->getLocalGhostID( ndx[i] )];
        } else {
            AMP_ERROR( "Tried to get a non-ghost ghost value" );
        }
    }
}
void Vector::getGhostAddValuesByGlobalID( int numVals, size_t *ndx, double *vals ) const
{
    INCREMENT_COUNT( "Virtual" );
    for ( int i = 0; i < numVals; i++ ) {
        if ( ( ndx[i] < getLocalStartID() ) ||
             ( ndx[i] >= ( getLocalStartID() + getLocalMaxID() ) ) ) {
            vals[i] = ( *d_AddBuffer )[d_CommList->getLocalGhostID( ndx[i] )];
        } else {
            AMP_ERROR( "Tried to get a non-ghost ghost value" );
        }
    }
}


/****************************************************************
* Functions to initalize the data                               *
****************************************************************/
void Vector::zero()
{
    setToScalar( 0.0 );
    for ( size_t i       = 0; i != d_Ghosts->size(); i++ )
        ( *d_Ghosts )[i] = 0.0;
}
void Vector::setToScalar( double alpha )
{
    iterator curMe = begin();
    iterator last  = end();
    while ( curMe != last ) {
        *curMe = alpha;
        ++curMe;
    }
    dataChanged();
    for ( size_t i            = 0; i != d_Ghosts->size(); i++ )
        ( *d_Ghosts )[i]      = alpha;
    ( *getUpdateStatusPtr() ) = UNCHANGED;
}
void Vector::setRandomValues()
{
    RandomVariable<double> r( 0., 1., getDefaultRNG() );
    iterator curMe = begin();
    iterator last  = end();
    while ( curMe != last ) {
        double curRand = r;
        *curMe         = curRand;
        ++curMe;
    }
    dataChanged();
    this->makeConsistent( CONSISTENT_SET );
}
void Vector::setRandomValues( RNG::shared_ptr rng )
{
    RandomVariable<double> r( 0., 1., rng );
    iterator curMe = begin();
    iterator last  = end();
    while ( curMe != last ) {
        *curMe = r;
        ++curMe;
    }
    dataChanged();
    this->makeConsistent( CONSISTENT_SET );
}


/****************************************************************
* Basic linear algebra                                          *
****************************************************************/
void Vector::copyVector( Vector::const_shared_ptr rhs )
{
    if ( rhs->getLocalSize() != getLocalSize() )
        AMP_ERROR( "Destination vector and source vector not the same size" );
    VectorDataIterator cur1      = begin();
    VectorDataIterator end1      = end();
    ConstVectorDataIterator cur2 = rhs->begin();
    while ( cur1 != end1 ) {
        *cur1 = *cur2;
        ++cur1;
        ++cur2;
    }
    if ( isA<DataChangeFirer>() )
        castTo<DataChangeFirer>().fireDataChange();
    copyGhostValues( rhs );
    // Copy the consistency state from the rhs
    *d_UpdateState = *( rhs->getUpdateStatusPtr() );
}
void Vector::scale( double alpha )
{
    iterator curMe = begin();
    iterator last  = end();
    while ( curMe != last ) {
        *curMe *= alpha;
        ++curMe;
    }
    dataChanged();
}
void Vector::scale( double alpha, const VectorOperations &x )
{
    const Vector &t_c = x.castTo<Vector>();
    AMP_ASSERT( getLocalSize() == t_c.getLocalSize() );
    iterator curMe        = begin();
    iterator last         = end();
    const_iterator curRhs = t_c.begin();
    while ( curMe != last ) {
        *curMe = alpha * *curRhs;
        ++curRhs;
        ++curMe;
    }
    dataChanged();
}
void Vector::add( const VectorOperations &x, const VectorOperations &y )
{
    const Vector &t_x = x.castTo<const Vector>();
    const Vector &t_y = y.castTo<const Vector>();
    AMP_ASSERT( getLocalSize() == t_x.getLocalSize() );
    AMP_ASSERT( getLocalSize() == t_y.getLocalSize() );
    iterator curMe         = begin();
    iterator last          = end();
    const_iterator curXRhs = t_x.begin();
    const_iterator curYRhs = t_y.begin();
    while ( curMe != last ) {
        *curMe = *curXRhs + *curYRhs;
        ++curXRhs;
        ++curYRhs;
        ++curMe;
    }
    dataChanged();
}
void Vector::subtract( const VectorOperations &x, const VectorOperations &y )
{
    const Vector &t_x = x.castTo<const Vector>();
    const Vector &t_y = y.castTo<const Vector>();
    AMP_ASSERT( getLocalSize() == t_x.getLocalSize() );
    AMP_ASSERT( getLocalSize() == t_y.getLocalSize() );
    iterator curMe         = begin();
    iterator last          = end();
    const_iterator curXRhs = t_x.begin();
    const_iterator curYRhs = t_y.begin();
    while ( curMe != last ) {
        *curMe = *curXRhs - *curYRhs;
        ++curXRhs;
        ++curYRhs;
        ++curMe;
    }
    dataChanged();
}
void Vector::multiply( const VectorOperations &x, const VectorOperations &y )
{
    const Vector &t_x = x.castTo<const Vector>();
    const Vector &t_y = y.castTo<const Vector>();
    AMP_ASSERT( getLocalSize() == t_x.getLocalSize() );
    AMP_ASSERT( getLocalSize() == t_y.getLocalSize() );
    iterator curMe         = begin();
    iterator last          = end();
    const_iterator curXRhs = t_x.begin();
    const_iterator curYRhs = t_y.begin();
    while ( curMe != last ) {
        *curMe = *curXRhs * *curYRhs;
        ++curXRhs;
        ++curYRhs;
        ++curMe;
    }
    dataChanged();
}
void Vector::divide( const VectorOperations &x, const VectorOperations &y )
{
    const Vector &t_x = x.castTo<const Vector>();
    const Vector &t_y = y.castTo<const Vector>();
    AMP_ASSERT( getLocalSize() == t_x.getLocalSize() );
    AMP_ASSERT( getLocalSize() == t_y.getLocalSize() );
    iterator curMe         = begin();
    iterator last          = end();
    const_iterator curXRhs = t_x.begin();
    const_iterator curYRhs = t_y.begin();
    while ( curMe != last ) {
        *curMe = *curXRhs / *curYRhs;
        ++curXRhs;
        ++curYRhs;
        ++curMe;
    }
    dataChanged();
}
void Vector::reciprocal( const VectorOperations &x )
{
    const Vector &t_c = x.castTo<Vector>();
    AMP_ASSERT( getLocalSize() == t_c.getLocalSize() );
    iterator curMe        = begin();
    iterator last         = end();
    const_iterator curRhs = t_c.begin();
    while ( curMe != last ) {
        *curMe = 1. / *curRhs;
        ++curRhs;
        ++curMe;
    }
    dataChanged();
}
void Vector::linearSum( double alpha,
                        const VectorOperations &x,
                        double beta,
                        const VectorOperations &y )
{
    const Vector &t_x = x.castTo<const Vector>();
    const Vector &t_y = y.castTo<const Vector>();
    AMP_ASSERT( getLocalSize() == t_x.getLocalSize() );
    AMP_ASSERT( getLocalSize() == t_y.getLocalSize() );
    iterator curMe         = begin();
    iterator last          = end();
    const_iterator curXRhs = t_x.begin();
    const_iterator curYRhs = t_y.begin();
    while ( curMe != last ) {
        *curMe = alpha * *curXRhs + beta * *curYRhs;
        ++curXRhs;
        ++curYRhs;
        ++curMe;
    }
    dataChanged();
}
void Vector::axpy( double alpha, const VectorOperations &x, const VectorOperations &y )
{
    const Vector &t_x = x.castTo<const Vector>();
    const Vector &t_y = y.castTo<const Vector>();
    AMP_ASSERT( getLocalSize() == t_x.getLocalSize() );
    AMP_ASSERT( getLocalSize() == t_y.getLocalSize() );
    iterator curMe         = begin();
    iterator last          = end();
    const_iterator curXRhs = t_x.begin();
    const_iterator curYRhs = t_y.begin();
    while ( curMe != last ) {
        *curMe = alpha * *curXRhs + *curYRhs;
        ++curXRhs;
        ++curYRhs;
        ++curMe;
    }
    dataChanged();
}
void Vector::axpby( double alpha, double beta, const VectorOperations &x )
{
    const Vector &t_x = x.castTo<const Vector>();
    AMP_ASSERT( getLocalSize() == t_x.getLocalSize() );
    iterator curMe         = begin();
    iterator last          = end();
    const_iterator curXRhs = t_x.begin();
    while ( curMe != last ) {
        *curMe = alpha * *curXRhs + beta * *curMe;
        ++curXRhs;
        ++curMe;
    }
    dataChanged();
}
void Vector::abs( const VectorOperations &x )
{
    const Vector &t_x = x.castTo<const Vector>();
    AMP_ASSERT( getLocalSize() == t_x.getLocalSize() );
    iterator curMe         = begin();
    iterator last          = end();
    const_iterator curXRhs = t_x.begin();
    while ( curMe != last ) {
        *curMe = fabs( *curXRhs );
        ++curXRhs;
        ++curMe;
    }
    dataChanged();
}


/****************************************************************
* min, max, norms, etc.                                         *
* Note: these routines require communication                    *
****************************************************************/
double Vector::min() const
{
    double ans = localMin();
    if ( getCommunicationList() )
        ans = getComm().minReduce( ans );
    return ans;
}
double Vector::max() const
{
    double ans = localMax();
    if ( getCommunicationList() )
        ans = getComm().maxReduce( ans );
    return ans;
}
double Vector::dot( const VectorOperations &x ) const
{
    AMP::shared_ptr<const Vector> vec = x.castTo<const Vector>().shared_from_this();
    double ans                        = localDot( vec );
    if ( getCommunicationList() )
        ans = getComm().sumReduce( ans );
    return ans;
}
double Vector::L1Norm( void ) const
{
    double ans = localL1Norm();
    if ( getCommunicationList() )
        ans = getComm().sumReduce( ans );
    return ans;
}
double Vector::maxNorm( void ) const
{
    double ans = localMaxNorm();
    if ( getCommunicationList() )
        ans = getComm().maxReduce( ans );
    return ans;
}
double Vector::L2Norm( void ) const
{
    double ans = localL2Norm();
    if ( getCommunicationList() )
        ans = sqrt( getComm().sumReduce( ans * ans ) );
    return ans;
}
double Vector::localMin( void ) const
{
    size_t N_blocks = numberOfDataBlocks();
    double ans      = 1e300;
    for ( size_t i = 0; i < N_blocks; i++ ) {
        size_t size        = sizeOfDataBlock( i );
        const double *data = getRawDataBlock<double>( i );
        for ( size_t j = 0; j < size; j++ )
            ans = std::min( data[j], ans );
    }
    return ans;
}
double Vector::localMax( void ) const
{
    size_t N_blocks = numberOfDataBlocks();
    double ans      = -1e300;
    for ( size_t i = 0; i < N_blocks; i++ ) {
        size_t size        = sizeOfDataBlock( i );
        const double *data = getRawDataBlock<double>( i );
        for ( size_t j = 0; j < size; j++ )
            ans = std::max( data[j], ans );
    }
    return ans;
}
double Vector::localL1Norm( void ) const
{
    size_t N_blocks = numberOfDataBlocks();
    double ans      = 0.0;
    for ( size_t i = 0; i < N_blocks; i++ ) {
        size_t size        = sizeOfDataBlock( i );
        const double *data = getRawDataBlock<double>( i );
        for ( size_t j = 0; j < size; j++ )
            ans += fabs( data[j] );
    }
    return ans;
}
double Vector::localL2Norm( void ) const
{
    size_t N_blocks = numberOfDataBlocks();
    double ans      = 0.0;
    for ( size_t i = 0; i < N_blocks; i++ ) {
        size_t size        = sizeOfDataBlock( i );
        const double *data = getRawDataBlock<double>( i );
        for ( size_t j = 0; j < size; j++ )
            ans += data[j] * data[j];
    }
    return sqrt( ans );
}
double Vector::localMaxNorm( void ) const
{
    size_t N_blocks = numberOfDataBlocks();
    double ans      = 0.0;
    for ( size_t i = 0; i < N_blocks; i++ ) {
        size_t size        = sizeOfDataBlock( i );
        const double *data = getRawDataBlock<double>( i );
        for ( size_t j = 0; j < size; j++ )
            ans = std::max( fabs( data[j] ), ans );
    }
    return ans;
}
double Vector::localDot( AMP::shared_ptr<const Vector> x ) const
{
    AMP_ASSERT( getLocalSize() == x->getLocalSize() );
    const_iterator curMe   = begin();
    const_iterator last    = end();
    const_iterator curXRhs = x->begin();
    double ans             = 0.0;
    while ( curMe != last ) {
        ans += *curMe * *curXRhs;
        ++curXRhs;
        ++curMe;
    }
    return ans;
}


void Vector::setCommunicationList( CommunicationList::shared_ptr comm )
{
    AMP_ASSERT( comm );
    d_CommList = comm;
    if ( comm ) {
        addCommunicationListToParameters( comm );
        d_Ghosts = AMP::shared_ptr<std::vector<double>>(
            new std::vector<double>( d_CommList->getVectorReceiveBufferSize() ) );
        d_AddBuffer = AMP::shared_ptr<std::vector<double>>(
            new std::vector<double>( d_CommList->getVectorReceiveBufferSize() ) );
    }
}


bool Vector::equals( Vector const &rhs, double tol ) const
{
    int RetVal = 0;
    if ( ( getGlobalSize() == rhs.getGlobalSize() ) && ( getLocalSize() == rhs.getLocalSize() ) ) {
        ConstVectorDataIterator cur1 = begin();
        ConstVectorDataIterator cur2 = rhs.begin();
        ConstVectorDataIterator last = end();
        bool failed                  = false;
        while ( cur1 != last ) {
            double v1 = *cur1;
            double v2 = *cur2;
            if ( fabs( v1 - v2 ) > tol ) {
                failed = true;
                break;
            }
            ++cur1;
            ++cur2;
        }
        if ( !failed ) {
            RetVal = 1;
        }
    }

    int ans;
    if ( d_CommList ) {
        ans = d_CommList->getComm().minReduce( RetVal );
    } else {
        ans = RetVal;
    }

    return ans == 1 ? true : false;
}


// The following two functions are Jungho's
// This will fail when y_i = 0... Needs to be fixed
double Vector::minQuotient( const VectorOperations &x, const VectorOperations &y )
{
    const Vector &x_vec         = x.castTo<Vector>();
    const Vector &y_vec         = y.castTo<Vector>();
    Vector::const_iterator curx = x_vec.begin();
    Vector::const_iterator cury = y_vec.begin();
    Vector::const_iterator endx = x_vec.end();
    while ( curx != endx ) {
        if ( *cury != 0.0 )
            break;
        ++curx;
        ++cury;
    }
    // Probably should do a test across processors, not just the one
    AMP_INSIST( curx != endx, "denominator is the zero vector on an entire process" );
    double myRetVal = ( *curx ) / ( *cury );
    while ( curx != endx ) {
        if ( *cury != 0.0 ) {
            myRetVal = std::min( myRetVal, ( *curx ) / ( *cury ) );
        }
        ++curx;
        ++cury;
    }
    double retVal = x_vec.getComm().minReduce( myRetVal );
    return retVal;
}


double Vector::wrmsNorm( const VectorOperations &x, const VectorOperations &y )
{
    double dot_prod             = 0.0;
    int global_size             = 0;
    Vector::shared_ptr temp_vec = x.castTo<Vector>().cloneVector();
    temp_vec->multiply( x, y );

    dot_prod    = temp_vec->dot( temp_vec );
    global_size = temp_vec->getGlobalSize();

    return ( sqrt( dot_prod / global_size ) );
}


double Vector::wrmsNormMask( const VectorOperations &x,
                             const VectorOperations &y,
                             const VectorOperations &mask )
{
    double dot_prod     = 0.0;
    const Vector &x_vec = x.castTo<Vector>();
    const Vector &y_vec = y.castTo<Vector>();
    const Vector &m_vec = mask.castTo<Vector>();

    Vector::const_iterator curx = x_vec.begin();
    Vector::const_iterator endx = x_vec.end();
    Vector::const_iterator cury = y_vec.begin();
    Vector::const_iterator curm = m_vec.begin();
    while ( curx != endx ) {
        if ( *curm > 0.0 ) {
            dot_prod += ( *curx ) * ( *curx ) * ( *cury ) * ( *cury );
        }
        ++curx;
        ++cury;
        ++curm;
    }
    AMP_ASSERT( cury == y_vec.end() );
    AMP_ASSERT( curm == m_vec.end() );
    double all_dot = x_vec.getComm().sumReduce( dot_prod );
    return sqrt( all_dot / (double) x_vec.getGlobalSize() );
}


void Vector::makeConsistent( ScatterType t )
{
    if ( t == CONSISTENT_ADD ) {
        AMP_ASSERT( *d_UpdateState != SETTING );
        std::vector<double> send_vec_add( d_CommList->getVectorReceiveBufferSize() );
        std::vector<double> recv_vec_add( d_CommList->getVectorSendBufferSize() );
        d_CommList->packReceiveBuffer( send_vec_add, *this );
        d_CommList->scatter_add( send_vec_add, recv_vec_add );
        d_CommList->unpackSendBufferAdd( recv_vec_add, *this );
        for ( auto &elem : *d_AddBuffer ) {
            elem = 0.0;
        }
    }
    *d_UpdateState = SETTING;
    std::vector<double> send_vec( d_CommList->getVectorSendBufferSize() );
    std::vector<double> recv_vec( d_CommList->getVectorReceiveBufferSize() );
    d_CommList->packSendBuffer( send_vec, *this );
    d_CommList->scatter_set( send_vec, recv_vec );
    d_CommList->unpackReceiveBufferSet( recv_vec, *this );
    *d_UpdateState = UNCHANGED;
    this->setUpdateStatus( UNCHANGED );
}


void Vector::copyGhostValues( const AMP::shared_ptr<const Vector> &rhs )
{
    if ( getGhostSize() == 0 ) {
        // No ghosts to fill, we don't need to do anything
    } else if ( getGhostSize() == rhs->getGhostSize() ) {
        // The ghosts in the src vector match the current vector
        // Copy the ghosts from the rhs
        std::vector<size_t> ghostIDs = getCommunicationList()->getGhostIDList();
        std::vector<double> values( ghostIDs.size() );
        rhs->getGhostValuesByGlobalID( ghostIDs.size(), &ghostIDs[0], &values[0] );
        this->setGhostValuesByGlobalID( ghostIDs.size(), &ghostIDs[0], &values[0] );
        // Copy the consistency state from the rhs
        *d_UpdateState = *( rhs->getUpdateStatusPtr() );
    } else {
        // We can't copy the ghosts from the rhs
        // Use makeConsistent to fill the ghosts
        // Note: this will incure global communication
        makeConsistent( CONSISTENT_SET );
    }
}


void Vector::dumpOwnedData( std::ostream &out, size_t GIDoffset, size_t LIDoffset ) const
{
    const_iterator curElement = begin();
    size_t gid                = GIDoffset;
    if ( getCommunicationList() )
        gid += getCommunicationList()->getStartGID();
    size_t lid = LIDoffset;
    while ( curElement != end() ) {
        out << "  GID: " << gid << "  LID: " << lid << "  Value: " << *curElement << "\n";
        ++curElement;
        ++gid;
        ++lid;
    }
}


void Vector::dumpGhostedData( std::ostream &out, size_t offset ) const
{
    if ( !getCommunicationList() )
        return;
    const std::vector<size_t> &ghosts = getCommunicationList()->getGhostIDList();
    auto curVal                       = d_Ghosts->begin();
    for ( auto &ghost : ghosts ) {
        out << "  GID: " << ( ghost + offset ) << "  Value: " << ( *curVal ) << "\n";
        ++curVal;
    }
}


std::ostream &operator<<( std::ostream &out, const Vector &v )
{
    out << "Vector type: " << v.type() << "\n";
    if ( v.getVariable() ) {
        out << "Variable name: " << v.getVariable()->getName() << "\n";
    }
    if ( v.getCommunicationList() ) {
        int rank = v.getComm().getRank();
        out << "Processor: " << rank << "\n";
    }
    out << "\n"
        << "Number of owned elements: " << v.getLocalSize() << "\n"
        << "Number of ghosted elements: " << v.getGhostSize() << "\n";
    if ( v.numberOfDataBlocks() > 1 )
        out << "Number of sub-vectors: " << v.numberOfDataBlocks() << "\n";
    out << "\n"
        << "Data block pointers: \n";

    for ( size_t i = 0; i != v.numberOfDataBlocks(); i++ )
        out << "                     " << v.getRawDataBlock<double>( i ) << "  ( "
            << v.sizeOfDataBlock( i ) << " elements )\n";

    out << "\n"
        << "Local Data:\n";
    v.dumpOwnedData( out );

    out << "\n"
        << "Ghosted Data:\n";
    v.dumpGhostedData( out );

    return out;
}
}
}
