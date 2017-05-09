#include "utils/Utilities.h"

#include "EpetraVectorEngine.h"

#ifdef USE_EXT_MPI
#include <Epetra_MpiComm.h>
#else
#include <Epetra_SerialComm.h>
#endif

namespace AMP {
namespace LinearAlgebra {


static inline double *getBufferPtr( VectorEngine::BufferPtr buf )
{
    if ( buf->empty() )
        return nullptr;
    return &buf->operator[]( 0 );
}


static inline Epetra_Vector& getEpetraVector( VectorOperations& vec )
{
    auto epetra = dynamic_cast<EpetraVectorEngine*>( &vec );
    AMP_INSIST( epetra != nullptr, "Not an EpetraVectorEngine" );
    return epetra->getEpetra_Vector();
}
static inline const Epetra_Vector& getEpetraVector( const VectorOperations& vec )
{
    auto epetra = dynamic_cast<const EpetraVectorEngine*>( &vec );
    AMP_INSIST( epetra != nullptr, "Not an EpetraVectorEngine" );
    return epetra->getEpetra_Vector();
}


/********************************************************
* EpetraVectorEngineParameters constructors             *
********************************************************/
EpetraVectorEngineParameters::EpetraVectorEngineParameters( size_t local_size,
                                                            size_t global_size,
                                                            AMP_MPI comm )
    : VectorEngineParameters( local_size, global_size, comm )
{
}
EpetraVectorEngineParameters::EpetraVectorEngineParameters(
    size_t local_size, size_t global_size, AMP::shared_ptr<Epetra_Map> emap, AMP_MPI ecomm )
    : VectorEngineParameters( local_size, global_size, ecomm ), d_emap( emap )
{
}
EpetraVectorEngineParameters::~EpetraVectorEngineParameters() {}


/********************************************************
* Function to return (and create) the Epetra_Map        *
********************************************************/
Epetra_Map &EpetraVectorEngineParameters::getEpetraMap()
{
    if ( d_emap.get() != nullptr )
        return *d_emap;
// Create the epetra map
#ifdef USE_EXT_MPI
    Epetra_MpiComm comm = d_comm.getCommunicator();
#else
    Epetra_SerialComm comm;
#endif
    AMP_INSIST( d_global < 0x80000000,
                "Epetra does not support vectors with global size greater than 2^31" );
    size_t local_size = d_end - d_begin;
    // std::vector<int> ids(local_size,0);
    // for (size_t i=0; i<local_size; i++)
    //    ids[i] = (int) (i+d_begin);
    // d_emap = AMP::shared_ptr<Epetra_Map> ( new Epetra_Map ( (int) d_global, (int) local_size,
    // &ids[0], 0, comm ) );
    d_emap =
        AMP::shared_ptr<Epetra_Map>( new Epetra_Map( (int) d_global, (int) local_size, 0, comm ) );
    // Check the map to make sure it is correct
    AMP_ASSERT( local_size == (size_t) d_emap->NumMyPoints() );
    AMP_ASSERT( d_global == (size_t) d_emap->NumGlobalPoints() );
    AMP_ASSERT( d_begin == (size_t) d_emap->MinMyGID() );
    AMP_ASSERT( d_end - 1 == (size_t) d_emap->MaxMyGID() );
    AMP_ASSERT( 0 == (size_t) d_emap->MinAllGID() );
    AMP_ASSERT( d_global - 1 == (size_t) d_emap->MaxAllGID() );
    AMP_ASSERT( 0 == (size_t) d_emap->MinLID() );
    if ( local_size == 0 )
        AMP_ASSERT( 0 == (size_t) d_emap->MaxLID() );
    else
        AMP_ASSERT( local_size - 1 == (size_t) d_emap->MaxLID() );
    return *d_emap;
}


/********************************************************
* Constructor                                           *
********************************************************/
EpetraVectorEngine::EpetraVectorEngine( VectorEngineParameters::shared_ptr alias, BufferPtr buf )
    : d_epetraVector(
          View, 
          dynamic_pointer_cast<EpetraVectorEngineParameters>( alias )->getEpetraMap(), 
          getBufferPtr( buf ) ),
      d_iLocalSize( 0 ),
      d_iGlobalSize( 0 )
{
    d_Params = alias;
}


VectorEngine::BufferPtr EpetraVectorEngine::getNewBuffer()
{
    BufferPtr retval( new std::vector<double>( getLocalSize() ) );
    //    d_epetraVector.ResetView ( &*(retval->begin()) );
    return retval;
}

void EpetraVectorEngine::swapEngines( VectorEngine::shared_ptr x )
{
    double *my_pointer;
    double *oth_pointer;
    getEpetra_Vector().ExtractView( &my_pointer );
    getEpetraVector( *x ).ExtractView( &oth_pointer );
    getEpetraVector( *x ).ResetView( my_pointer );
    getEpetra_Vector().ResetView( oth_pointer );
}

void EpetraVectorEngine::setToScalar( const double alpha )
{
    getEpetra_Vector().PutScalar( alpha );
}

void EpetraVectorEngine::scale( double alpha, const VectorOperations &x )
{
    getEpetra_Vector().Scale( alpha, getEpetraVector( x ) );
}

void EpetraVectorEngine::scale( double alpha ) { getEpetra_Vector().Scale( alpha ); }

VectorEngine::shared_ptr EpetraVectorEngine::cloneEngine( BufferPtr p ) const
{
    return shared_ptr( new EpetraVectorEngine( d_Params, p ) );
}

void EpetraVectorEngine::add( const VectorOperations &x, const VectorOperations &y )
{
    getEpetra_Vector().Update( 1., getEpetraVector( x ), 1., getEpetraVector( y ),
                               0. );
}

void EpetraVectorEngine::subtract( const VectorOperations &x, const VectorOperations &y )
{
    getEpetra_Vector().Update( 1., getEpetraVector( x ), -1., getEpetraVector( y ), 0. );
}

void EpetraVectorEngine::multiply( const VectorOperations &x, const VectorOperations &y )
{
    getEpetra_Vector().Multiply( 1., getEpetraVector( x ), getEpetraVector( y ), 0. );
}

void EpetraVectorEngine::divide( const VectorOperations &x, const VectorOperations &y )
{
    getEpetra_Vector().ReciprocalMultiply( 1., getEpetraVector( y ), getEpetraVector( x ), 0. );
}

void EpetraVectorEngine::reciprocal( const VectorOperations &x )
{
    getEpetra_Vector().Reciprocal( getEpetraVector( x ) );
}

void EpetraVectorEngine::linearSum( double alpha,
                                    const VectorOperations &x,
                                    double beta,
                                    const VectorOperations &y )
{
    getEpetra_Vector().Update( alpha, getEpetraVector( x ), beta, getEpetraVector( y ), 0. );
}

void EpetraVectorEngine::axpy( double alpha, const VectorOperations &x, const VectorOperations &y )
{
    linearSum( alpha, x, 1., y );
}

void EpetraVectorEngine::axpby( double alpha, double beta, const VectorOperations &x )
{
    linearSum( alpha, x, beta, *this );
}

void EpetraVectorEngine::abs( const VectorOperations &x )
{
    getEpetra_Vector().Abs( getEpetraVector( x ) );
}

double EpetraVectorEngine::min( void ) const
{
    double retVal;
    getEpetra_Vector().MinValue( &retVal );
    return retVal;
}

double EpetraVectorEngine::max( void ) const
{
    double retVal;
    getEpetra_Vector().MaxValue( &retVal );
    return retVal;
}

void EpetraVectorEngine::setRandomValues( void )
{
    getEpetra_Vector().Random();
    this->abs( *this );
}


void EpetraVectorEngine::setValuesByLocalID( int num, size_t *indices, const double *vals )
{
    INCREMENT_COUNT( "Virtual" );
    for ( int i                        = 0; i != num; i++ )
        getEpetra_Vector()[indices[i]] = vals[i];
}


void EpetraVectorEngine::setLocalValuesByGlobalID( int num, size_t *indices, const double *vals )
{
    INCREMENT_COUNT( "Virtual" );
    if ( num == 0 )
        return;
    AMP_ASSERT( getGlobalSize() < 0x80000000 );
    std::vector<int> indices2( num, 0 );
    for ( int i     = 0; i < num; i++ )
        indices2[i] = (int) indices[i];
    getEpetra_Vector().ReplaceGlobalValues( num, const_cast<double *>( vals ), &indices2[0] );
}

void EpetraVectorEngine::addValuesByLocalID( int num, size_t *indices, const double *vals )
{
    INCREMENT_COUNT( "Virtual" );
    if ( num == 0 )
        return;
    for ( int i = 0; i != num; i++ )
        getEpetra_Vector()[indices[i]] += vals[i];
}

void EpetraVectorEngine::addLocalValuesByGlobalID( int num, size_t *indices, const double *vals )
{
    INCREMENT_COUNT( "Virtual" );
    if ( num == 0 )
        return;
    AMP_ASSERT( getGlobalSize() < 0x80000000 );
    std::vector<int> indices2( num, 0 );
    for ( int i     = 0; i < num; i++ )
        indices2[i] = (int) indices[i];
    getEpetra_Vector().SumIntoGlobalValues( num, const_cast<double *>( vals ), &indices2[0] );
}

void EpetraVectorEngine::getValuesByLocalID( int, size_t *, double * ) const
{
    INCREMENT_COUNT( "Virtual" );
    AMP_ERROR( "This shouldn't be called" );
}

void EpetraVectorEngine::getLocalValuesByGlobalID( int, size_t *, double * ) const
{
    INCREMENT_COUNT( "Virtual" );
    AMP_ERROR( "This shouldn't be called" );
}
double EpetraVectorEngine::L1Norm( void ) const
{
    double retVal;
    getEpetra_Vector().Norm1( &retVal );
    return retVal;
}

double EpetraVectorEngine::L2Norm( void ) const
{
    double retVal;
    getEpetra_Vector().Norm2( &retVal );
    return retVal;
}


double EpetraVectorEngine::maxNorm( void ) const
{
    double retVal;
    getEpetra_Vector().NormInf( &retVal );
    return retVal;
}

double EpetraVectorEngine::dot( const VectorOperations &x ) const
{
    double retVal;
    getEpetra_Vector().Dot( getEpetraVector( x ), &retVal );
    return retVal;
}

double EpetraVectorEngine::localDot( const VectorOperations & ) const
{
    AMP_ERROR("Not implimented");
    return 0;
}

void EpetraVectorEngine::putRawData( const double *in )
{
    double *p;
    getEpetra_Vector().ExtractView( &p );
    size_t N = getLocalSize();
    memcpy( p, in, N * sizeof( double ) );
}

void EpetraVectorEngine::copyOutRawData( double *out ) const
{
    getEpetra_Vector().ExtractCopy( out );
}


double EpetraVectorEngine::localMin( void ) const
{
    AMP_ERROR("Not implimented");
    return 0;
}
double EpetraVectorEngine::localMax( void ) const
{
    AMP_ERROR("Not implimented");
    return 0;
}
double EpetraVectorEngine::localL1Norm( void ) const
{
    AMP_ERROR("Not implimented");
    return 0;
}
double EpetraVectorEngine::localL2Norm( void ) const
{
    AMP_ERROR("Not implimented");
    return 0;
}
double EpetraVectorEngine::localMaxNorm( void ) const
{
    AMP_ERROR("Not implimented");
    return 0;
}


}
}
