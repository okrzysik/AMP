#include "AMP/vectors/trilinos/epetra/EpetraVectorEngine.h"
#include "AMP/utils/Utilities.h"
#include "AMP/vectors/data/VectorDataCPU.h"

#ifdef USE_EXT_MPI
#include <Epetra_MpiComm.h>

#include <utility>

#else
#include <Epetra_SerialComm.h>
#endif

namespace AMP {
namespace LinearAlgebra {


static inline double *getBufferPtr( std::shared_ptr<VectorData> buf )
{
    size_t N_blocks = buf->numberOfDataBlocks();
    if ( N_blocks == 0 )
        return nullptr;
    if ( N_blocks > 1 )
        AMP_ERROR( "More than 1 data block detected" );
    return buf->getRawDataBlock<double>( 0 );
}


static inline Epetra_Vector &getEpetraVector( VectorOperations &vec )
{
    auto epetra = dynamic_cast<EpetraVectorEngine *>( &vec );
    AMP_INSIST( epetra != nullptr, "Not an EpetraVectorEngine" );
    return epetra->getEpetra_Vector();
}


/********************************************************
 * EpetraVectorEngineParameters constructors             *
 ********************************************************/
EpetraVectorEngineParameters::EpetraVectorEngineParameters( size_t local_size,
                                                            size_t global_size,
                                                            const AMP_MPI &comm )
    : VectorEngineParameters( local_size, global_size, comm )
{
}
EpetraVectorEngineParameters::EpetraVectorEngineParameters( size_t local_size,
                                                            size_t global_size,
                                                            std::shared_ptr<Epetra_Map> emap,
                                                            const AMP_MPI &ecomm )
    : VectorEngineParameters( local_size, global_size, ecomm ), d_emap( std::move( emap ) )
{
}
EpetraVectorEngineParameters::~EpetraVectorEngineParameters() = default;


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
    // d_emap = std::shared_ptr<Epetra_Map> ( new Epetra_Map ( (int) d_global, (int) local_size,
    // &ids[0], 0, comm ) );
    d_emap = std::make_shared<Epetra_Map>( (int) d_global, (int) local_size, 0, comm );
    // Check the map to make sure it is correct
    AMP_ASSERT( local_size == (size_t) d_emap->NumMyPoints() );
    AMP_ASSERT( d_global == (size_t) d_emap->NumGlobalPoints() );
    AMP_ASSERT( d_begin == (size_t) d_emap->MinMyGID() );
    AMP_ASSERT( d_end - 1 == (size_t) d_emap->MaxMyGID() );
    AMP_ASSERT( 0 == (size_t) d_emap->MinAllGID() );
    AMP_ASSERT( d_global - 1 == (size_t) d_emap->MaxAllGID() );
    AMP_ASSERT( 0 == (size_t) d_emap->MinLID() );
    AMP_ASSERT( d_emap->LinearMap() );
    if ( local_size == 0 )
        AMP_ASSERT( 0 == (size_t) d_emap->MaxLID() );
    else
        AMP_ASSERT( local_size - 1 == (size_t) d_emap->MaxLID() );
    return *d_emap;
}


/********************************************************
 * Constructor                                           *
 ********************************************************/
EpetraVectorEngine::EpetraVectorEngine( std::shared_ptr<VectorEngineParameters> alias,
                                        std::shared_ptr<VectorData> buf )
    : EpetraVectorData(
          View,
          std::dynamic_pointer_cast<EpetraVectorEngineParameters>( alias )->getEpetraMap(),
          getBufferPtr( buf ),
          std::dynamic_pointer_cast<EpetraVectorEngineParameters>( alias )->beginDOF(),
          std::dynamic_pointer_cast<EpetraVectorEngineParameters>( alias )->getLocalSize(),
          std::dynamic_pointer_cast<EpetraVectorEngineParameters>( alias )->getGlobalSize() )
{
    d_Params     = alias;
    d_buf_scope  = buf;
    d_VectorData = dynamic_cast<VectorData *>( this );
}

std::shared_ptr<VectorEngine> EpetraVectorEngine::cloneEngine( std::shared_ptr<VectorData> p ) const
{
  std::cout << "EpetraVectorEngine::cloneEngine with "
	    << " d_startIndex " << std::dynamic_pointer_cast<EpetraVectorEngineParameters>( d_Params )->beginDOF()
	    << " d_localSize " << std::dynamic_pointer_cast<EpetraVectorEngineParameters>( d_Params )->getLocalSize()
	    << " d_globalSize " << std::dynamic_pointer_cast<EpetraVectorEngineParameters>( d_Params )->getGlobalSize()
	    << std::endl;
    return std::make_shared<EpetraVectorEngine>( d_Params, p );
}

Vector::shared_ptr EpetraVectorEngine::cloneVector( const Variable::shared_ptr name ) const
{
    auto params = std::dynamic_pointer_cast<EpetraVectorEngineParameters>( d_Params );
    auto buffer = std::make_shared<VectorDataCPU<double>>(
        params->beginDOF(), params->getLocalSize(), params->getGlobalSize() );

    auto retVal = std::make_shared<EpetraVectorEngine>( d_Params, buffer );
    retVal->setVariable( name );
    return retVal;
}

void EpetraVectorEngine::swapVectors( Vector &other )
{
    double *my_pointer;
    double *oth_pointer;
    getEpetra_Vector().ExtractView( &my_pointer );
    getEpetraVector( other ).ExtractView( &oth_pointer );
    getEpetraVector( other ).ResetView( my_pointer );
    getEpetra_Vector().ResetView( oth_pointer );
}

void EpetraVectorEngine::aliasVector( Vector & )
{
    AMP_ERROR( "EpetraVectorEngine::aliasVector not implemented" );
}

void EpetraVectorEngine::assemble() {}

} // namespace LinearAlgebra
} // namespace AMP
