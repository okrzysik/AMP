#include "vectors/trilinos/epetra/EpetraVectorEngine.h"
#include "utils/Utilities.h"
#include "vectors/data/VectorDataCPU.h"

#ifdef USE_EXT_MPI
#include <Epetra_MpiComm.h>

#include <utility>

#else
#include <Epetra_SerialComm.h>
#endif

namespace AMP {
namespace LinearAlgebra {


static inline double *getBufferPtr( AMP::shared_ptr<VectorData> buf )
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
                                                            AMP_MPI comm )
    : VectorEngineParameters( local_size, global_size, comm )
{
}
EpetraVectorEngineParameters::EpetraVectorEngineParameters( size_t local_size,
                                                            size_t global_size,
                                                            AMP::shared_ptr<Epetra_Map> emap,
                                                            AMP_MPI ecomm )
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
    // d_emap = AMP::shared_ptr<Epetra_Map> ( new Epetra_Map ( (int) d_global, (int) local_size,
    // &ids[0], 0, comm ) );
    d_emap = AMP::make_shared<Epetra_Map>( (int) d_global, (int) local_size, 0, comm );
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
EpetraVectorEngine::EpetraVectorEngine( AMP::shared_ptr<VectorEngineParameters> alias,
                                        AMP::shared_ptr<VectorData> buf )
    : EpetraVectorData(
          View,
          dynamic_pointer_cast<EpetraVectorEngineParameters>( alias )->getEpetraMap(),
          getBufferPtr( buf ),
          dynamic_pointer_cast<EpetraVectorEngineParameters>( alias )->beginDOF(),
          dynamic_pointer_cast<EpetraVectorEngineParameters>( alias )->getLocalSize(),
          dynamic_pointer_cast<EpetraVectorEngineParameters>( alias )->getGlobalSize() )
{
    d_Params    = alias;
    d_buf_scope = buf;
}


AMP::shared_ptr<VectorData> EpetraVectorEngine::getNewBuffer()
{
    auto buffer = AMP::make_shared<VectorDataCPU<double>>(
        getLocalStartID(), getLocalSize(), getGlobalSize() );
    return buffer;
}

void EpetraVectorEngine::swapEngines( AMP::shared_ptr<VectorEngine> x )
{
    double *my_pointer;
    double *oth_pointer;
    getEpetra_Vector().ExtractView( &my_pointer );
    getEpetraVector( *x ).ExtractView( &oth_pointer );
    getEpetraVector( *x ).ResetView( my_pointer );
    getEpetra_Vector().ResetView( oth_pointer );
}

AMP::shared_ptr<VectorEngine> EpetraVectorEngine::cloneEngine( AMP::shared_ptr<VectorData> p ) const
{
    return AMP::make_shared<EpetraVectorEngine>( d_Params, p );
}


} // namespace LinearAlgebra
} // namespace AMP
