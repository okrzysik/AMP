#include "AMP/vectors/trilinos/epetra/NativeEpetraVector.h"
#include "AMP/vectors/Vector.h"
#include "AMP/vectors/data/VectorDataCPU.h"
#include "AMP/utils/Utilities.h"

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
    auto epetra = dynamic_cast<NativeEpetraVector *>( &vec );
    AMP_INSIST( epetra != nullptr, "Not an NativeEpetraVector" );
    return epetra->getEpetra_Vector();
}


/********************************************************
 * NativeEpetraVectorParameters constructors             *
 ********************************************************/
NativeEpetraVectorParameters::NativeEpetraVectorParameters( size_t local_size,
                                                            size_t global_size,
                                                            const AMP_MPI &comm )
  : VectorParameters(), d_local_size{local_size}, d_global_size{global_size}, d_comm{comm}
{
    d_comm.sumScan( &local_size, &d_end, 1 );
    d_begin = d_end - local_size;
}
NativeEpetraVectorParameters::NativeEpetraVectorParameters( size_t local_size,
                                                            size_t global_size,
                                                            std::shared_ptr<Epetra_Map> emap,
                                                            const AMP_MPI &ecomm )
  : VectorParameters(), d_local_size{local_size}, d_global_size{global_size}, d_comm{ecomm}, d_emap( std::move( emap ))
{
    d_comm.sumScan( &local_size, &d_end, 1 );
    d_begin = d_end - local_size;
}
												    
NativeEpetraVectorParameters::~NativeEpetraVectorParameters() = default;


/********************************************************
 * Function to return (and create) the Epetra_Map        *
 ********************************************************/
Epetra_Map &NativeEpetraVectorParameters::getEpetraMap()
{
    if ( d_emap.get() != nullptr )
        return *d_emap;
// Create the epetra map
#ifdef USE_EXT_MPI
    Epetra_MpiComm comm = d_comm.getCommunicator();
#else
    Epetra_SerialComm comm;
#endif
    AMP_INSIST( d_global_size < 0x80000000,
                "Epetra does not support vectors with global size greater than 2^31" );
    size_t local_size = d_end - d_begin;
    // std::vector<int> ids(local_size,0);
    // for (size_t i=0; i<local_size; i++)
    //    ids[i] = (int) (i+d_begin);
    // d_emap = std::shared_ptr<Epetra_Map> ( new Epetra_Map ( (int) d_global_size, (int) local_size,
    // &ids[0], 0, comm ) );
    d_emap = std::make_shared<Epetra_Map>( (int) d_global_size, (int) local_size, 0, comm );
    // Check the map to make sure it is correct
    AMP_ASSERT( local_size == (size_t) d_emap->NumMyPoints() );
    AMP_ASSERT( d_global_size == (size_t) d_emap->NumGlobalPoints() );
    AMP_ASSERT( d_begin == (size_t) d_emap->MinMyGID() );
    AMP_ASSERT( d_end - 1 == (size_t) d_emap->MaxMyGID() );
    AMP_ASSERT( 0 == (size_t) d_emap->MinAllGID() );
    AMP_ASSERT( d_global_size - 1 == (size_t) d_emap->MaxAllGID() );
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
NativeEpetraVector::NativeEpetraVector( std::shared_ptr<VectorParameters> alias,
                                        std::shared_ptr<VectorData> buf )
    : EpetraVectorData(
          View,
          std::dynamic_pointer_cast<NativeEpetraVectorParameters>( alias )->getEpetraMap(),
          getBufferPtr( buf ),
          std::dynamic_pointer_cast<NativeEpetraVectorParameters>( alias )->beginDOF(),
          std::dynamic_pointer_cast<NativeEpetraVectorParameters>( alias )->getLocalSize(),
          std::dynamic_pointer_cast<NativeEpetraVectorParameters>( alias )->getGlobalSize() ),
      d_Params{ alias}
{
    d_buf_scope = buf;
    CommunicationListParameters::shared_ptr params( new CommunicationListParameters() );
    params->d_comm      = getComm();
    auto local_size = std::dynamic_pointer_cast<NativeEpetraVectorParameters>( alias )->getLocalSize();
    params->d_localsize = local_size;
    setCommunicationList( std::make_shared<CommunicationList>( params ) );
    d_DOFManager = std::make_shared<AMP::Discretization::DOFManager>( local_size,
                                                                      getComm() );
  
}

NativeEpetraVector::~NativeEpetraVector() {}

Epetra_Vector &NativeEpetraVector::getEpetra_Vector() { return d_epetraVector; }

const Epetra_Vector &NativeEpetraVector::getEpetra_Vector() const { return d_epetraVector; }

AMP_MPI NativeEpetraVector::getComm() const {
  
  return std::dynamic_pointer_cast<NativeEpetraVectorParameters>(d_Params)->getComm();
}

Vector::shared_ptr NativeEpetraVector::cloneVector( const Variable::shared_ptr name ) const
{
  auto params = std::dynamic_pointer_cast<NativeEpetraVectorParameters>( d_Params );
  auto buffer = std::make_shared<VectorDataCPU<double>>( params->beginDOF(),
							 params->getLocalSize(),
							 params->getGlobalSize());
  
  auto retVal = std::make_shared<NativeEpetraVector>(d_Params, buffer);
  retVal->setVariable( name );
  return retVal;
}

void NativeEpetraVector::swapVectors( Vector &other )
{
    double *my_pointer;
    double *oth_pointer;
    getEpetra_Vector().ExtractView( &my_pointer );
    getEpetraVector( other ).ExtractView( &oth_pointer );
    getEpetraVector( other ).ResetView( my_pointer );
    getEpetra_Vector().ResetView( oth_pointer );
}

void NativeEpetraVector::aliasVector( Vector & )
{
  AMP_ERROR("NativeEpetraVector::aliasVector not implemented");
}

void NativeEpetraVector::assemble()
{
}

} // namespace LinearAlgebra
} // namespace AMP
