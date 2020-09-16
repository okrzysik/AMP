#include "AMP/vectors/trilinos/epetra/EpetraVectorData.h"
#include "AMP/vectors/trilinos/epetra/EpetraVectorEngine.h"

#ifdef USE_EXT_MPI
#include <Epetra_MpiComm.h>

#include <utility>

#else
#include <Epetra_SerialComm.h>
#endif


namespace AMP::LinearAlgebra {


/********************************************************
 * EpetraVectorEngineParameters constructors             *
 ********************************************************/
EpetraVectorEngineParameters::EpetraVectorEngineParameters(
    std::shared_ptr<CommunicationList> commList,
    std::shared_ptr<AMP::Discretization::DOFManager> dofManager )
{
    d_CommList   = commList;
    d_DOFManager = dofManager;
    AMP_INSIST( d_DOFManager, "Requires a non null DOFManager" );
    auto local_size = d_DOFManager->numLocalDOF();
    d_global        = d_DOFManager->numGlobalDOF();
    d_comm          = d_DOFManager->getComm();
    d_comm.sumScan( &local_size, &d_end, 1 );
    d_begin = d_end - local_size;
}

EpetraVectorEngineParameters::EpetraVectorEngineParameters( size_t local_size,
                                                            size_t global_size,
                                                            const AMP_MPI &comm )
    : d_begin{ 0 }, d_end{ 0 }, d_global{ global_size }, d_comm{ comm }
{
    d_comm.sumScan( &local_size, &d_end, 1 );
    d_begin = d_end - local_size;
}

EpetraVectorEngineParameters::EpetraVectorEngineParameters( size_t local_size,
                                                            size_t global_size,
                                                            std::shared_ptr<Epetra_Map> emap,
                                                            const AMP_MPI &ecomm )
    : d_begin{ 0 },
      d_end{ 0 },
      d_global{ global_size },
      d_comm{ ecomm },
      d_emap( std::move( emap ) )
{
    d_comm.sumScan( &local_size, &d_end, 1 );
    d_begin = d_end - local_size;
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


} // namespace AMP::LinearAlgebra
