#include "matrices/trilinos/EpetraMatrix.h"
#include "matrices/trilinos/ManagedEpetraMatrix.h"

#include <EpetraExt_Transpose_RowMatrix.h>

#ifdef USE_EXT_MPI
#include <Epetra_MpiComm.h>
#else
#include <Epetra_SerialComm.h>
#endif


namespace AMP {
namespace LinearAlgebra {


Matrix::shared_ptr EpetraMatrix::transpose() const
{
    EpetraExt::RowMatrix_Transpose transposer;
    return Matrix::shared_ptr( new ManagedEpetraMatrix(
        dynamic_cast<Epetra_CrsMatrix *>( &transposer( *d_epetraMatrix ) ), true ) );
}


EpetraMatrix::EpetraMatrix( Epetra_Map &map, Epetra_Map *, int *entities )
{
    // if ( colMap )
    //     d_epetraMatrix = new Epetra_FECrsMatrix ( Copy , map , *colMap , entities , false );
    // else
    d_epetraMatrix = new Epetra_FECrsMatrix( Copy, map, entities, false );
    d_DeleteMatrix = true;
}


Matrix::shared_ptr EpetraMatrix::createView( shared_ptr in_matrix )
{
    if ( in_matrix->isA<ManagedEpetraMatrix>() ) return in_matrix;

    AMP_ERROR( "Managed memory matrix is not well defined" );
    return Matrix::shared_ptr();
}


void EpetraMatrix::setEpetraMaps( Vector::shared_ptr range, Vector::shared_ptr domain )
{
    if ( range ) {
#ifdef USE_EXT_MPI
        Epetra_MpiComm comm = range->getComm().getCommunicator();
#else
        Epetra_SerialComm comm;
#endif
        AMP_INSIST( range->getGlobalSize() < 0x80000000,
                    "Epetra does not support vectors with global size greater than 2^31" );
        int N_global = static_cast<int>( range->getGlobalSize() );
        int N_local  = static_cast<int>( range->getLocalSize() );
        d_RangeMap   = AMP::shared_ptr<Epetra_Map>( new Epetra_Map( N_global, N_local, 0, comm ) );
        if ( domain ) {
            AMP_INSIST( domain->getGlobalSize() < 0x80000000,
                        "Epetra does not support vectors with global size greater than 2^31" );
            N_global = static_cast<int>( domain->getGlobalSize() );
            N_local  = static_cast<int>( domain->getLocalSize() );
            d_DomainMap =
                AMP::shared_ptr<Epetra_Map>( new Epetra_Map( N_global, N_local, 0, comm ) );
        }
    }
}


void EpetraMatrix::fillComplete()
{
    if ( d_RangeMap ) {
        if ( d_DomainMap )
            d_epetraMatrix->FillComplete( *d_DomainMap, *d_RangeMap );
        else
            d_epetraMatrix->FillComplete( *d_RangeMap, *d_RangeMap );
    }
    else {
        d_epetraMatrix->FillComplete();
    }
}
}
} // end namespace
