#include "AMP/matrices/trilinos/EpetraMatrix.h"
#include "AMP/matrices/trilinos/ManagedEpetraMatrix.h"

DISABLE_WARNINGS
#include "Epetra_CrsMatrix.h"
#include "Epetra_FECrsMatrix.h"
#include <EpetraExt_Transpose_RowMatrix.h>
ENABLE_WARNINGS

#ifdef USE_EXT_MPI
#include <Epetra_MpiComm.h>
#else
#include <Epetra_SerialComm.h>
#endif


namespace AMP {
namespace LinearAlgebra {


void EpetraMatrix::VerifyEpetraReturn( int err, const char *func ) const
{
    std::stringstream error;
    error << func << ": " << err;
    if ( err < 0 )
        AMP_ERROR( error.str() );
    if ( err > 0 )
        AMP_ERROR( error.str() );
}

EpetraMatrix::EpetraMatrix( Epetra_CrsMatrix *inMatrix, bool dele )
    : d_epetraMatrix( inMatrix ), d_DeleteMatrix( dele )
{
}

EpetraMatrix::~EpetraMatrix()
{
    if ( d_DeleteMatrix )
        delete d_epetraMatrix;
}

Epetra_CrsMatrix &EpetraMatrix::getEpetra_CrsMatrix() { return *d_epetraMatrix; }

const Epetra_CrsMatrix &EpetraMatrix::getEpetra_CrsMatrix() const { return *d_epetraMatrix; }

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


std::shared_ptr<EpetraMatrix> EpetraMatrix::createView( shared_ptr in_matrix )
{
    auto mat = std::dynamic_pointer_cast<EpetraMatrix>( in_matrix );
    if ( !mat )
        AMP_ERROR( "Managed memory matrix is not well defined" );
    return mat;
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
        auto N_global = static_cast<int>( range->getGlobalSize() );
        auto N_local  = static_cast<int>( range->getLocalSize() );
        d_RangeMap    = std::make_shared<Epetra_Map>( N_global, N_local, 0, comm );
        if ( domain ) {
            AMP_INSIST( domain->getGlobalSize() < 0x80000000,
                        "Epetra does not support vectors with global size greater than 2^31" );
            N_global    = static_cast<int>( domain->getGlobalSize() );
            N_local     = static_cast<int>( domain->getLocalSize() );
            d_DomainMap = std::make_shared<Epetra_Map>( N_global, N_local, 0, comm );
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
    } else {
        d_epetraMatrix->FillComplete();
    }
}
} // namespace LinearAlgebra
} // namespace AMP
