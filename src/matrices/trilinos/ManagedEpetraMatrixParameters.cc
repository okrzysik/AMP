#include "matrices/trilinos/ManagedEpetraMatrixParameters.h"
#include "utils/AMP_MPI.h"
#include "utils/Utilities.h"

#include "EpetraExt_MatrixMatrix.h"
#include "Epetra_Map.h"
#ifdef USE_EXT_MPI
#include <Epetra_MpiComm.h>
#else
#include <Epetra_SerialComm.h>
#endif

namespace AMP {
namespace LinearAlgebra {


template <class T>
static inline T *getPtr( std::vector<T> &x )
{
    if ( !x.empty() ) return &x[0];
    return NULL;
}


ManagedEpetraMatrixParameters::ManagedEpetraMatrixParameters(
    AMP::Discretization::DOFManager::shared_ptr left,
    AMP::Discretization::DOFManager::shared_ptr right,
    AMP_MPI comm )
    : MatrixParameters( left, right, comm )
{
    d_vEntriesPerRow.resize( getLocalNumberOfRows() );
}


Epetra_Map &ManagedEpetraMatrixParameters::getEpetraRowMap()
{
#ifdef USE_EXT_MPI
    Epetra_MpiComm comm = d_comm.getCommunicator();
#else
    Epetra_SerialComm comm;
#endif
    AMP_ASSERT( d_DOFManagerLeft.get() != NULL );
    AMP_ASSERT( d_DOFManagerRight.get() != NULL );
    AMP_INSIST( d_DOFManagerLeft->numGlobalDOF() < 0x80000000,
                "Epetra does not support vectors with global size greater than 2^31" );
    int N_row_local  = static_cast<int>( d_DOFManagerLeft->numLocalDOF() );
    int N_row_global = static_cast<int>( d_DOFManagerLeft->numGlobalDOF() );
    if ( d_eRowMap.get() == 0 ) {
        d_eRowMap =
            AMP::shared_ptr<Epetra_Map>( new Epetra_Map( N_row_global, N_row_local, 0, comm ) );
    }
    return *d_eRowMap;
}


Epetra_Map *ManagedEpetraMatrixParameters::getEpetraColMap()
{
#ifdef USE_EXT_MPI
    Epetra_MpiComm comm = d_comm.getCommunicator();
#else
    Epetra_SerialComm comm;
#endif
    return NULL;
    /*AMP_ASSERT(d_DOFManagerLeft.get()!=NULL);
    AMP_ASSERT(d_DOFManagerRight.get()!=NULL);
    AMP_INSIST(d_DOFManagerLeft->numGlobalDOF()<0x80000000,"Epetra does not support vectors with
    global size greater
    than 2^31");
    int N_col_local = static_cast<int>( d_DOFManagerRight->numLocalDOF() );
    int N_col_global = static_cast<int>( d_DOFManagerRight->numGlobalDOF() );
    AMP_ASSERT(d_ColGlobal==N_col_global);
    if ( d_eColMap.get() == 0 ) {
        std::vector<int> cols;
        cols.reserve ( d_sColumns.size() );
        for ( std::set<int>::iterator curCol = d_sColumns.begin(); curCol != d_sColumns.end() ;
    curCol++ )
          cols.push_back ( *curCol );
        d_eColMap = AMP::shared_ptr<Epetra_Map>( new Epetra_Map( -1, cols.size(), getPtr(cols), 0,
    comm ) );
        //d_eColMap = AMP::shared_ptr<Epetra_Map>( new Epetra_Map ( N_col_global, N_col_local,
    d_ColBase, comm ) );
    }
    return d_eColMap.get();*/
}


void ManagedEpetraMatrixParameters::addColumns( int a, int *b )
{
    for ( int i = 0; i != a; i++ ) d_sColumns.insert( b[i] );
}


const int *ManagedEpetraMatrixParameters::entryList() const { return &*d_vEntriesPerRow.begin(); }

int *ManagedEpetraMatrixParameters::entryList() { return &( d_vEntriesPerRow[0] ); }

void ManagedEpetraMatrixParameters::setEntriesInRow( int row, int entries )
{
    d_vEntriesPerRow[row] = entries;
}

int &ManagedEpetraMatrixParameters::entriesInRow( int i ) { return d_vEntriesPerRow[i]; }

int ManagedEpetraMatrixParameters::entriesInRow( int i ) const { return d_vEntriesPerRow[i]; }

int ManagedEpetraMatrixParameters::maxEntitiesInRow() const
{
    return *std::max_element( d_vEntriesPerRow.begin(), d_vEntriesPerRow.end() );
}

bool ManagedEpetraMatrixParameters::isSquare()
{
    AMP_ASSERT( d_DOFManagerLeft.get() != NULL );
    AMP_ASSERT( d_DOFManagerRight.get() != NULL );
    return d_DOFManagerLeft->numLocalDOF() == d_DOFManagerRight->numLocalDOF() &&
           d_DOFManagerLeft->numGlobalDOF() == d_DOFManagerRight->numGlobalDOF();
}

AMP::shared_ptr<Epetra_Map> ManagedEpetraMatrixParameters::getEpetraRowMapPtr()
{
    return d_eRowMap;
}

AMP::shared_ptr<Epetra_Map> ManagedEpetraMatrixParameters::getEpetraColMapPtr()
{
    return d_eColMap;
}

AMP_MPI ManagedEpetraMatrixParameters::getEpetraComm() { return d_comm; }
}
}
