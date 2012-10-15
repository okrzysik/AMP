#include "matrices/trilinos/ManagedEpetraMatrixParameters.h"
#include "utils/Utilities.h"
#include "utils/AMP_MPI.h"

#include "EpetraExt_MatrixMatrix.h"
#ifdef USE_EXT_MPI
    #include <Epetra_MpiComm.h>
#else
    #include <Epetra_SerialComm.h>
#endif

namespace AMP {
namespace LinearAlgebra {




ManagedEpetraMatrixParameters::ManagedEpetraMatrixParameters ( int local_size , int global_size , int first_dof , AMP_MPI comm )
        : MatrixParameters () ,
          d_comm ( comm ) ,
          d_vEntriesPerRow ( local_size ) ,
          d_ColGlobal ( global_size ) ,
          d_ColBase ( first_dof ) ,
          d_RowBase ( first_dof )
{
}


ManagedEpetraMatrixParameters::ManagedEpetraMatrixParameters ( int local_size , int global_size , int row_first_dof , int col_global , int col_first_dof , AMP_MPI comm )
        : MatrixParameters () ,
          d_comm ( comm ) ,
          d_vEntriesPerRow ( local_size ) ,
          d_ColGlobal ( col_global ) ,
          d_ColBase ( col_first_dof ) ,
          d_RowBase ( row_first_dof )
{
}


Epetra_Map  &ManagedEpetraMatrixParameters::getEpetraRowMap ()
{
    #ifdef USE_EXT_MPI
        Epetra_MpiComm  comm = d_comm.getCommunicator();
    #else
        Epetra_SerialComm  comm;
    #endif
    AMP_ASSERT(d_DOFManagerLeft.get()!=NULL);
    AMP_ASSERT(d_DOFManagerRight.get()!=NULL);
    AMP_INSIST(d_DOFManagerLeft->numGlobalDOF()<0x80000000,"Epetra does not support vectors with global size greater than 2^31");
    int N_row_local = static_cast<int>( d_DOFManagerLeft->numLocalDOF() );
    int N_row_global = static_cast<int>( d_DOFManagerLeft->numGlobalDOF() );
    if ( d_eRowMap.get() == 0 ) {
        d_eRowMap = boost::shared_ptr<Epetra_Map>( new Epetra_Map ( N_row_global, N_row_local, d_RowBase, comm ) );
    }
    return *d_eRowMap;
}


Epetra_Map  &ManagedEpetraMatrixParameters::getEpetraColMap ()
{
    #ifdef USE_EXT_MPI
        Epetra_MpiComm  comm = d_comm.getCommunicator();
    #else
        Epetra_SerialComm  comm;
    #endif
    AMP_ASSERT(d_DOFManagerLeft.get()!=NULL);
    AMP_ASSERT(d_DOFManagerRight.get()!=NULL);
    AMP_INSIST(d_DOFManagerLeft->numGlobalDOF()<0x80000000,"Epetra does not support vectors with global size greater than 2^31");
    int N_col_local = static_cast<int>( d_DOFManagerRight->numLocalDOF() );
    int N_col_global = static_cast<int>( d_DOFManagerRight->numGlobalDOF() );
    if ( d_eColMap.get() == 0 ) {
        /*if ( d_ColBase < 0 ) return 0;
        std::vector<int> cols;
        cols.reserve ( d_sColumns.size() );
        for ( std::set<int>::iterator curCol = d_sColumns.begin(); curCol != d_sColumns.end() ; curCol++ )
          cols.push_back ( *curCol );
        if ( d_eColMap.get() == 0 )
          d_eColMap = boost::shared_ptr<Epetra_Map>( new Epetra_Map( d_ColGlobal, cols.size(), &cols[0], 0, comm ) );*/
        d_eColMap = boost::shared_ptr<Epetra_Map>( new Epetra_Map ( N_col_global, N_col_local, d_ColBase, comm ) );
    }
    return *d_eColMap;
}


void ManagedEpetraMatrixParameters::addColumns ( int a , int *b )
{
    for ( int i = 0 ; i != a ; i++ )
      d_sColumns.insert ( b[i] );
}


const int* ManagedEpetraMatrixParameters::entryList () const       
{ 
    return &*d_vEntriesPerRow.begin(); 
}

int* ManagedEpetraMatrixParameters::entryList()              
{ 
    return &(d_vEntriesPerRow[0]); 
}

void ManagedEpetraMatrixParameters::setEntriesInRow ( int row , int entries ) 
{ 
    d_vEntriesPerRow[row] = entries; 
}

int& ManagedEpetraMatrixParameters::entriesInRow ( int i )       
{ 
    return d_vEntriesPerRow[i]; 
}

int ManagedEpetraMatrixParameters::entriesInRow ( int i ) const 
{ 
    return d_vEntriesPerRow[i]; 
}

int    ManagedEpetraMatrixParameters::maxEntitiesInRow () const 
{ 
    return *std::max_element ( d_vEntriesPerRow.begin() , d_vEntriesPerRow.end() ); 
}

bool  ManagedEpetraMatrixParameters::isSquare () 
{ 
    return d_ColBase < 0; 
}

boost::shared_ptr < Epetra_Map >   ManagedEpetraMatrixParameters::getEpetraRowMapPtr () 
{ 
    return d_eRowMap; 
}

boost::shared_ptr < Epetra_Map >   ManagedEpetraMatrixParameters::getEpetraColMapPtr () 
{ 
    return d_eColMap; 
}

AMP_MPI   ManagedEpetraMatrixParameters::getEpetraComm () 
{ 
    return d_comm; 
}


}
}

