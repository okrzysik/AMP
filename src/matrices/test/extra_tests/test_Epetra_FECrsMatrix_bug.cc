/*
 *  This test tries to replicate a bug with the Epetra_FECrsMatrix using makeConsistent
 *  when the communication pattern.  Note: this test will fail (crash) if the
 *  -D_GLIBCXX_DEBUG -D_GLIBCXX_DEBUG_PEDANTIC flags are set.
 *
 */

#include "mpi.h"
#include <stdio.h>

#include <Epetra_FECrsMatrix.h>
#include <Epetra_Map.h>
#include <Epetra_MpiComm.h>


int main( int argc, char *argv[] )
{
    int rank, size;
    MPI_Init( &argc, &argv );
    MPI_Comm_size( MPI_COMM_WORLD, &size );
    MPI_Comm_rank( MPI_COMM_WORLD, &rank );
    if ( size != 2 ) {
        std::cout << "Test is only designed for 2 processors\n";
        return -1;
    }

    // Create the matrix
    Epetra_MpiComm comm = MPI_COMM_WORLD;
    Epetra_Map map( 4, 2, 0, comm );
    int entities[4] = { 4, 4, 4, 4 };
    Epetra_FECrsMatrix matrix( Copy, map, entities, false );

    // Zero out the matrix (using processor 0)
    if ( rank == 0 ) {
        int cols[4]      = { 1, 1, 1, 1 };
        double values[4] = { 0, 0, 0, 0 };
        for ( int i = 0; i < 4; i++ ) matrix.ReplaceGlobalValues( i, 4, values, cols );
    }
    matrix.GlobalAssemble( true );

    // Modify the diagonal (by the local owner)
    double one = 1.0;
    if ( rank == 0 ) {
        int rows[2] = { 0, 1 };
        matrix.ReplaceGlobalValues( rows[0], 1, &one, &rows[0] );
        matrix.ReplaceGlobalValues( rows[1], 1, &one, &rows[1] );
    }
    else {
        int rows[2] = { 2, 3 };
        matrix.ReplaceGlobalValues( rows[0], 1, &one, &rows[0] );
        matrix.ReplaceGlobalValues( rows[1], 1, &one, &rows[1] );
    }
    matrix.GlobalAssemble( true );

    MPI_Barrier( MPI_COMM_WORLD );
    MPI_Finalize();
    return 0;
}
