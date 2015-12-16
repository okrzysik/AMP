
#include <cstdio>
#include <cstdlib>
#include <iostream>

int main( int argc, char **argv )
{

    if ( argc == 1 ) {
        std::cout << "Usage: argv[0] output_file " << std::endl;
        exit( 0 );
    }

    double xarr[8] = { 0.0, 1.0, 1.0, 0.0, 0.0, 1.0, 1.0, 0.0 };
    double yarr[8] = { 0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 1.0, 1.0 };
    double zarr[8] = { 0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0 };

    FILE *fp = fopen( argv[1], "w" );

    fprintf( fp, "\n\n" );

    fprintf( fp, "number_of_ids = 8 \n\n" );

    for ( int i = 0; i < 8; i++ ) {
        fprintf( fp, "id_%d = %d \n", i, ( i + 1 ) );
        fprintf( fp, "number_of_dofs_%d = 3 \n", i );
        fprintf( fp, "dof_%d_0 = 0 \n", i );
        fprintf( fp, "dof_%d_1 = 1 \n", i );
        fprintf( fp, "dof_%d_2 = 2 \n", i );
        fprintf( fp,
                 "value_%d_0 = %lf \n",
                 i,
                 ( 0.5 * 1.0e-3 * ( ( 2.0 * xarr[i] ) + yarr[i] + zarr[i] ) ) );
        fprintf( fp,
                 "value_%d_1 = %lf \n",
                 i,
                 ( 0.5 * 1.0e-3 * ( xarr[i] + ( 2.0 * yarr[i] ) + zarr[i] ) ) );
        fprintf( fp,
                 "value_%d_2 = %lf \n",
                 i,
                 ( 0.5 * 1.0e-3 * ( xarr[i] + yarr[i] + ( 2.0 * zarr[i] ) ) ) );
        fprintf( fp, "\n" );
    } // end for i

    fprintf( fp, "\n\n" );

    fclose( fp );
}
