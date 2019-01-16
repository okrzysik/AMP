
#include <cstdio>
#include <cstdlib>
#include <iostream>

#define __NODE__( xi, yi, zi, nx, ny ) \
    ( ( xi ) + ( ( yi ) * ( nx ) ) + ( ( zi ) * ( nx ) * ( ny ) ) )

int writeCookMesh( int argc, char **argv )
{

    if ( argc == 1 ) {
        std::cout << "Usage: " << argv[0] << " nx ny nz output_file " << std::endl;
        return 1;
    }

    int nx = atoi( argv[1] );
    int ny = atoi( argv[2] );
    int nz = atoi( argv[3] );

    double hx = 48.0 / static_cast<double>( nx - 1 );
    double hy = 1.0 / static_cast<double>( ny - 1 );

    FILE *fp = fopen( argv[4], "w" );

    fprintf( fp, "Mesh { \n" );

    int numPts = nx * ny * nz;
    fprintf( fp, "NumberOfNodes = %d \n", numPts );

    for ( int zi = 0, pi = 0; zi < nz; zi++ ) {
        for ( int yi = 0; yi < ny; yi++ ) {
            for ( int xi = 0; xi < nx; xi++, pi++ ) {
                double px = ( static_cast<double>( xi ) ) * hx;
                double py = ( static_cast<double>( yi ) ) * hy;
                double zl = 44.0 * px / 48.0;
                double zu = 44.0 + ( 16.0 * px / 48.0 );
                double pz = zl + ( static_cast<double>( zi ) * ( zu - zl ) /
                                   static_cast<double>( nz - 1 ) );
                fprintf( fp, "Point%d = %lf, %lf, %lf \n", pi, px, py, pz );
            } // end for xi
        }     // end for yi
    }         // end for zi

    fprintf( fp, "\n" );

    int numElem = ( nx - 1 ) * ( ny - 1 ) * ( nz - 1 );
    fprintf( fp, "NumberOfElements = %d \n", numElem );

    for ( int zi = 0, pi = 0; zi < ( nz - 1 ); zi++ ) {
        for ( int yi = 0; yi < ( ny - 1 ); yi++ ) {
            for ( int xi = 0; xi < ( nx - 1 ); xi++, pi++ ) {
                int p0 = __NODE__( xi, yi, zi, nx, ny );
                int p1 = __NODE__( ( xi + 1 ), yi, zi, nx, ny );
                int p2 = __NODE__( ( xi + 1 ), ( yi + 1 ), zi, nx, ny );
                int p3 = __NODE__( xi, ( yi + 1 ), zi, nx, ny );
                int p4 = __NODE__( xi, yi, ( zi + 1 ), nx, ny );
                int p5 = __NODE__( ( xi + 1 ), yi, ( zi + 1 ), nx, ny );
                int p6 = __NODE__( ( xi + 1 ), ( yi + 1 ), ( zi + 1 ), nx, ny );
                int p7 = __NODE__( xi, ( yi + 1 ), ( zi + 1 ), nx, ny );
                fprintf( fp,
                         "Elem%d = %d, %d, %d, %d, %d, %d, %d, %d \n",
                         pi,
                         p0,
                         p1,
                         p2,
                         p3,
                         p4,
                         p5,
                         p6,
                         p7 );
            } // end for xi
        }     // end for yi
    }         // end for zi

    fprintf( fp, "\n" );

    fprintf( fp, "NumberOfBoundaryNodeIds = 2 \n\n" );

    // x = 0 surface
    fprintf( fp, "NumberOfBoundaryNodes1 = %d \n", ( ny * nz ) );
    fprintf( fp, "BoundaryNodeId1 = " );
    for ( int zi = 0, cnt = 0; zi < nz; zi++ ) {
        for ( int yi = 0; yi < ny; yi++, cnt++ ) {
            fprintf( fp, "%d", __NODE__( 0, yi, zi, nx, ny ) );
            if ( cnt < ( ( ny * nz ) - 1 ) ) {
                fprintf( fp, ", " );
            }
        } // end for yi
    }     // end for zi
    fprintf( fp, "\n\n" );

    // x = Lx surface
    fprintf( fp, "NumberOfBoundaryNodes2 = %d \n", ( ny * nz ) );
    fprintf( fp, "BoundaryNodeId2 = " );
    for ( int zi = 0, cnt = 0; zi < nz; zi++ ) {
        for ( int yi = 0; yi < ny; yi++, cnt++ ) {
            fprintf( fp, "%d", __NODE__( ( nx - 1 ), yi, zi, nx, ny ) );
            if ( cnt < ( ( ny * nz ) - 1 ) ) {
                fprintf( fp, ", " );
            }
        } // end for yi
    }     // end for zi
    fprintf( fp, "\n\n" );

    fprintf( fp, "NumberOfBoundarySideIds = 0 \n\n" );

    fprintf( fp, "} \n" );

    fclose( fp );
}
