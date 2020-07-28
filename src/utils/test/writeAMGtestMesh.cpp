
#include <cstdio>
#include <cstdlib>
#include <iostream>

#define __NODE__( xi, yi, zi, nx, ny ) \
    ( ( xi ) + ( ( yi ) * ( nx ) ) + ( ( zi ) * ( nx ) * ( ny ) ) )

int main( int argc, char **argv )
{

    if ( argc < 8 ) {
        std::cout << "Usage: argv[0] nx ny nz Lx Ly Lz output_file" << std::endl;
        return 1;
    }

    int nx = atoi( argv[1] );
    int ny = atoi( argv[2] );
    int nz = atoi( argv[3] );

    double Lx = atof( argv[4] );
    double Ly = atof( argv[5] );
    double Lz = atof( argv[6] );

    double hx = Lx / ( static_cast<double>( nx - 1 ) );
    double hy = Ly / ( static_cast<double>( ny - 1 ) );
    double hz = Lz / ( static_cast<double>( nz - 1 ) );

    FILE *fp;

    fp = fopen( argv[7], "w" );

    fprintf( fp, "Mesh { \n" );

    int numPts = nx * ny * nz;
    fprintf( fp, "NumberOfNodes = %d \n", numPts );

    for ( int zi = 0, pi = 0; zi < nz; zi++ ) {
        for ( int yi = 0; yi < ny; yi++ ) {
            for ( int xi = 0; xi < nx; xi++, pi++ ) {
                double p[3];
                p[0] = ( static_cast<double>( xi ) ) * hx;
                p[1] = ( static_cast<double>( yi ) ) * hy;
                p[2] = ( static_cast<double>( zi ) ) * hz;
                fprintf( fp, "Point%d = %lf, %lf, %lf \n", pi, p[0], p[1], p[2] );
            } // end for xi
        }     // end for yi
    }         // end for zi

    fprintf( fp, "\n" );

    int numElem = ( nx - 1 ) * ( ny - 1 ) * ( nz - 1 );
    fprintf( fp, "NumberOfElements = %d \n", numElem );

    for ( int zi = 0, pi = 0; zi < ( nz - 1 ); zi++ ) {
        for ( int yi = 0; yi < ( ny - 1 ); yi++ ) {
            for ( int xi = 0; xi < ( nx - 1 ); xi++, pi++ ) {
                int p[8];
                p[0] = __NODE__( xi, yi, zi, nx, ny );
                p[1] = __NODE__( ( xi + 1 ), yi, zi, nx, ny );
                p[2] = __NODE__( ( xi + 1 ), ( yi + 1 ), zi, nx, ny );
                p[3] = __NODE__( xi, ( yi + 1 ), zi, nx, ny );
                p[4] = __NODE__( xi, yi, ( zi + 1 ), nx, ny );
                p[5] = __NODE__( ( xi + 1 ), yi, ( zi + 1 ), nx, ny );
                p[6] = __NODE__( ( xi + 1 ), ( yi + 1 ), ( zi + 1 ), nx, ny );
                p[7] = __NODE__( xi, ( yi + 1 ), ( zi + 1 ), nx, ny );
                fprintf( fp,
                         "Elem%d = %d, %d, %d, %d, %d, %d, %d, %d \n",
                         pi,
                         p[0],
                         p[1],
                         p[2],
                         p[3],
                         p[4],
                         p[5],
                         p[6],
                         p[7] );
            } // end for xi
        }     // end for yi
    }         // end for zi

    fprintf( fp, "\n" );

    {
        int num = 6;
        fprintf( fp, "NumberOfBoundaryNodeIds = %d \n\n", num );
    }

    // x = 0 surface
    {
        int num = ( ny * nz );
        fprintf( fp, "NumberOfBoundaryNodes1 = %d \n", num );
        fprintf( fp, "BoundaryNodeId1 = " );
    }

    for ( int zi = 0, cnt = 0; zi < nz; zi++ ) {
        for ( int yi = 0; yi < ny; yi++, cnt++ ) {
            int num = __NODE__( 0, yi, zi, nx, ny );
            fprintf( fp, "%d", num );
            if ( cnt < ( ( ny * nz ) - 1 ) ) {
                fprintf( fp, ", " );
            }
        } // end for yi
    }     // end for zi

    fprintf( fp, "\n\n" );

    // x = Lx surface
    {
        int num = ( ny * nz );
        fprintf( fp, "NumberOfBoundaryNodes2 = %d \n", num );
        fprintf( fp, "BoundaryNodeId2 = " );
    }

    for ( int zi = 0, cnt = 0; zi < nz; zi++ ) {
        for ( int yi = 0; yi < ny; yi++, cnt++ ) {
            int num = __NODE__( ( nx - 1 ), yi, zi, nx, ny );
            fprintf( fp, "%d", num );
            if ( cnt < ( ( ny * nz ) - 1 ) ) {
                fprintf( fp, ", " );
            }
        } // end for yi
    }     // end for zi

    fprintf( fp, "\n\n" );

    // y = 0 surface
    {
        int num = ( nx * nz );
        fprintf( fp, "NumberOfBoundaryNodes3 = %d \n", num );
        fprintf( fp, "BoundaryNodeId3 = " );
    }

    for ( int zi = 0, cnt = 0; zi < nz; zi++ ) {
        for ( int xi = 0; xi < nx; xi++, cnt++ ) {
            int num = __NODE__( xi, 0, zi, nx, ny );
            fprintf( fp, "%d", num );
            if ( cnt < ( ( nx * nz ) - 1 ) ) {
                fprintf( fp, ", " );
            }
        } // end for yi
    }     // end for zi

    fprintf( fp, "\n\n" );

    // y = Ly surface
    {
        int num = ( nx * nz );
        fprintf( fp, "NumberOfBoundaryNodes4 = %d \n", num );
        fprintf( fp, "BoundaryNodeId4 = " );
    }

    for ( int zi = 0, cnt = 0; zi < nz; zi++ ) {
        for ( int xi = 0; xi < nx; xi++, cnt++ ) {
            int num = __NODE__( xi, ( ny - 1 ), zi, nx, ny );
            fprintf( fp, "%d", num );
            if ( cnt < ( ( nx * nz ) - 1 ) ) {
                fprintf( fp, ", " );
            }
        } // end for yi
    }     // end for zi

    fprintf( fp, "\n\n" );

    // z = 0 surface
    {
        int num = ( ny * nx );
        fprintf( fp, "NumberOfBoundaryNodes5 = %d \n", num );
        fprintf( fp, "BoundaryNodeId5 = " );
    }

    for ( int yi = 0, cnt = 0; yi < ny; yi++ ) {
        for ( int xi = 0; xi < nx; xi++, cnt++ ) {
            int num = __NODE__( xi, yi, 0, nx, ny );
            fprintf( fp, "%d", num );
            if ( cnt < ( ( ny * nx ) - 1 ) ) {
                fprintf( fp, ", " );
            }
        } // end for xi
    }     // end for yi

    fprintf( fp, "\n\n" );

    // z = Lz surface
    {
        int num = ( ny * nx );
        fprintf( fp, "NumberOfBoundaryNodes6 = %d \n", num );
        fprintf( fp, "BoundaryNodeId6 = " );
    }

    for ( int yi = 0, cnt = 0; yi < ny; yi++ ) {
        for ( int xi = 0; xi < nx; xi++, cnt++ ) {
            int num = __NODE__( xi, yi, ( nz - 1 ), nx, ny );
            fprintf( fp, "%d", num );
            if ( cnt < ( ( ny * nx ) - 1 ) ) {
                fprintf( fp, ", " );
            }
        } // end for xi
    }     // end for yi

    fprintf( fp, "\n\n" );

    {
        int num = 0;
        fprintf( fp, "NumberOfBoundarySideIds = %d \n\n", num );
    }

    fprintf( fp, "} \n" );

    fclose( fp );
}
