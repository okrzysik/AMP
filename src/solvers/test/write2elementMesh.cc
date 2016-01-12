
#include <cstdio>
#include <cstdlib>
#include <iostream>

#define __NODE__( xi, yi, zi, ny ) ( ( xi ) + ( 3 * ( yi ) ) + ( 3 * ( zi ) * ( ny ) ) )

int main( int argc, char **argv )
{

    if ( argc == 1 ) {
        std::cout << "Usage: argv[0] a ny nz Lx Ly Lz output_file " << std::endl;
        return 1;
    }

    double a = atof( argv[1] );

    int ny = atoi( argv[2] );
    int nz = atoi( argv[3] );

    double Lx = atof( argv[4] );
    double Ly = atof( argv[5] );
    double Lz = atof( argv[6] );

    double hx = Lx / 2.0;
    double hy = Ly / ( static_cast<double>( ny - 1 ) );
    double hz = Lz / ( static_cast<double>( nz - 1 ) );

    FILE *fp = fopen( argv[7], "w" );

    fprintf( fp, "Mesh { \n" );

    int numPts = 3 * ny * nz;
    fprintf( fp, "NumberOfNodes = %d \n", numPts );

    for ( int zi = 0, pi = 0; zi < nz; zi++ ) {
        for ( int yi = 0; yi < ny; yi++ ) {
            {
                int xi    = 0;
                double px = ( static_cast<double>( xi ) ) * hx;
                double py = ( static_cast<double>( yi ) ) * hy;
                double pz = ( static_cast<double>( zi ) ) * hz;
                fprintf( fp, "Point%d = %lf, %lf, %lf \n", pi, px, py, pz );
                pi++;
            }

            {
                int xi      = 1;
                double pz   = ( static_cast<double>( zi ) ) * hz;
                double zTmp = ( 0.5 * Lz ) - pz;
                double px   = ( ( static_cast<double>( xi ) ) * hx ) + ( 2.0 * zTmp * a / Lz );
                double py   = ( static_cast<double>( yi ) ) * hy;
                fprintf( fp, "Point%d = %lf, %lf, %lf \n", pi, px, py, pz );
                pi++;
            }

            {
                int xi    = 2;
                double px = ( static_cast<double>( xi ) ) * hx;
                double py = ( static_cast<double>( yi ) ) * hy;
                double pz = ( static_cast<double>( zi ) ) * hz;
                fprintf( fp, "Point%d = %lf, %lf, %lf \n", pi, px, py, pz );
                pi++;
            }
        } // end for yi
    }     // end for zi

    fprintf( fp, "\n" );

    int numElem = 2 * ( ny - 1 ) * ( nz - 1 );
    fprintf( fp, "NumberOfElements = %d \n", numElem );

    for ( int zi = 0, pi = 0; zi < ( nz - 1 ); zi++ ) {
        for ( int yi = 0; yi < ( ny - 1 ); yi++ ) {
            for ( int xi = 0; xi < 2; xi++, pi++ ) {
                int p0 = __NODE__( xi, yi, zi, ny );
                int p1 = __NODE__( ( xi + 1 ), yi, zi, ny );
                int p2 = __NODE__( ( xi + 1 ), ( yi + 1 ), zi, ny );
                int p3 = __NODE__( xi, ( yi + 1 ), zi, ny );
                int p4 = __NODE__( xi, yi, ( zi + 1 ), ny );
                int p5 = __NODE__( ( xi + 1 ), yi, ( zi + 1 ), ny );
                int p6 = __NODE__( ( xi + 1 ), ( yi + 1 ), ( zi + 1 ), ny );
                int p7 = __NODE__( xi, ( yi + 1 ), ( zi + 1 ), ny );
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

    fprintf( fp, "NumberOfBoundaryNodeIds = 4 \n\n" );

    // x = 0, z = 0 edge
    fprintf( fp, "NumberOfBoundaryNodes1 = %d \n", ny );
    fprintf( fp, "BoundaryNodeId1 = " );
    for ( int yi = 0; yi < ny; yi++ ) {
        fprintf( fp, "%d", __NODE__( 0, yi, 0, ny ) );
        if ( yi < ( ny - 1 ) ) {
            fprintf( fp, ", " );
        }
    } // end for yi
    fprintf( fp, "\n\n" );

    // x = 0, z = (nz - 1) edge
    fprintf( fp, "NumberOfBoundaryNodes2 = %d \n", ny );
    fprintf( fp, "BoundaryNodeId2 = " );
    for ( int yi = 0; yi < ny; yi++ ) {
        fprintf( fp, "%d", __NODE__( 0, yi, ( nz - 1 ), ny ) );
        if ( yi < ( ny - 1 ) ) {
            fprintf( fp, ", " );
        }
    } // end for yi
    fprintf( fp, "\n\n" );

    // x = 2, z = 0 edge
    fprintf( fp, "NumberOfBoundaryNodes3 = %d \n", ny );
    fprintf( fp, "BoundaryNodeId3 = " );
    for ( int yi = 0; yi < ny; yi++ ) {
        fprintf( fp, "%d", __NODE__( 2, yi, 0, ny ) );
        if ( yi < ( ny - 1 ) ) {
            fprintf( fp, ", " );
        }
    } // end for yi
    fprintf( fp, "\n\n" );

    // x = 2, z = (nz - 1) edge
    fprintf( fp, "NumberOfBoundaryNodes4 = %d \n", ny );
    fprintf( fp, "BoundaryNodeId4 = " );
    for ( int yi = 0; yi < ny; yi++ ) {
        fprintf( fp, "%d", __NODE__( 2, yi, ( nz - 1 ), ny ) );
        if ( yi < ( ny - 1 ) ) {
            fprintf( fp, ", " );
        }
    } // end for yi
    fprintf( fp, "\n\n" );

    fprintf( fp, "NumberOfBoundarySideIds = 0 \n\n" );

    fprintf( fp, "} \n" );

    fprintf( fp, "MeshParams { \n" );
    fprintf( fp, "a = %lf \n", a );
    fprintf( fp, "ny = %d \n", ny );
    fprintf( fp, "nz = %d \n", nz );
    fprintf( fp, "Lx = %lf \n", Lx );
    fprintf( fp, "Ly = %lf \n", Ly );
    fprintf( fp, "Lz = %lf \n", Lz );
    fprintf( fp, "} \n\n" );

    fclose( fp );
}
