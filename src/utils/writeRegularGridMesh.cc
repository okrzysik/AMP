
#include <cstdio>
#include <cstdlib>
#include <iostream>

#define __NODE__( xi, yi, zi, nx, ny ) \
    ( ( xi ) + ( ( yi ) * ( nx ) ) + ( ( zi ) * ( nx ) * ( ny ) ) )

int main( int argc, char **argv )
{

    if ( argc < 9 ) {
        std::cout << "Usage: argv[0] nx ny nz Lx Ly Lz use_binary_format output_file" << std::endl;
        exit( 0 );
    }

    int nx = atoi( argv[1] );
    int ny = atoi( argv[2] );
    int nz = atoi( argv[3] );

    double Lx = atof( argv[4] );
    double Ly = atof( argv[5] );
    double Lz = atof( argv[6] );

    int useBinary = atoi( argv[7] );

    double hx = Lx / ( static_cast<double>( nx - 1 ) );
    double hy = Ly / ( static_cast<double>( ny - 1 ) );
    double hz = Lz / ( static_cast<double>( nz - 1 ) );

    FILE *fp;

    if ( useBinary ) {
        fp = fopen( argv[8], "wb" );
    }
    else {
        fp = fopen( argv[8], "w" );
    }

    if ( !useBinary ) {
        fprintf( fp, "Mesh { \n" );
    }

    int numPts = nx * ny * nz;
    if ( useBinary ) {
        fwrite( &numPts, sizeof( int ), 1, fp );
    }
    else {
        fprintf( fp, "NumberOfNodes = %d \n", numPts );
    }

    for ( int zi = 0, pi = 0; zi < nz; zi++ ) {
        for ( int yi = 0; yi < ny; yi++ ) {
            for ( int xi = 0; xi < nx; xi++, pi++ ) {
                double p[3];
                p[0] = ( static_cast<double>( xi ) ) * hx;
                p[1] = ( static_cast<double>( yi ) ) * hy;
                p[2] = ( static_cast<double>( zi ) ) * hz;
                if ( useBinary ) {
                    fwrite( p, sizeof( double ), 3, fp );
                }
                else {
                    fprintf( fp, "Point%d = %lf, %lf, %lf \n", pi, p[0], p[1], p[2] );
                }
            } // end for xi
        }     // end for yi
    }         // end for zi

    if ( !useBinary ) {
        fprintf( fp, "\n" );
    }

    int numElem = ( nx - 1 ) * ( ny - 1 ) * ( nz - 1 );
    if ( useBinary ) {
        fwrite( &numElem, sizeof( int ), 1, fp );
    }
    else {
        fprintf( fp, "NumberOfElements = %d \n", numElem );
    }

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
                if ( useBinary ) {
                    fwrite( p, sizeof( int ), 8, fp );
                }
                else {
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
                }
            } // end for xi
        }     // end for yi
    }         // end for zi

    if ( !useBinary ) {
        fprintf( fp, "\n" );
    }

    {
        int num = 6;
        if ( useBinary ) {
            fwrite( &num, sizeof( int ), 1, fp );
        }
        else {
            fprintf( fp, "NumberOfBoundaryNodeIds = %d \n\n", num );
        }
    }

    // x = 0 surface
    {
        int num = ( ny * nz );
        if ( useBinary ) {
            fwrite( &num, sizeof( int ), 1, fp );
        }
        else {
            fprintf( fp, "NumberOfBoundaryNodes1 = %d \n", num );
            fprintf( fp, "BoundaryNodeId1 = " );
        }
    }

    for ( int zi = 0, cnt = 0; zi < nz; zi++ ) {
        for ( int yi = 0; yi < ny; yi++, cnt++ ) {
            int num = __NODE__( 0, yi, zi, nx, ny );
            if ( useBinary ) {
                fwrite( &num, sizeof( int ), 1, fp );
            }
            else {
                fprintf( fp, "%d", num );
                if ( cnt < ( ( ny * nz ) - 1 ) ) {
                    fprintf( fp, ", " );
                }
            }
        } // end for yi
    }     // end for zi

    if ( !useBinary ) {
        fprintf( fp, "\n\n" );
    }

    // x = Lx surface
    {
        int num = ( ny * nz );
        if ( useBinary ) {
            fwrite( &num, sizeof( int ), 1, fp );
        }
        else {
            fprintf( fp, "NumberOfBoundaryNodes2 = %d \n", num );
            fprintf( fp, "BoundaryNodeId2 = " );
        }
    }

    for ( int zi = 0, cnt = 0; zi < nz; zi++ ) {
        for ( int yi = 0; yi < ny; yi++, cnt++ ) {
            int num = __NODE__( ( nx - 1 ), yi, zi, nx, ny );
            if ( useBinary ) {
                fwrite( &num, sizeof( int ), 1, fp );
            }
            else {
                fprintf( fp, "%d", num );
                if ( cnt < ( ( ny * nz ) - 1 ) ) {
                    fprintf( fp, ", " );
                }
            }
        } // end for yi
    }     // end for zi

    if ( !useBinary ) {
        fprintf( fp, "\n\n" );
    }

    // x = Lx surface Face Nodes
    {
        int num = ( ( ny - 2 ) * ( nz - 2 ) );
        if ( useBinary ) {
            fwrite( &num, sizeof( int ), 1, fp );
        }
        else {
            fprintf( fp, "NumberOfBoundaryNodes3 = %d \n", num );
        }
    }

    if ( ( ny > 2 ) && ( nz > 2 ) ) {
        if ( !useBinary ) {
            fprintf( fp, "BoundaryNodeId3 = " );
        }

        for ( int zi = 1, cnt = 0; zi < ( nz - 1 ); zi++ ) {
            for ( int yi = 1; yi < ( ny - 1 ); yi++, cnt++ ) {
                int num = __NODE__( ( nx - 1 ), yi, zi, nx, ny );
                if ( useBinary ) {
                    fwrite( &num, sizeof( int ), 1, fp );
                }
                else {
                    fprintf( fp, "%d", num );
                    if ( cnt < ( ( ( ny - 2 ) * ( nz - 2 ) ) - 1 ) ) {
                        fprintf( fp, ", " );
                    }
                }
            } // end for yi
        }     // end for zi
    }

    if ( !useBinary ) {
        fprintf( fp, "\n\n" );
    }

    // x = Lx surface Edge Nodes Type 1
    {
        int num = ( 2 * ( ny - 2 ) );
        if ( useBinary ) {
            fwrite( &num, sizeof( int ), 1, fp );
        }
        else {
            fprintf( fp, "NumberOfBoundaryNodes4 = %d \n", num );
        }
    }

    if ( ny > 2 ) {
        if ( !useBinary ) {
            fprintf( fp, "BoundaryNodeId4 = " );
        }

        for ( int yi = 1; yi < ( ny - 1 ); yi++ ) {
            int num = __NODE__( ( nx - 1 ), yi, 0, nx, ny );
            if ( useBinary ) {
                fwrite( &num, sizeof( int ), 1, fp );
            }
            else {
                fprintf( fp, "%d, ", num );
            }
        } // end for yi

        for ( int yi = 1; yi < ( ny - 2 ); yi++ ) {
            int num = __NODE__( ( nx - 1 ), yi, ( nz - 1 ), nx, ny );
            if ( useBinary ) {
                fwrite( &num, sizeof( int ), 1, fp );
            }
            else {
                fprintf( fp, "%d, ", num );
            }
        } // end for yi

        int num = __NODE__( ( nx - 1 ), ( ny - 2 ), ( nz - 1 ), nx, ny );
        if ( useBinary ) {
            fwrite( &num, sizeof( int ), 1, fp );
        }
        else {
            fprintf( fp, "%d ", num );
        }
    }

    if ( !useBinary ) {
        fprintf( fp, "\n\n" );
    }

    // x = Lx surface z Edge Nodes Type 2
    {
        int num = ( 2 * ( nz - 2 ) );
        if ( useBinary ) {
            fwrite( &num, sizeof( int ), 1, fp );
        }
        else {
            fprintf( fp, "NumberOfBoundaryNodes5 = %d \n", num );
        }
    }

    if ( nz > 2 ) {
        if ( !useBinary ) {
            fprintf( fp, "BoundaryNodeId5 = " );
        }

        for ( int zi = 1; zi < ( nz - 1 ); zi++ ) {
            int num = __NODE__( ( nx - 1 ), 0, zi, nx, ny );
            if ( useBinary ) {
                fwrite( &num, sizeof( int ), 1, fp );
            }
            else {
                fprintf( fp, "%d, ", num );
            }
        } // end for zi

        for ( int zi = 1; zi < ( nz - 2 ); zi++ ) {
            int num = __NODE__( ( nx - 1 ), ( ny - 1 ), zi, nx, ny );
            if ( useBinary ) {
                fwrite( &num, sizeof( int ), 1, fp );
            }
            else {
                fprintf( fp, "%d, ", num );
            }
        } // end for zi

        int num = __NODE__( ( nx - 1 ), ( ny - 1 ), ( nz - 2 ), nx, ny );
        if ( useBinary ) {
            fwrite( &num, sizeof( int ), 1, fp );
        }
        else {
            fprintf( fp, "%d ", num );
        }
    }

    if ( !useBinary ) {
        fprintf( fp, "\n\n" );
    }

    // x = Lx surface Corner Nodes
    {
        int num = 4;
        if ( useBinary ) {
            fwrite( &num, sizeof( int ), 1, fp );
        }
        else {
            fprintf( fp, "NumberOfBoundaryNodes6 = %d \n", num );
            fprintf( fp, "BoundaryNodeId6 = " );
        }
    }

    {
        int num[4];
        num[0] = __NODE__( ( nx - 1 ), 0, 0, nx, ny );
        num[1] = __NODE__( ( nx - 1 ), ( ny - 1 ), 0, nx, ny );
        num[2] = __NODE__( ( nx - 1 ), 0, ( nz - 1 ), nx, ny );
        num[3] = __NODE__( ( nx - 1 ), ( ny - 1 ), ( nz - 1 ), nx, ny );
        if ( useBinary ) {
            fwrite( num, sizeof( int ), 4, fp );
        }
        else {
            fprintf( fp, "%d, %d, %d, %d ", num[0], num[1], num[2], num[3] );
        }
    }

    if ( !useBinary ) {
        fprintf( fp, "\n\n" );
    }

    {
        int num = 0;
        if ( useBinary ) {
            fwrite( &num, sizeof( int ), 1, fp );
        }
        else {
            fprintf( fp, "NumberOfBoundarySideIds = %d \n\n", num );
        }
    }

    if ( !useBinary ) {
        fprintf( fp, "} \n" );
    }

    fclose( fp );
}
