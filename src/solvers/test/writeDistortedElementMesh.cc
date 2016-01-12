
#include <cstdio>
#include <cstdlib>
#include <iostream>

int main( int argc, char **argv )
{

    if ( argc == 1 ) {
        std::cout << "Usage: argv[0] output_file " << std::endl;
        return 1;
    }

    FILE *fp = fopen( argv[1], "w" );

    fprintf( fp, "\n Mesh { \n" );

    fprintf( fp, "NumberOfNodes = 8 \n" );

    fprintf( fp, "Point0 = 0.0, 0.0, 0.0 \n" );
    fprintf( fp, "Point1 = 1.0, 0.0, 0.0 \n" );
    fprintf( fp, "Point2 = 0.84, 1.567, 0.0 \n" );
    fprintf( fp, "Point3 = 0.0, 1.23, 0.0 \n" );
    fprintf( fp, "Point4 = 0.0, 0.1, 1.0 \n" );
    fprintf( fp, "Point5 = 1.0, 0.0, 1.0 \n" );
    fprintf( fp, "Point6 = 1.0, 1.0, 1.0 \n" );
    fprintf( fp, "Point7 = 0.0, 1.0, 1.0 \n" );

    fprintf( fp, "\n" );

    fprintf( fp, "NumberOfElements = 1 \n" );

    fprintf( fp, "Elem0 = 0, 1, 2, 3, 4, 5, 6, 7 \n" );

    fprintf( fp, "\n" );

    fprintf( fp, "NumberOfBoundaryNodeIds = 0 \n\n" );

    fprintf( fp, "NumberOfBoundarySideIds = 1 \n\n" );

    fprintf( fp, "NumberOfBoundarySides1 = 6 \n\n" );

    fprintf( fp, "BoundarySideId1 = 0, 0, 0, 1, 0, 2, 0, 3, 0, 4, 0, 5 \n\n" );

    fprintf( fp, "} \n" );

    fclose( fp );
}
