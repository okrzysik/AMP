#include "AMP/mesh/testHelpers/meshWriters.h"

#include <iostream>
#include <string>
#include <vector>


void run( const std::string &exe, const std::vector<const char *> &args )
{
    auto CHECK_ARGS = [&args]( int N, const char *msg ) {
        if ( (int) args.size() != N ) {
            std::cerr << msg;
            exit( 1 );
        }
    };

    if ( exe == "write2elementMesh" ) {
        CHECK_ARGS( 7, "Usage: <exe> write2elmenetMesh a ny nz Lx Ly Lz output_file\n" );
        double a  = atof( args[0] );
        int ny    = atoi( args[1] );
        int nz    = atoi( args[2] );
        double Lx = atof( args[3] );
        double Ly = atof( args[4] );
        double Lz = atof( args[5] );
        AMP::Mesh::MeshWriters::write2elementMesh( a, ny, nz, Lx, Ly, Lz, args[6] );
    } else if ( exe == "write3pointConstrainedBox" ) {
        CHECK_ARGS( 7, "Usage: <exe> write3pointConstrainedBox nx ny nz Lx Ly Lz output_file\n" );
        int nx    = atoi( args[0] );
        int ny    = atoi( args[1] );
        int nz    = atoi( args[2] );
        double Lx = atof( args[3] );
        double Ly = atof( args[4] );
        double Lz = atof( args[5] );
        AMP::Mesh::MeshWriters::writeConstrainedMesh( nx, ny, nz, Lx, Ly, Lz, args[6] );
    } else if ( exe == "write7elementMesh-1" ) {
        CHECK_ARGS( 1, "Usage: <exe> write7elementMesh-1 output_file\n" );
        AMP::Mesh::MeshWriters::write7elementMesh( 2, args[0] );
    } else if ( exe == "write7elementMesh-2" ) {
        CHECK_ARGS( 1, "Usage: <exe> write7elementMesh-2 output_file\n" );
        AMP::Mesh::MeshWriters::write7elementMesh( 8, args[0] );
    } else if ( exe == "writeBox" ) {
        CHECK_ARGS( 7, "Usage: <exe> writeBox nx ny nz Lx Ly Lz output_file\n" );
        int nx    = atoi( args[0] );
        int ny    = atoi( args[1] );
        int nz    = atoi( args[2] );
        double Lx = atof( args[3] );
        double Ly = atof( args[4] );
        double Lz = atof( args[5] );
        AMP::Mesh::MeshWriters::writeBox( nx, ny, nz, Lx, Ly, Lz, args[6] );
    } else if ( exe == "writeDispValsForPatchTest-2" ) {
        CHECK_ARGS( 1, "Usage: <exe> writeDispValsForPatchTest-2 output_file\n" );
        AMP::Mesh::MeshWriters::writeDispValsForPatchTest( args[0] );
    } else if ( exe == "writeDistortedElementMesh" ) {
        CHECK_ARGS( 1, "Usage: <exe> writeDistortedElementMesh output_file\n" );
        AMP::Mesh::MeshWriters::writeDistortedElement( args[0] );
    } else if ( exe == "writePlateWithHole" ) {
        CHECK_ARGS( 9, "Usage: <exe> writePlateWithHole le me ne pe a b c r output_file\n" );
        int le   = atoi( args[0] );
        int me   = atoi( args[1] );
        int ne   = atoi( args[2] );
        int pe   = atoi( args[3] );
        double a = atof( args[4] ); // X-dimension
        double b = atof( args[5] ); // Z-dimension
        double c = atof( args[6] ); // Y-dimension
        double r = atof( args[7] );
        AMP::Mesh::MeshWriters::writePlateWithHole( le, me, ne, pe, a, b, c, r, args[8] );
    } else if ( exe == "writeCookMesh" ) {
        CHECK_ARGS( 4, "Usage: <exe> writeCookMesh nx ny nz filename\n" );
        int nx = atoi( args[0] );
        int ny = atoi( args[1] );
        int nz = atoi( args[2] );
        AMP::Mesh::MeshWriters::writeCookMesh( nx, ny, nz, args[3] );
    } else if ( exe == "writeAMGMesh" ) {
        CHECK_ARGS( 7, "Usage: <exe> writeAMGMesh nx ny nz Lx Ly Lz filename\n" );
        int nx    = atoi( args[0] );
        int ny    = atoi( args[1] );
        int nz    = atoi( args[2] );
        double Lx = atof( args[3] );
        double Ly = atof( args[4] );
        double Lz = atof( args[5] );
        AMP::Mesh::MeshWriters::writeAMGMesh( nx, ny, nz, Lx, Ly, Lz, args[6] );
    } else {
        std::cerr << "Unknown mesh: " << exe << std::endl;
        exit( 1 );
    }
}


int main( int argc, char **argv )
{
    if ( argc == 1 ) {
        // clang-format off
        run( "write2elementMesh", { "0.5", "5", "10", "1.0", "2.0", "3.0", "outWrite2elementMesh" } );
        run( "write3pointConstrainedBox", { "4", "7", "10", "1.0", "2.0", "3.0", "outWrite3pointConstrainedBox" } );
        run( "writeBox", { "4", "7", "10", "1.0", "2.0", "3.0", "outWriteBox" } );
        run( "writeDistortedElementMesh", { "outDistortedElementMesh" } );
        run( "writeDispValsForPatchTest-2", { "outDispValsForPatchTest-2" } );
        // Meshes stored in AMP-Data and used by some tests
        run( "writeCookMesh", { "9", "2", "9", "cookMesh0" } );
        run( "writeCookMesh", { "17", "3", "17", "cookMesh1" } );
        run( "writeCookMesh", { "33", "5", "33", "cookMesh2" } );
        run( "writeCookMesh", { "65", "9", "65", "cookMesh3" } );
        run( "writeCookMesh", { "129", "17", "129", "cookMesh4" } );
        run( "writePlateWithHole", { "5", "15", "15", "15", "5.0", "10.0", "1.0", "1.0", "regPlateWithHole1" } );
        run( "writePlateWithHole", { "10", "30", "30", "30", "5.0", "10.0", "1.0", "1.0", "regPlateWithHole2" } );
        run( "write7elementMesh-1", { "mesh7elem-1" } );
        run( "write7elementMesh-2", { "mesh7elem-2" } );
        // clang-format on
    } else {
        std::string exe( argv[1] );
        std::vector<const char *> args2( &argv[2], &argv[2] + ( argc - 2 ) );
        run( exe, args2 );
    }
    return 0;
}
