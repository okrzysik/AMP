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

    if ( exe == "write2elmenetMesh" ) {
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
    } else if ( exe == "writePlateWithHoleMesh" ) {
        CHECK_ARGS( 9, "Usage: <exe> writePlateWithHoleMesh le me ne pe a b c r output_file\n" );
        int const le   = atoi( args[0] );
        int const me   = atoi( args[1] );
        int const ne   = atoi( args[2] );
        int const pe   = atoi( args[3] );
        double const a = atof( args[4] ); // X-dimension
        double const b = atof( args[5] ); // Z-dimension
        double const c = atof( args[6] ); // Y-dimension
        double const r = atof( args[7] );
        AMP::Mesh::MeshWriters::writePlateWithHole( le, me, ne, pe, a, b, c, r, args[8] );
    } else if ( exe == "writeCookMesh" ) {
        CHECK_ARGS( 4, "Usage: <exe> writeCookMesh nx ny nz filename\n" );
        int nx = atoi( args[0] );
        int ny = atoi( args[1] );
        int nz = atoi( args[2] );
        AMP::Mesh::MeshWriters::writeCookMesh( nx, ny, nz, args[3] );
    } else {
        std::cerr << "Unknown mesh: " << exe << std::endl;
        exit( 1 );
    }
}


int main( int argc, char **argv )
{
    if ( argc == 1 ) {
        run( "write7elementMesh-1", { "out7elementMesh-1" } );
        run( "write7elementMesh-2", { "out7elementMesh-2" } );
        run( "writeDistortedElementMesh", { "outDistortedElementMesh" } );
        run( "writeDispValsForPatchTest-2", { "outDispValsForPatchTest-2" } );
        run( "writeCookMesh", { "10", "10", "10", "test.dat" } );
        return 0;
    } else {
        std::string exe( argv[1] );
        std::vector<const char *> args2( &argv[2], &argv[2] + ( argc - 2 ) );
        run( exe, args2 );
    }
    return 0;
}
