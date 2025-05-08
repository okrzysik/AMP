#include "AMP/mesh/testHelpers/meshWriters.h"
#include "AMP/utils/Database.h"
#include "AMP/utils/UtilityMacros.h"

#include "AMP/AMP_TPLs.h"
#ifdef AMP_USE_LIBMESH
    #include "AMP/mesh/libmesh/libmeshMesh.h"
    #include "AMP/utils/AMP_MPI.h"
DISABLE_WARNINGS
    #include "libmesh/libmesh_config.h"
    #undef LIBMESH_ENABLE_REFERENCE_COUNTING
    #include "libmesh/boundary_info.h"
    #include "libmesh/cell_hex8.h"
    #include "libmesh/elem.h"
    #include "libmesh/mesh.h"
    #include "libmesh/mesh_communication.h"
ENABLE_WARNINGS
#endif

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <sstream>
#include <vector>

#include "ProfilerApp.h"


static constexpr double PI = 3.1415926535897932;


namespace AMP::Mesh::MeshWriters {


static DatabasePtr convertDatabase( DatabasePtr db );

template<class TYPE>
static inline void fread( TYPE *data, size_t N, FILE *fp )
{
    size_t n = fread( data, sizeof( TYPE ), N, fp );
    AMP_INSIST( n == N, "Error while reading the file" );
}
template<class TYPE>
static inline void fwrite( const TYPE *data, size_t N, FILE *fp )
{
    size_t n = fwrite( data, sizeof( TYPE ), N, fp );
    AMP_INSIST( n == N, "Error while writing the file" );
}


/********************************************************
 * Put/Get data to the database                          *
 ********************************************************/
template<class TYPE, size_t SIZE>
static void
putVector( AMP::Database &db, std::string_view key, const std::vector<std::array<TYPE, SIZE>> &x )
{
    if ( x.empty() )
        db.putArray<TYPE>( key, {} );
    AMP::Array<TYPE> data2;
    data2.viewRaw( { SIZE, x.size() }, const_cast<TYPE *>( x[0].data() ) );
    db.putArray( key, data2 );
}
static void
putVector( AMP::Database &db, std::string_view key, const std::vector<std::vector<int>> x )
{
    db.putScalar<int>( key, x.size() );
    for ( size_t i = 0; i < x.size(); i++ ) {
        auto key2 = std::string( key ) + std::to_string( i + 1 );
        db.putVector( key2, x[i] );
    }
}
template<class TYPE, size_t SIZE>
static std::vector<std::array<TYPE, SIZE>> getVector( const AMP::Database &db,
                                                      std::string_view key )
{
    auto x = db.getArray<TYPE>( key );
    AMP_ASSERT( x.size( 0 ) == SIZE );
    std::vector<std::array<TYPE, SIZE>> y( x.size( 1 ) );
    memcpy( y.data(), x.data(), sizeof( TYPE ) * x.length() );
    return y;
}
std::vector<std::vector<int>> getVector( const AMP::Database &db, std::string_view key )
{
    int N = db.getScalar<int>( key );
    std::vector<std::vector<int>> x( N );
    for ( size_t i = 0; i < x.size(); i++ ) {
        auto key2 = std::string( key ) + std::to_string( i + 1 );
        x[i]      = db.getVector<int>( key2 );
    }
    return x;
}


/********************************************************
 * Generate the database for a known mesh                *
 ********************************************************/
std::shared_ptr<AMP::Database> generateTestMesh( const std::string &name )
{
    PROFILE( "generateTestMesh" );
    if ( name == "distortedElementMesh" ) {
        return createDistortedElement();
    } else if ( name == "cookMesh0" ) {
        return createCookMesh( 9, 2, 9 );
    } else if ( name == "cookMesh1" ) {
        return createCookMesh( 17, 3, 17 );
    } else if ( name == "cookMesh2" ) {
        return createCookMesh( 33, 5, 33 );
    } else if ( name == "cookMesh3" ) {
        return createCookMesh( 65, 9, 65 );
    } else if ( name == "cookMesh4" ) {
        return createCookMesh( 129, 17, 129 );
    } else if ( name == "regPlateWithHole1" ) {
        return createPlateWithHole( 5, 15, 15, 15, 5.0, 10.0, 1.0, 1.0 );
    } else if ( name == "regPlateWithHole2" ) {
        return createPlateWithHole( 10, 30, 30, 30, 5.0, 10.0, 1.0, 1.0 );
    } else if ( name == "mesh7elem-1" ) {
        return create7elementMesh( 2 );
    } else if ( name == "mesh7elem-2" ) {
        return create7elementMesh( 8 );
    } else if ( name == "boxMesh-1" ) {
        return createAMGMesh( 2, 2, 2, 10, 10, 10 );
    } else if ( name == "boxMesh-2" ) {
        return createAMGMesh( 3, 3, 3, 10, 10, 10 );
    } else if ( name == "boxMesh-3" ) {
        return createAMGMesh( 5, 5, 5, 10, 10, 10 );
    } else if ( name == "boxMesh-4" ) {
        return createAMGMesh( 9, 9, 9, 10, 10, 10 );
    } else if ( name == "boxMesh-5" ) {
        return createAMGMesh( 17, 17, 17, 10, 10, 10 );
    }
    return {};
}


/********************************************************
 * Create and write all known test meshes                *
 ********************************************************/
void generateAll()
{
    const char *ascii[] = { "distortedElementMesh",
                            "cookMesh0",
                            "cookMesh1",
                            "cookMesh2",
                            "cookMesh3",
                            "cookMesh4",
                            "regPlateWithHole1",
                            "regPlateWithHole2",
                            "mesh7elem-1",
                            "mesh7elem-2",
                            "boxMesh-1",
                            "boxMesh-2",
                            "boxMesh-3",
                            "boxMesh-4",
                            "boxMesh-5" };
    // const char* binary[] = {};
    for ( auto name : ascii ) {
        auto db = generateTestMesh( name );
        AMP_ASSERT( db );
        writeTestMesh( *db, name );
    }
    /*for ( auto name : binary ) {
        auto db = generateTestMesh( name );
        AMP_ASSERT( db );
        writeBinaryTestMesh( *db, name );
    }*/
}


/********************************************************
 * create a box mesh                                     *
 ********************************************************/
DatabasePtr createBox( int nx, int ny, int nz, double Lx, double Ly, double Lz )
{
    auto NODE = [nx, ny]( int xi, int yi, int zi ) { return xi + yi * nx + zi * nx * ny; };
    double hx = Lx / ( static_cast<double>( nx - 1 ) );
    double hy = Ly / ( static_cast<double>( ny - 1 ) );
    double hz = Lz / ( static_cast<double>( nz - 1 ) );

    int numNode = nx * ny * nz;
    std::vector<std::array<double, 3>> nodes( numNode );
    for ( int zi = 0, k = 0, pi = 0; zi < nz; zi++ ) {
        for ( int yi = 0; yi < ny; yi++ ) {
            for ( int xi = 0; xi < nx; xi++, pi++, k++ ) {
                nodes[k] = { xi * hx, yi * hy, zi * hz };
            }
        }
    }

    int numElem = ( nx - 1 ) * ( ny - 1 ) * ( nz - 1 );
    std::vector<std::array<double, 8>> elem( numElem );
    for ( int zi = 0, k = 0, pi = 0; zi < ( nz - 1 ); zi++ ) {
        for ( int yi = 0; yi < ( ny - 1 ); yi++ ) {
            for ( int xi = 0; xi < ( nx - 1 ); xi++, pi++, k++ ) {
                elem[k][0] = NODE( xi, yi, zi );
                elem[k][1] = NODE( ( xi + 1 ), yi, zi );
                elem[k][2] = NODE( ( xi + 1 ), ( yi + 1 ), zi );
                elem[k][3] = NODE( xi, ( yi + 1 ), zi );
                elem[k][4] = NODE( xi, yi, ( zi + 1 ) );
                elem[k][5] = NODE( ( xi + 1 ), yi, ( zi + 1 ) );
                elem[k][6] = NODE( ( xi + 1 ), ( yi + 1 ), ( zi + 1 ) );
                elem[k][7] = NODE( xi, ( yi + 1 ), ( zi + 1 ) );
            }
        }
    }

    // z = Lz surface
    std::vector<std::vector<int>> BoundaryNodes( 2 );
    BoundaryNodes[0] = std::vector<int>( nx * ny );
    auto ids         = BoundaryNodes[0].data();
    for ( int yi = 0, k = 0; yi < ny; yi++ ) {
        for ( int xi = 0; xi < nx; xi++, k++ ) {
            ids[k] = NODE( xi, yi, ( nz - 1 ) );
        }
    }

    // z = 0 surface
    BoundaryNodes[1] = std::vector<int>( nx * ny );
    ids              = BoundaryNodes[1].data();
    for ( int yi = 0, k = 0; yi < ny; yi++ ) {
        for ( int xi = 0; xi < nx; xi++, k++ ) {
            ids[k] = NODE( xi, yi, 0 );
        }
    }

    std::vector<std::vector<int>> SideIds;

    auto db   = std::make_unique<AMP::Database>();
    auto mesh = db->createAddDatabase( "Mesh" );
    putVector( *mesh, "Nodes", nodes );
    putVector( *mesh, "Elems", elem );
    putVector( *mesh, "BoundaryNodes", BoundaryNodes );
    putVector( *mesh, "SideIds", SideIds );
    db->putScalar( "Lx", Lx );
    db->putScalar( "Ly", Ly );
    db->putScalar( "Lz", Lz );
    db->putScalar( "nx", nx );
    db->putScalar( "ny", ny );
    db->putScalar( "nz", nz );
    return db;
}


/********************************************************
 * write displacements for patch test                    *
 ********************************************************/
void writeDispValsForPatchTest( const std::string &file )
{
    const double x[8] = { 0.0, 1.0, 1.0, 0.0, 0.0, 1.0, 1.0, 0.0 };
    const double y[8] = { 0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 1.0, 1.0 };
    const double z[8] = { 0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0 };

    FILE *fp = fopen( file.data(), "w" );
    fprintf( fp, "\n\n" );
    fprintf( fp, "number_of_ids = 8 \n\n" );
    for ( int i = 0; i < 8; i++ ) {
        fprintf( fp, "id_%d = %d \n", i, ( i + 1 ) );
        fprintf( fp, "number_of_dofs_%d = 3 \n", i );
        fprintf( fp, "dof_%d_0 = 0 \n", i );
        fprintf( fp, "dof_%d_1 = 1 \n", i );
        fprintf( fp, "dof_%d_2 = 2 \n", i );
        fprintf( fp, "value_%d_0 = %lf \n", i, ( 0.0005 * ( ( 2 * x[i] ) + y[i] + z[i] ) ) );
        fprintf( fp, "value_%d_1 = %lf \n", i, ( 0.0005 * ( x[i] + ( 2 * y[i] ) + z[i] ) ) );
        fprintf( fp, "value_%d_2 = %lf \n", i, ( 0.0005 * ( x[i] + y[i] + ( 2 * z[i] ) ) ) );
        fprintf( fp, "\n" );
    }
}


/********************************************************
 * write a distorted element mesh                        *
 ********************************************************/
DatabasePtr createDistortedElement()
{
    std::vector<std::array<double, 3>> nodes( 8 );
    nodes[0] = { 0.0, 0.0, 0.0 };
    nodes[1] = { 1.0, 0.0, 0.0 };
    nodes[2] = { 0.84, 1.567, 0.0 };
    nodes[3] = { 0.0, 1.23, 0.0 };
    nodes[4] = { 0.0, 0.1, 1.0 };
    nodes[5] = { 1.0, 0.0, 1.0 };
    nodes[6] = { 1.0, 1.0, 1.0 };
    nodes[7] = { 0.0, 1.0, 1.0 };
    std::vector<std::array<double, 8>> elem( 1 );
    elem[0] = { 0, 1, 2, 3, 4, 5, 6, 7 };
    std::vector<std::vector<int>> BoundaryNodes;
    std::vector<std::vector<int>> SideIds = { { 0, 0, 0, 1, 0, 2, 0, 3, 0, 4, 0, 5 } };
    auto db                               = std::make_unique<AMP::Database>();
    auto mesh                             = db->createAddDatabase( "Mesh" );
    putVector( *mesh, "Nodes", nodes );
    putVector( *mesh, "Elems", elem );
    putVector( *mesh, "BoundaryNodes", BoundaryNodes );
    putVector( *mesh, "SideIds", SideIds );
    return db;
}

/********************************************************
 * write a plate with a hole                             *
 ********************************************************/
DatabasePtr
createPlateWithHole( int le, int me, int ne, int pe, double a, double b, double c, double r )
{
    int const ze[] = { ne, pe, me, ne, me, ne, ne, pe };
    int const xe[] = { me, ne, ne, pe, ne, pe, me, ne };

    std::vector<double> lYarr( le + 1 );
    std::vector<double> mXarr( me + 1 );
    std::vector<double> nXarr( ne + 1 );
    std::vector<double> nZarr( ne + 1 );
    std::vector<double> pZarr( pe + 1 );
    std::vector<double> rMxArr( me + 1 );
    std::vector<double> rMzArr( me + 1 );
    std::vector<double> rPxArr( pe + 1 );
    std::vector<double> rPzArr( pe + 1 );

    for ( int li = 0; li <= le; li++ ) {
        lYarr[li] = static_cast<double>( li ) * c / static_cast<double>( le );
    } // end for li

    for ( int mi = 0; mi <= me; mi++ ) {
        mXarr[mi] = static_cast<double>( mi ) * a / static_cast<double>( me );
        double th = ( PI / 2.0 ) -
                    ( static_cast<double>( mi ) * ( PI ) / ( 4.0 * static_cast<double>( me ) ) );
        rMxArr[mi] = r * cos( th );
        rMzArr[mi] = r * sin( th );
    } // end for mi

    for ( int ni = 0; ni <= ne; ni++ ) {
        nXarr[ni] = r + ( static_cast<double>( ni ) * ( a - r ) / static_cast<double>( ne ) );
        nZarr[ni] = r + ( static_cast<double>( ni ) * ( b - r ) / static_cast<double>( ne ) );
    } // end for ni

    for ( int pi = 0; pi <= pe; pi++ ) {
        pZarr[pi]  = static_cast<double>( pi ) * b / static_cast<double>( pe );
        double th  = static_cast<double>( pi ) * ( PI ) / ( 4.0 * static_cast<double>( pe ) );
        rPxArr[pi] = r * cos( th );
        rPzArr[pi] = r * sin( th );
    } // end for pi


    std::vector<std::vector<std::vector<std::vector<int>>>> uniqueNodeId( 8 );
    for ( int ei = 0; ei < 8; ei++ ) {
        uniqueNodeId[ei].resize( le + 1 );
        for ( int li = 0; li <= le; li++ ) {
            uniqueNodeId[ei][li].resize( ze[ei] + 1 );
            for ( int k = 0; k < ( ze[ei] + 1 ); k++ ) {
                uniqueNodeId[ei][li][k].resize( xe[ei] + 1 );
            } // end for k
        } // end for li
    } // end for ei

    int nodeCnt = 0;
    int numPts  = 4 * ( le + 1 ) * ( ne + 1 ) * ( me + pe );
    std::vector<std::array<double, 3>> nodes( numPts );
    for ( int li = 0; li <= le; li++ ) {

        // Node zone 1
        for ( int ni = 0; ni <= ne; ni++ ) {
            uniqueNodeId[0][li][ni][0] = nodeCnt;
            uniqueNodeId[4][li][0][ni] = nodeCnt;
            nodes[nodeCnt++]           = { 0.0, lYarr[li], nZarr[ni] };
        } // end for ni

        // Node zone 2
        for ( int ni = 0; ni <= ne; ni++ ) {
            uniqueNodeId[2][li][0][ni] = nodeCnt;
            uniqueNodeId[6][li][ni][0] = nodeCnt;
            nodes[nodeCnt++]           = { 0.0, lYarr[li], -nZarr[ni] };
        } // end for ni

        // Node zone 3
        for ( int ni = 0; ni <= ne; ni++ ) {
            uniqueNodeId[1][li][0][ni] = nodeCnt;
            uniqueNodeId[3][li][ni][0] = nodeCnt;
            nodes[nodeCnt++]           = { nXarr[ni], lYarr[li], 0.0 };
        } // end for ni

        // Node zone 4
        for ( int ni = 0; ni <= ne; ni++ ) {
            uniqueNodeId[5][li][ni][0] = nodeCnt;
            uniqueNodeId[7][li][0][ni] = nodeCnt;
            nodes[nodeCnt++]           = { -nXarr[ni], lYarr[li], 0.0 };
        } // end for ni

        // Node zone 5
        for ( int mi = 1; mi <= me; mi++ ) {
            for ( int ni = 0; ni <= ne; ni++ ) {
                uniqueNodeId[0][li][ni][mi] = nodeCnt;
                if ( mi == me ) {
                    uniqueNodeId[1][li][pe][ni] = nodeCnt;
                }
                double xPos =
                    rMxArr[mi] + ( ( mXarr[mi] - rMxArr[mi] ) * static_cast<double>( ni ) /
                                   static_cast<double>( ne ) );
                double zPos      = rMzArr[mi] + ( ( b - rMzArr[mi] ) * static_cast<double>( ni ) /
                                             static_cast<double>( ne ) );
                nodes[nodeCnt++] = { xPos, lYarr[li], zPos };
            } // end for ni
        } // end for mi

        // Node zone 6
        for ( int pi = 1; pi < pe; pi++ ) {
            for ( int ni = 0; ni <= ne; ni++ ) {
                uniqueNodeId[1][li][pi][ni] = nodeCnt;
                double xPos = rPxArr[pi] + ( ( a - rPxArr[pi] ) * static_cast<double>( ni ) /
                                             static_cast<double>( ne ) );
                double zPos =
                    rPzArr[pi] + ( ( pZarr[pi] - rPzArr[pi] ) * static_cast<double>( ni ) /
                                   static_cast<double>( ne ) );
                nodes[nodeCnt++] = { xPos, lYarr[li], zPos };
            } // end for ni
        } // end for pi

        // Node zone 7
        for ( int mi = 1; mi <= me; mi++ ) {
            for ( int ni = 0; ni <= ne; ni++ ) {
                uniqueNodeId[2][li][mi][ni] = nodeCnt;
                if ( mi == me )
                    uniqueNodeId[3][li][ni][pe] = nodeCnt;
                double xPos =
                    rMxArr[mi] + ( ( mXarr[mi] - rMxArr[mi] ) * static_cast<double>( ni ) /
                                   static_cast<double>( ne ) );
                double zPos      = rMzArr[mi] + ( ( b - rMzArr[mi] ) * static_cast<double>( ni ) /
                                             static_cast<double>( ne ) );
                nodes[nodeCnt++] = { xPos, lYarr[li], -zPos };
            } // end for ni
        } // end for mi

        // Node zone 8
        for ( int pi = 1; pi < pe; pi++ ) {
            for ( int ni = 0; ni <= ne; ni++ ) {
                uniqueNodeId[3][li][ni][pi] = nodeCnt;
                double xPos = rPxArr[pi] + ( ( a - rPxArr[pi] ) * static_cast<double>( ni ) /
                                             static_cast<double>( ne ) );
                double zPos =
                    rPzArr[pi] + ( ( pZarr[pi] - rPzArr[pi] ) * static_cast<double>( ni ) /
                                   static_cast<double>( ne ) );
                nodes[nodeCnt++] = { xPos, lYarr[li], -zPos };
            } // end for ni
        } // end for pi

        // Node zone 9
        for ( int mi = 1; mi <= me; mi++ ) {
            for ( int ni = 0; ni <= ne; ni++ ) {
                uniqueNodeId[4][li][mi][ni] = nodeCnt;
                if ( mi == me ) {
                    uniqueNodeId[5][li][ni][pe] = nodeCnt;
                }
                double xPos =
                    rMxArr[mi] + ( ( mXarr[mi] - rMxArr[mi] ) * static_cast<double>( ni ) /
                                   static_cast<double>( ne ) );
                double zPos      = rMzArr[mi] + ( ( b - rMzArr[mi] ) * static_cast<double>( ni ) /
                                             static_cast<double>( ne ) );
                nodes[nodeCnt++] = { -xPos, lYarr[li], zPos };
            } // end for ni
        } // end for mi

        // Node zone 10
        for ( int pi = 1; pi < pe; pi++ ) {
            for ( int ni = 0; ni <= ne; ni++ ) {
                uniqueNodeId[5][li][ni][pi] = nodeCnt;
                double xPos = rPxArr[pi] + ( ( a - rPxArr[pi] ) * static_cast<double>( ni ) /
                                             static_cast<double>( ne ) );
                double zPos =
                    rPzArr[pi] + ( ( pZarr[pi] - rPzArr[pi] ) * static_cast<double>( ni ) /
                                   static_cast<double>( ne ) );
                nodes[nodeCnt++] = { -xPos, lYarr[li], zPos };
            } // end for ni
        } // end for pi

        // Node zone 11
        for ( int mi = 1; mi <= me; mi++ ) {
            for ( int ni = 0; ni <= ne; ni++ ) {
                uniqueNodeId[6][li][ni][mi] = nodeCnt;
                if ( mi == me ) {
                    uniqueNodeId[7][li][pe][ni] = nodeCnt;
                }
                double xPos =
                    rMxArr[mi] + ( ( mXarr[mi] - rMxArr[mi] ) * static_cast<double>( ni ) /
                                   static_cast<double>( ne ) );
                double zPos      = rMzArr[mi] + ( ( b - rMzArr[mi] ) * static_cast<double>( ni ) /
                                             static_cast<double>( ne ) );
                nodes[nodeCnt++] = { -xPos, lYarr[li], -zPos };
            } // end for ni
        } // end for mi

        // Node zone 12
        for ( int pi = 1; pi < pe; pi++ ) {
            for ( int ni = 0; ni <= ne; ni++ ) {
                uniqueNodeId[7][li][pi][ni] = nodeCnt;
                double xPos = rPxArr[pi] + ( ( a - rPxArr[pi] ) * static_cast<double>( ni ) /
                                             static_cast<double>( ne ) );
                double zPos =
                    rPzArr[pi] + ( ( pZarr[pi] - rPzArr[pi] ) * static_cast<double>( ni ) /
                                   static_cast<double>( ne ) );
                nodes[nodeCnt++] = { -xPos, lYarr[li], -zPos };
            } // end for ni
        } // end for pi

    } // end for li


    int numElem = 4 * le * ne * ( me + pe );
    std::vector<std::array<int, 8>> elem( numElem );

    int elemCnt = 0;
    for ( int li = 0; li < le; li++ ) {
        for ( int ei = 0; ei < 8; ei++ ) {
            for ( int zi = 0; zi < ze[ei]; zi++ ) {
                for ( int xi = 0; xi < xe[ei]; xi++ ) {
                    int p[8];
                    p[0]            = uniqueNodeId[ei][li][zi][xi];
                    p[1]            = uniqueNodeId[ei][li][zi][xi + 1];
                    p[2]            = uniqueNodeId[ei][li + 1][zi][xi + 1];
                    p[3]            = uniqueNodeId[ei][li + 1][zi][xi];
                    p[4]            = uniqueNodeId[ei][li][zi + 1][xi];
                    p[5]            = uniqueNodeId[ei][li][zi + 1][xi + 1];
                    p[6]            = uniqueNodeId[ei][li + 1][zi + 1][xi + 1];
                    p[7]            = uniqueNodeId[ei][li + 1][zi + 1][xi];
                    elem[elemCnt++] = { p[0], p[1], p[2], p[3], p[4], p[5], p[6], p[7] };
                } // end for xi
            } // end for zi
        } // end for ei
    } // end for li


    std::vector<std::vector<int>> BoundaryNodes( 2 );

    // Top
    for ( int li = 0; li <= le; li++ ) {
        // Elem zone 1
        for ( int mi = 0; mi <= me; mi++ )
            BoundaryNodes[0].push_back( uniqueNodeId[0][li][ne][mi] );

        // Elem zone 5
        for ( int mi = 1; mi < me; mi++ )
            BoundaryNodes[0].push_back( uniqueNodeId[4][li][mi][ne] );

        // Elem zone 5: mi = me
        if ( li < le ) {
            BoundaryNodes[0].push_back( uniqueNodeId[4][li][me][ne] );
        } else {
            BoundaryNodes[0].push_back( uniqueNodeId[4][li][me][ne] );
        }
    } // end for li


    // Bottom
    for ( int li = 0; li <= le; li++ ) {
        // Elem zone 3
        for ( int mi = 0; mi <= me; mi++ )
            BoundaryNodes[1].push_back( uniqueNodeId[2][li][mi][ne] );

        // Elem zone 7
        for ( int mi = 1; mi < me; mi++ )
            BoundaryNodes[1].push_back( uniqueNodeId[6][li][ne][mi] );

        // Elem zone 7: mi = me
        if ( li < le ) {
            BoundaryNodes[1].push_back( uniqueNodeId[6][li][ne][me] );
        } else {
            BoundaryNodes[1].push_back( uniqueNodeId[6][li][ne][me] );
        }
    } // end for li

    std::vector<std::vector<int>> SideIds;

    auto db   = std::make_unique<AMP::Database>();
    auto mesh = db->createAddDatabase( "Mesh" );
    putVector( *mesh, "Nodes", nodes );
    putVector( *mesh, "Elems", elem );
    putVector( *mesh, "BoundaryNodes", BoundaryNodes );
    putVector( *mesh, "SideIds", SideIds );
    db->putScalar( "a", a );
    db->putScalar( "b", b );
    db->putScalar( "c", c );
    db->putScalar( "r", r );
    db->putScalar( "le", le );
    db->putScalar( "me", me );
    db->putScalar( "ne", ne );
    db->putScalar( "pe", pe );
    return db;
}


/********************************************************
 * write a 2 element mesh                                *
 ********************************************************/
DatabasePtr create2elementMesh( double a, int ny, int nz, double Lx, double Ly, double Lz )
{
    double hx = Lx / 2.0;
    double hy = Ly / ( static_cast<double>( ny - 1 ) );
    double hz = Lz / ( static_cast<double>( nz - 1 ) );

    auto NODE = [ny]( int xi, int yi, int zi ) {
        return ( ( xi ) + ( 3 * ( yi ) ) + ( 3 * ( zi ) * ( ny ) ) );
    };

    int numPts = 3 * ny * nz;
    std::vector<std::array<double, 3>> nodes( numPts );
    for ( int zi = 0, pi = 0; zi < nz; zi++ ) {
        for ( int yi = 0; yi < ny; yi++ ) {
            {
                int xi      = 0;
                double px   = xi * hx;
                double py   = yi * hy;
                double pz   = zi * hz;
                nodes[pi++] = { px, py, pz };
            }

            {
                int xi      = 1;
                double pz   = zi * hz;
                double zTmp = ( 0.5 * Lz ) - pz;
                double px   = ( xi * hx ) + ( 2.0 * zTmp * a / Lz );
                double py   = yi * hy;
                nodes[pi++] = { px, py, pz };
            }

            {
                int xi      = 2;
                double px   = xi * hx;
                double py   = yi * hy;
                double pz   = zi * hz;
                nodes[pi++] = { px, py, pz };
            }
        } // end for yi
    } // end for zi

    int numElem = 2 * ( ny - 1 ) * ( nz - 1 );
    std::vector<std::array<int, 8>> elem( numElem );
    for ( int zi = 0, pi = 0; zi < ( nz - 1 ); zi++ ) {
        for ( int yi = 0; yi < ( ny - 1 ); yi++ ) {
            for ( int xi = 0; xi < 2; xi++ ) {
                int p0     = NODE( xi, yi, zi );
                int p1     = NODE( ( xi + 1 ), yi, zi );
                int p2     = NODE( ( xi + 1 ), ( yi + 1 ), zi );
                int p3     = NODE( xi, ( yi + 1 ), zi );
                int p4     = NODE( xi, yi, ( zi + 1 ) );
                int p5     = NODE( ( xi + 1 ), yi, ( zi + 1 ) );
                int p6     = NODE( ( xi + 1 ), ( yi + 1 ), ( zi + 1 ) );
                int p7     = NODE( xi, ( yi + 1 ), ( zi + 1 ) );
                elem[pi++] = { p0, p1, p2, p3, p4, p5, p6, p7 };
            } // end for xi
        } // end for yi
    } // end for zi

    std::vector<std::vector<int>> BoundaryNodes( 4 );

    // x = 0, z = 0 edge
    BoundaryNodes[0] = std::vector<int>( ny );
    for ( int yi = 0; yi < ny; yi++ )
        BoundaryNodes[0][yi] = NODE( 0, yi, 0 );

    // x = 0, z = (nz - 1) edge
    BoundaryNodes[1] = std::vector<int>( ny );
    for ( int yi = 0; yi < ny; yi++ )
        BoundaryNodes[0][yi] = NODE( 0, yi, ( nz - 1 ) );

    // x = 2, z = 0 edge
    BoundaryNodes[2] = std::vector<int>( ny );
    for ( int yi = 0; yi < ny; yi++ )
        BoundaryNodes[0][yi] = NODE( 2, yi, 0 );

    // x = 2, z = (nz - 1) edge
    BoundaryNodes[3] = std::vector<int>( ny );
    for ( int yi = 0; yi < ny; yi++ )
        BoundaryNodes[0][yi] = NODE( 2, yi, ( nz - 1 ) );

    std::vector<std::vector<int>> SideIds;

    auto db   = std::make_unique<AMP::Database>();
    auto mesh = db->createAddDatabase( "Mesh" );
    putVector( *mesh, "Nodes", nodes );
    putVector( *mesh, "Elems", elem );
    putVector( *mesh, "BoundaryNodes", BoundaryNodes );
    putVector( *mesh, "SideIds", SideIds );
    auto params = db->createAddDatabase( "MeshParams" );
    params->putScalar( "a", a );
    params->putScalar( "ny", ny );
    params->putScalar( "nz", nz );
    params->putScalar( "Lx", Lx );
    params->putScalar( "Ly", Ly );
    params->putScalar( "Lz", Lz );
    return db;
}


/********************************************************
 * write a 7 element mesh                                *
 ********************************************************/
DatabasePtr create7elementMesh( int NumBoundaryNodeIds )
{
    std::vector<std::array<double, 3>> nodes( 16 );
    nodes[0]  = { 0.0, 0.0, 0.0 };
    nodes[1]  = { 1.0, 0.0, 0.0 };
    nodes[2]  = { 1.0, 1.0, 0.0 };
    nodes[3]  = { 0.0, 1.0, 0.0 };
    nodes[4]  = { 0.0, 0.0, 1.0 };
    nodes[5]  = { 1.0, 0.0, 1.0 };
    nodes[6]  = { 1.0, 1.0, 1.0 };
    nodes[7]  = { 0.0, 1.0, 1.0 };
    nodes[8]  = { 0.249, 0.342, 0.192 };
    nodes[9]  = { 0.826, 0.288, 0.288 };
    nodes[10] = { 0.850, 0.649, 0.263 };
    nodes[11] = { 0.273, 0.750, 0.230 };
    nodes[12] = { 0.320, 0.186, 0.643 };
    nodes[13] = { 0.677, 0.305, 0.683 };
    nodes[14] = { 0.788, 0.693, 0.644 };
    nodes[15] = { 0.165, 0.745, 0.702 };

    std::vector<std::array<int, 8>> elem( 7 );
    elem[0] = { 0, 1, 2, 3, 8, 9, 10, 11 };
    elem[1] = { 12, 13, 14, 15, 4, 5, 6, 7 };
    elem[2] = { 0, 8, 11, 3, 4, 12, 15, 7 };
    elem[3] = { 9, 1, 2, 10, 13, 5, 6, 14 };
    elem[4] = { 0, 1, 9, 8, 4, 5, 13, 12 };
    elem[5] = { 11, 10, 2, 3, 15, 14, 6, 7 };
    elem[6] = { 8, 9, 10, 11, 12, 13, 14, 15 };

    std::vector<std::vector<int>> BoundaryNodes( NumBoundaryNodeIds );
    if ( NumBoundaryNodeIds == 2 ) {
        BoundaryNodes[0] = { 0, 1, 2, 3 };
        BoundaryNodes[1] = { 4, 5, 6, 7 };
    } else if ( NumBoundaryNodeIds == 8 ) {
        for ( int i = 0; i < 8; i++ )
            BoundaryNodes[i] = { i };
    } else {
        AMP_ERROR( "Not valid" );
    }

    std::vector<std::vector<int>> SideIds;

    auto db   = std::make_unique<AMP::Database>();
    auto mesh = db->createAddDatabase( "Mesh" );
    putVector( *mesh, "Nodes", nodes );
    putVector( *mesh, "Elems", elem );
    putVector( *mesh, "BoundaryNodes", BoundaryNodes );
    putVector( *mesh, "SideIds", SideIds );
    return db;
}


/********************************************************
 * write constrained mesh                                *
 ********************************************************/
DatabasePtr createConstrainedMesh( int nx, int ny, int nz, double Lx, double Ly, double Lz )
{
    auto NODE = [nx, ny]( int xi, int yi, int zi ) { return xi + yi * nx + zi * nx * ny; };

    double hx = Lx / ( static_cast<double>( nx - 1 ) );
    double hy = Ly / ( static_cast<double>( ny - 1 ) );
    double hz = Lz / ( static_cast<double>( nz - 1 ) );

    std::vector<std::array<double, 3>> nodes( nx * ny * nz );
    for ( int zi = 0, pi = 0; zi < nz; zi++ ) {
        for ( int yi = 0; yi < ny; yi++ ) {
            for ( int xi = 0; xi < nx; xi++, pi++ ) {
                nodes[pi] = { xi * hx, yi * hy, zi * hz };
            }
        }
    }


    int numElem = ( nx - 1 ) * ( ny - 1 ) * ( nz - 1 );
    std::vector<std::array<int, 8>> elem( numElem );
    for ( int zi = 0, pi = 0; zi < ( nz - 1 ); zi++ ) {
        for ( int yi = 0; yi < ( ny - 1 ); yi++ ) {
            for ( int xi = 0; xi < ( nx - 1 ); xi++, pi++ ) {
                int p[8];
                p[0]     = NODE( xi, yi, zi );
                p[1]     = NODE( ( xi + 1 ), yi, zi );
                p[2]     = NODE( ( xi + 1 ), ( yi + 1 ), zi );
                p[3]     = NODE( xi, ( yi + 1 ), zi );
                p[4]     = NODE( xi, yi, ( zi + 1 ) );
                p[5]     = NODE( ( xi + 1 ), yi, ( zi + 1 ) );
                p[6]     = NODE( ( xi + 1 ), ( yi + 1 ), ( zi + 1 ) );
                p[7]     = NODE( xi, ( yi + 1 ), ( zi + 1 ) );
                elem[pi] = { p[0], p[1], p[2], p[3], p[4], p[5], p[6], p[7] };
            }
        }
    }

    std::vector<std::vector<int>> BoundaryNodes( 4 );

    // x = Lx surface
    BoundaryNodes[0].resize( ny * nz );
    for ( int zi = 0, cnt = 0; zi < nz; zi++ ) {
        for ( int yi = 0; yi < ny; yi++, cnt++ ) {
            BoundaryNodes[0][cnt] = NODE( ( nx - 1 ), yi, zi );
        }
    }
    // Center of top surface
    BoundaryNodes[1] = { NODE( ( ( nx - 1 ) / 2 ), ( ( ny - 1 ) / 2 ), ( nz - 1 ) ) };
    // Center of bottom surface
    BoundaryNodes[2] = { NODE( ( ( nx - 1 ) / 2 ), ( ( ny - 1 ) / 2 ), 0 ) };
    // Center of bottom front edge
    BoundaryNodes[3] = { NODE( ( ( nx - 1 ) / 2 ), 0, 0 ) };

    std::vector<std::vector<int>> SideIds;

    auto db   = std::make_unique<AMP::Database>();
    auto mesh = db->createAddDatabase( "Mesh" );
    putVector( *mesh, "Nodes", nodes );
    putVector( *mesh, "Elems", elem );
    putVector( *mesh, "BoundaryNodes", BoundaryNodes );
    putVector( *mesh, "SideIds", SideIds );
    db->putScalar( "nx", nx );
    db->putScalar( "ny", ny );
    db->putScalar( "nz", nz );
    db->putScalar( "Lx", Lx );
    db->putScalar( "Ly", Ly );
    db->putScalar( "Lz", Lz );
    return db;
}


/********************************************************
 * write cook mesh (mechanics)                           *
 ********************************************************/
DatabasePtr createCookMesh( int nx, int ny, int nz )
{
    auto NODE = [nx, ny]( int xi, int yi, int zi ) { return xi + yi * nx + zi * nx * ny; };

    double hx = 48.0 / static_cast<double>( nx - 1 );
    double hy = 1.0 / static_cast<double>( ny - 1 );

    int numPts = nx * ny * nz;
    std::vector<std::array<double, 3>> nodes( numPts );
    for ( int zi = 0, pi = 0; zi < nz; zi++ ) {
        for ( int yi = 0; yi < ny; yi++ ) {
            for ( int xi = 0; xi < nx; xi++, pi++ ) {
                double px = xi * hx;
                double py = yi * hy;
                double zl = 44.0 * px / 48.0;
                double zu = 44.0 + ( 16.0 * px / 48.0 );
                double pz = zl + ( zi * ( zu - zl ) / ( nz - 1 ) );
                nodes[pi] = { px, py, pz };
            }
        }
    }

    int numElem = ( nx - 1 ) * ( ny - 1 ) * ( nz - 1 );
    std::vector<std::array<int, 8>> elem( numElem );
    for ( int zi = 0, pi = 0; zi < ( nz - 1 ); zi++ ) {
        for ( int yi = 0; yi < ( ny - 1 ); yi++ ) {
            for ( int xi = 0; xi < ( nx - 1 ); xi++, pi++ ) {
                int p0   = NODE( xi, yi, zi );
                int p1   = NODE( ( xi + 1 ), yi, zi );
                int p2   = NODE( ( xi + 1 ), ( yi + 1 ), zi );
                int p3   = NODE( xi, ( yi + 1 ), zi );
                int p4   = NODE( xi, yi, ( zi + 1 ) );
                int p5   = NODE( ( xi + 1 ), yi, ( zi + 1 ) );
                int p6   = NODE( ( xi + 1 ), ( yi + 1 ), ( zi + 1 ) );
                int p7   = NODE( xi, ( yi + 1 ), ( zi + 1 ) );
                elem[pi] = { p0, p1, p2, p3, p4, p5, p6, p7 };
            }
        }
    }

    std::vector<std::vector<int>> BoundaryNodes( 2 );

    // x = 0 surface
    BoundaryNodes[0].resize( ny * nz );
    for ( int zi = 0, cnt = 0; zi < nz; zi++ ) {
        for ( int yi = 0; yi < ny; yi++, cnt++ ) {
            BoundaryNodes[0][cnt] = NODE( 0, yi, zi );
        }
    }

    // x = Lx surface
    BoundaryNodes[1].resize( ny * nz );
    for ( int zi = 0, cnt = 0; zi < nz; zi++ ) {
        for ( int yi = 0; yi < ny; yi++, cnt++ ) {
            BoundaryNodes[1][cnt] = NODE( ( nx - 1 ), yi, zi );
        } // end for yi
    } // end for zi

    std::vector<std::vector<int>> SideIds;

    auto db   = std::make_unique<AMP::Database>();
    auto mesh = db->createAddDatabase( "Mesh" );
    putVector( *mesh, "Nodes", nodes );
    putVector( *mesh, "Elems", elem );
    putVector( *mesh, "BoundaryNodes", BoundaryNodes );
    putVector( *mesh, "SideIds", SideIds );
    return db;
}


/********************************************************
 * write AMG test mesh                                   *
 ********************************************************/
DatabasePtr createAMGMesh( int nx, int ny, int nz, double Lx, double Ly, double Lz )
{
    auto NODE = [nx, ny]( int xi, int yi, int zi ) { return xi + ( yi * nx ) + ( zi * nx * ny ); };

    double hx = Lx / ( static_cast<double>( nx - 1 ) );
    double hy = Ly / ( static_cast<double>( ny - 1 ) );
    double hz = Lz / ( static_cast<double>( nz - 1 ) );

    int numPts = nx * ny * nz;
    std::vector<std::array<double, 3>> nodes( numPts );
    for ( int zi = 0, pi = 0; zi < nz; zi++ ) {
        for ( int yi = 0; yi < ny; yi++ ) {
            for ( int xi = 0; xi < nx; xi++, pi++ ) {
                nodes[pi] = { xi * hx, yi * hy, zi * hz };
            }
        }
    }

    int numElem = ( nx - 1 ) * ( ny - 1 ) * ( nz - 1 );
    std::vector<std::array<int, 8>> elem( numElem );
    for ( int zi = 0, pi = 0; zi < ( nz - 1 ); zi++ ) {
        for ( int yi = 0; yi < ( ny - 1 ); yi++ ) {
            for ( int xi = 0; xi < ( nx - 1 ); xi++, pi++ ) {
                int p[8];
                p[0]     = NODE( xi, yi, zi );
                p[1]     = NODE( ( xi + 1 ), yi, zi );
                p[2]     = NODE( ( xi + 1 ), ( yi + 1 ), zi );
                p[3]     = NODE( xi, ( yi + 1 ), zi );
                p[4]     = NODE( xi, yi, ( zi + 1 ) );
                p[5]     = NODE( ( xi + 1 ), yi, ( zi + 1 ) );
                p[6]     = NODE( ( xi + 1 ), ( yi + 1 ), ( zi + 1 ) );
                p[7]     = NODE( xi, ( yi + 1 ), ( zi + 1 ) );
                elem[pi] = { p[0], p[1], p[2], p[3], p[4], p[5], p[6], p[7] };
            }
        }
    }

    std::vector<std::vector<int>> BoundaryNodes( 6 );
    // x = 0 surface
    BoundaryNodes[0].resize( ny * nz );
    for ( int zi = 0, cnt = 0; zi < nz; zi++ ) {
        for ( int yi = 0; yi < ny; yi++, cnt++ ) {
            BoundaryNodes[0][cnt] = NODE( 0, yi, zi );
        }
    }
    // x = Lx surface
    BoundaryNodes[1].resize( ny * nz );
    for ( int zi = 0, cnt = 0; zi < nz; zi++ ) {
        for ( int yi = 0; yi < ny; yi++, cnt++ ) {
            BoundaryNodes[1][cnt] = NODE( ( nx - 1 ), yi, zi );
        }
    }
    // y = 0 surface
    BoundaryNodes[2].resize( nx * nz );
    for ( int zi = 0, cnt = 0; zi < nz; zi++ ) {
        for ( int xi = 0; xi < nx; xi++, cnt++ ) {
            BoundaryNodes[2][cnt] = NODE( xi, 0, zi );
        }
    }
    // y = Ly surface
    BoundaryNodes[3].resize( nx * nz );
    for ( int zi = 0, cnt = 0; zi < nz; zi++ ) {
        for ( int xi = 0; xi < nx; xi++, cnt++ ) {
            BoundaryNodes[3][cnt] = NODE( xi, ( ny - 1 ), zi );
        }
    }
    // z = 0 surface
    BoundaryNodes[4].resize( nx * ny );
    for ( int yi = 0, cnt = 0; yi < ny; yi++ ) {
        for ( int xi = 0; xi < nx; xi++, cnt++ ) {
            BoundaryNodes[4][cnt] = NODE( xi, yi, 0 );
        }
    }
    // z = Lz surface
    BoundaryNodes[5].resize( nx * ny );
    for ( int yi = 0, cnt = 0; yi < ny; yi++ ) {
        for ( int xi = 0; xi < nx; xi++, cnt++ ) {
            BoundaryNodes[5][cnt] = NODE( xi, yi, ( nz - 1 ) );
        }
    }

    std::vector<std::vector<int>> SideIds;

    auto db   = std::make_unique<AMP::Database>();
    auto mesh = db->createAddDatabase( "Mesh" );
    putVector( *mesh, "Nodes", nodes );
    putVector( *mesh, "Elems", elem );
    putVector( *mesh, "BoundaryNodes", BoundaryNodes );
    putVector( *mesh, "SideIds", SideIds );
    return db;
}


/********************************************************
 * Read a test mesh file                                 *
 * Note: we cache certain known meshes                   *
 ********************************************************/
std::shared_ptr<AMP::Database> readTestMesh( const std::string &filename, bool useGenerator )
{
    PROFILE( "readTestMesh" );
    if ( useGenerator ) {
        auto db = generateTestMesh( filename );
        if ( db )
            return db;
    }
    // Read the file
    auto db = AMP::Database::parseInputFile( filename );
    return convertDatabase( db );
}
DatabasePtr convertDatabase( std::shared_ptr<AMP::Database> input )
{
    if ( input->getDatabase( "Mesh" )->keyExists( "Nodes" ) )
        return input; // We are using the new format
    PROFILE( "convertDatabase" );
    // Need to convert the old format
    auto output = std::make_shared<AMP::Database>();
    for ( auto key : input->getAllKeys() ) {
        if ( key != "Mesh" )
            output->putData( key, input->getData( key )->clone() );
    }
    auto db1 = input->getDatabase( "Mesh" );
    auto db2 = output->createAddDatabase( "Mesh" );
    char key[100];
    // Convert the nodes
    std::vector<std::array<double, 3>> nodes( db1->getScalar<int>( "NumberOfNodes" ) );
    for ( int i = 0; i < (int) nodes.size(); i++ ) {
        snprintf( key, 100, "Point%d", i );
        auto x   = db1->getVector<double>( key );
        nodes[i] = { x[0], x[1], x[2] };
    }
    putVector( *db2, "Nodes", nodes );
    // Convert the elements
    std::vector<std::array<int, 8>> elem( db1->getScalar<int>( "NumberOfElements" ) );
    for ( int i = 0; i < (int) elem.size(); i++ ) {
        snprintf( key, 100, "Elem%d", i );
        auto x  = db1->getVector<int>( key );
        elem[i] = { x[0], x[1], x[2], x[3], x[4], x[5], x[6], x[7] };
    }
    putVector( *db2, "Elems", elem );
    // Convert the BoundaryNodes
    std::vector<std::vector<int>> BoundaryNodes( db1->getScalar<int>( "NumberOfBoundaryNodeIds" ) );
    for ( int i = 1; i <= (int) BoundaryNodes.size(); i++ ) {
        snprintf( key, 100, "BoundaryNodeId%d", i );
        BoundaryNodes[i - 1] = db1->getWithDefault<std::vector<int>>( key, {}, {} );
    }
    putVector( *db2, "BoundaryNodes", BoundaryNodes );
    // Convert the SideIds
    std::vector<std::vector<int>> SideIds( db1->getScalar<int>( "NumberOfBoundarySideIds" ) );
    for ( int i = 1; i <= (int) SideIds.size(); i++ ) {
        snprintf( key, 100, "BoundarySideId%d", i );
        SideIds[i - 1] = db1->getWithDefault<std::vector<int>>( key, {}, {} );
    }
    putVector( *db2, "SideIds", SideIds );
    return output;
}
std::shared_ptr<AMP::Database> readBinaryTestMesh( const std::string &filename, bool useGenerator )
{
    PROFILE( "readBinaryTestMesh" );

    if ( useGenerator ) {
        auto db = generateTestMesh( filename );
        if ( db )
            return db;
    }

    auto db   = std::make_unique<AMP::Database>();
    auto mesh = db->createAddDatabase( "Mesh" );

    FILE *fp = fopen( filename.c_str(), "rb" );

    int num_nodes;
    fread( &num_nodes, 1, fp );
    std::vector<std::array<double, 3>> nodes( num_nodes );
    fread( nodes[0].data(), 3 * num_nodes, fp );
    putVector( *mesh, "Nodes", nodes );
    int num_elem;
    fread( &num_elem, 1, fp );
    std::vector<std::array<int, 8>> elem( num_elem );
    fread( elem[0].data(), 8 * num_elem, fp );
    putVector( *mesh, "Elems", elem );

    int num_boundaryNodeIds;
    fread( &num_boundaryNodeIds, 1, fp );
    std::vector<std::vector<int>> BoundaryNodes( num_boundaryNodeIds );
    for ( int bid = 1; bid <= num_boundaryNodeIds; bid++ ) {
        int bndDofSize;
        fread( &bndDofSize, 1, fp );
        if ( bndDofSize > 0 ) {
            BoundaryNodes[bid - 1].resize( bndDofSize );
            fread( BoundaryNodes[bid - 1].data(), bndDofSize, fp );
        }
    }
    putVector( *mesh, "BoundaryNodes", BoundaryNodes );

    int num_boundarySideIds;
    fread( &num_boundarySideIds, 1, fp );
    std::vector<std::vector<int>> SideIds( num_boundarySideIds );
    for ( int bid = 1; bid <= num_boundarySideIds; bid++ ) {
        int bndDofSize;
        fread( &bndDofSize, 1, fp );
        if ( bndDofSize > 0 ) {
            SideIds[bid - 1].resize( 2 * bndDofSize );
            fread( SideIds[bid - 1].data(), 2 * bndDofSize, fp );
        }
    }
    putVector( *mesh, "SideIds", SideIds );

    fclose( fp );
    return db;
}


/********************************************************
 * write a test mesh database to file in the same format *
 * as the original files in AMP_DATA                     *
 ********************************************************/
void writeVector( FILE *fp, const char *key, int i, size_t N, const int *x )
{
    if ( N == 0 ) {
        fprintf( fp, "%s%d = \n", key, i );
    } else if ( N == 1 ) {
        fprintf( fp, "%s%d = %d \n", key, i, x[0] );
    } else {
        fprintf( fp, "%s%d = %d", key, i, x[0] );
        for ( size_t j = 1; j < N; j++ )
            fprintf( fp, ", %d", x[j] );
    }
}
void writeTestMesh( const AMP::Database &db0, const std::string &filename )
{
    PROFILE( "writeTestMesh" );
    FILE *fp = fopen( filename.data(), "w" );

    auto db = db0.getDatabase( "Mesh" );
    fprintf( fp, "Mesh { \n" );

    auto nodes = getVector<double, 3>( *db, "Nodes" );
    fprintf( fp, "NumberOfNodes = %i \n", (int) nodes.size() );
    for ( int i = 0; i < (int) nodes.size(); i++ ) {
        auto &x = nodes[i];
        fprintf( fp, "Point%d = %lf, %lf, %lf \n", i, x[0], x[1], x[2] );
    }
    fprintf( fp, "\n" );

    auto elem = getVector<int, 8>( *db, "Elems" );
    fprintf( fp, "NumberOfElements = %i \n", (int) elem.size() );
    for ( int i = 0; i < (int) elem.size(); i++ ) {
        auto &x = elem[i];
        writeVector( fp, "Elem", i, x.size(), x.data() );
        fprintf( fp, " \n" );
    }
    fprintf( fp, "\n" );

    auto BoundaryNodes = getVector( *db, "BoundaryNodes" );
    fprintf( fp, "NumberOfBoundaryNodeIds = %i \n\n", (int) BoundaryNodes.size() );
    for ( int i = 1; i <= (int) BoundaryNodes.size(); i++ ) {
        auto &x = BoundaryNodes[i - 1];
        fprintf( fp, "NumberOfBoundaryNodes%d = %d \n", i, (int) x.size() );
        writeVector( fp, "BoundaryNodeId", i, x.size(), x.data() );
        fprintf( fp, "\n\n" );
    }

    auto SideIds = getVector( *db, "SideIds" );
    fprintf( fp, "NumberOfBoundarySideIds = %i \n\n", (int) SideIds.size() );
    for ( int i = 1; i <= (int) SideIds.size(); i++ ) {
        auto &x = SideIds[i - 1];
        int N   = x.empty() ? 0 : x.size() / 2;
        fprintf( fp, "NumberOfBoundarySides%d = %d \n", i, N );
        writeVector( fp, "BoundarySideId", i, x.size(), x.data() );
        fprintf( fp, "\n\n" );
    }

    fprintf( fp, "} \n" );

    // Print any reference data (may be printed out of order)
    auto keys = db0.getAllKeys();
    if ( keys.size() > 1 )
        fprintf( fp, "\n" );
    for ( auto key : keys ) {
        if ( key == "Mesh" )
            continue;
        auto data = db0.getData( key );
        std::ostringstream os;
        data->print( os );
        fprintf( fp, "%s = %s\n", key.data(), os.str().data() );
    }
    fclose( fp );
}

void writeBinaryTestMesh( const AMP::Database &db0, const std::string &filename )
{
    PROFILE( "writeBinaryTestMesh" );

    auto db            = db0.getDatabase( "Mesh" );
    auto nodes         = getVector<double, 3>( *db, "Nodes" );
    auto elem          = getVector<int, 8>( *db, "Elems" );
    auto BoundaryNodes = getVector( *db, "BoundaryNodes" );
    auto SideIds       = getVector( *db, "SideIds" );

    FILE *fp = fopen( filename.data(), "wb" );

    int num_nodes = nodes.size();
    fwrite( &num_nodes, 1, fp );
    fwrite( nodes[0].data(), 3 * num_nodes, fp );

    int num_elem = elem.size();
    fwrite( &num_elem, 1, fp );
    fwrite( elem[0].data(), 8 * num_elem, fp );

    int num_boundaryNodeIds = BoundaryNodes.size();
    fwrite( &num_boundaryNodeIds, 1, fp );
    for ( int bid = 1; bid <= num_boundaryNodeIds; bid++ ) {
        int bndDofSize = BoundaryNodes[bid - 1].size();
        fwrite( &bndDofSize, 1, fp );
        if ( bndDofSize > 0 )
            fwrite( BoundaryNodes[bid - 1].data(), bndDofSize, fp );
    }

    int num_boundarySideIds = SideIds.size();
    fwrite( &num_boundarySideIds, 1, fp );
    for ( int bid = 1; bid <= num_boundarySideIds; bid++ ) {
        int bndDofSize = SideIds[bid - 1].size();
        fwrite( &bndDofSize, 1, fp );
        if ( bndDofSize > 0 )
            fwrite( SideIds[bid - 1].data(), 2 * bndDofSize, fp );
    }

    fclose( fp );
}


/********************************************************
 * LibMesh generators                                    *
 ********************************************************/
std::shared_ptr<libmeshMesh>
readBinaryTestMeshLibMesh( const std::string &file, const AMP_MPI &comm, const std::string &name )
{
    auto db = readBinaryTestMesh( file );
    return readTestMeshLibMesh( db, comm, name );
}
std::shared_ptr<libmeshMesh>
readTestMeshLibMesh( const std::string &file, const AMP_MPI &comm, const std::string &name )
{
    auto db = readTestMesh( file );
    return readTestMeshLibMesh( db, comm, name );
}
std::shared_ptr<libmeshMesh> readTestMeshLibMesh( std::shared_ptr<AMP::Database> db,
                                                  const AMP_MPI &comm,
                                                  const std::string &name )
{
    AMP_ASSERT( db );
#ifdef AMP_USE_LIBMESH
    [[maybe_unused]] auto libmeshInit =
        std::make_shared<AMP::Mesh::initializeLibMesh>( AMP_COMM_WORLD );
    auto libMeshComm = std::make_shared<libMesh::Parallel::Communicator>( comm.getCommunicator() );
    auto mesh        = std::make_shared<libMesh::Mesh>( *libMeshComm, 3 );
    if ( comm.getRank() == 0 ) {
        db = convertDatabase( db );
        if ( db->keyExists( "Mesh" ) )
            db = db->getDatabase( "Mesh" );
        auto nodes         = getVector<double, 3>( *db, "Nodes" );
        auto elems         = getVector<int, 8>( *db, "Elems" );
        auto BoundaryNodes = getVector( *db, "BoundaryNodes" );
        auto SideIds       = getVector( *db, "SideIds" );
        mesh->reserve_elem( elems.size() );
        mesh->reserve_nodes( nodes.size() );
        for ( int i = 0; i < (int) nodes.size(); i++ ) {
            auto &point = nodes[i];
            mesh->add_point( libMesh::Point( point[0], point[1], point[2] ), i );
        }
        for ( int i = 0; i < (int) elems.size(); i++ ) {
            libMesh::Elem *elem = mesh->add_elem( new libMesh::Hex8 );
            for ( int j = 0; j < 8; j++ )
                elem->set_node( j ) = mesh->node_ptr( elems[i][j] );
        }
        for ( int bid = 1; bid <= (int) BoundaryNodes.size(); bid++ ) {
            auto &bndDofIndices = BoundaryNodes[bid - 1];
            if ( !bndDofIndices.empty() ) {
                for ( int i = 0; i < (int) bndDofIndices.size(); i++ )
                    mesh->boundary_info->add_node( mesh->node_ptr( bndDofIndices[i] ), bid );
            }
        }
        for ( int bid = 1; bid <= (int) SideIds.size(); bid++ ) {
            auto bndDofIndices = SideIds[bid - 1];
            if ( !bndDofIndices.empty() ) {
                int N = bndDofIndices.empty() ? 0 : bndDofIndices.size() / 2;
                for ( int i = 0; i < N; i++ ) {
                    mesh->boundary_info->add_side(
                        mesh->elem_ptr( bndDofIndices[2 * i] ), bndDofIndices[( 2 * i ) + 1], bid );
                }
            }
        }
    }
    libMesh::MeshCommunication().broadcast( *mesh );
    mesh->prepare_for_use( false );
    return std::make_shared<AMP::Mesh::libmeshMesh>( mesh, name, libMeshComm );
#else
    return {};
#endif
}


} // namespace AMP::Mesh::MeshWriters
