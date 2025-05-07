#include "AMP/mesh/testHelpers/meshWriters.h"
#include "AMP/utils/UtilityMacros.h"

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <vector>


static constexpr double PI = 3.1415926535897932;


namespace AMP::Mesh::MeshWriters {


/********************************************************
 * write a box mesh                                      *
 ********************************************************/
void writeBox( int nx, int ny, int nz, double Lx, double Ly, double Lz, const std::string &file )
{

    auto NODE = [nx, ny]( int xi, int yi, int zi ) { return xi + yi * nx + zi * nx * ny; };
    double hx = Lx / ( static_cast<double>( nx - 1 ) );
    double hy = Ly / ( static_cast<double>( ny - 1 ) );
    double hz = Lz / ( static_cast<double>( nz - 1 ) );

    FILE *fp = fopen( file.data(), "w" );
    fprintf( fp, "Mesh { \n" );
    fprintf( fp, "NumberOfNodes = %d \n", ( nx * ny * nz ) );
    for ( int zi = 0, pi = 0; zi < nz; zi++ ) {
        for ( int yi = 0; yi < ny; yi++ ) {
            for ( int xi = 0; xi < nx; xi++, pi++ ) {
                double p[3];
                p[0] = ( static_cast<double>( xi ) ) * hx;
                p[1] = ( static_cast<double>( yi ) ) * hy;
                p[2] = ( static_cast<double>( zi ) ) * hz;
                fprintf( fp, "Point%d = %lf, %lf, %lf \n", pi, p[0], p[1], p[2] );
            } // end for xi
        } // end for yi
    } // end for zi

    fprintf( fp, "\n" );

    int numElem = ( nx - 1 ) * ( ny - 1 ) * ( nz - 1 );
    fprintf( fp, "NumberOfElements = %d \n", numElem );

    for ( int zi = 0, pi = 0; zi < ( nz - 1 ); zi++ ) {
        for ( int yi = 0; yi < ( ny - 1 ); yi++ ) {
            for ( int xi = 0; xi < ( nx - 1 ); xi++, pi++ ) {
                int p[8];
                p[0] = NODE( xi, yi, zi );
                p[1] = NODE( ( xi + 1 ), yi, zi );
                p[2] = NODE( ( xi + 1 ), ( yi + 1 ), zi );
                p[3] = NODE( xi, ( yi + 1 ), zi );
                p[4] = NODE( xi, yi, ( zi + 1 ) );
                p[5] = NODE( ( xi + 1 ), yi, ( zi + 1 ) );
                p[6] = NODE( ( xi + 1 ), ( yi + 1 ), ( zi + 1 ) );
                p[7] = NODE( xi, ( yi + 1 ), ( zi + 1 ) );
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
        }
    }
    fprintf( fp, "\n" );
    fprintf( fp, "NumberOfBoundaryNodeIds = 2 \n\n" );

    // z = Lz surface
    fprintf( fp, "NumberOfBoundaryNodes1 = %d \n", ( nx * ny ) );
    fprintf( fp, "BoundaryNodeId1 = " );
    for ( int yi = 0, cnt = 0; yi < ny; yi++ ) {
        for ( int xi = 0; xi < nx; xi++, cnt++ ) {
            int num = NODE( xi, yi, ( nz - 1 ) );
            fprintf( fp, "%d", num );
            if ( cnt < ( ( nx * ny ) - 1 ) ) {
                fprintf( fp, ", " );
            }
        }
    }
    fprintf( fp, "\n\n" );

    // z = 0 surface
    fprintf( fp, "NumberOfBoundaryNodes2 = %d \n", ( nx * ny ) );
    fprintf( fp, "BoundaryNodeId2 = " );
    for ( int yi = 0, cnt = 0; yi < ny; yi++ ) {
        for ( int xi = 0; xi < nx; xi++, cnt++ ) {
            int num = NODE( xi, yi, 0 );
            fprintf( fp, "%d", num );
            if ( cnt < ( ( nx * ny ) - 1 ) ) {
                fprintf( fp, ", " );
            }
        }
    }
    fprintf( fp, "\n\n" );
    fprintf( fp, "NumberOfBoundarySideIds = 0 \n\n" );
    fprintf( fp, "} \n" );
    fprintf( fp, " Lx = %lf\n Ly = %lf\n Lz = %lf\n", Lx, Ly, Lz );
    fprintf( fp, " nx = %d\n ny = %d\n nz = %d\n", nx, ny, nz );
    fprintf( fp, "\n\n" );
    fclose( fp );
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
    fprintf( fp, "\n\n" );
    fclose( fp );
}


/********************************************************
 * write a distorted element mesh                        *
 ********************************************************/
void writeDistortedElement( const std::string &name )
{
    FILE *fp = fopen( name.data(), "w" );
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


/********************************************************
 * write a plate with a hole                             *
 ********************************************************/
void writePlateWithHole( int le,
                         int me,
                         int ne,
                         int pe,
                         double a,
                         double b,
                         double c,
                         double r,
                         const std::string &file )
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

    FILE *fp = fopen( file.data(), "w" );

    fprintf( fp, "Mesh { \n" );

    int numPts = 4 * ( le + 1 ) * ( ne + 1 ) * ( me + pe );
    fprintf( fp, "NumberOfNodes = %d \n", numPts );

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
    for ( int li = 0; li <= le; li++ ) {

        // Node zone 1
        for ( int ni = 0; ni <= ne; ni++ ) {
            uniqueNodeId[0][li][ni][0] = nodeCnt;
            uniqueNodeId[4][li][0][ni] = nodeCnt;
            fprintf( fp, "Point%d = %lf, %lf, %lf \n", nodeCnt, 0.0, lYarr[li], nZarr[ni] );
            nodeCnt++;
        } // end for ni

        // Node zone 2
        for ( int ni = 0; ni <= ne; ni++ ) {
            uniqueNodeId[2][li][0][ni] = nodeCnt;
            uniqueNodeId[6][li][ni][0] = nodeCnt;
            fprintf( fp, "Point%d = %lf, %lf, %lf \n", nodeCnt, 0.0, lYarr[li], -nZarr[ni] );
            nodeCnt++;
        } // end for ni

        // Node zone 3
        for ( int ni = 0; ni <= ne; ni++ ) {
            uniqueNodeId[1][li][0][ni] = nodeCnt;
            uniqueNodeId[3][li][ni][0] = nodeCnt;
            fprintf( fp, "Point%d = %lf, %lf, %lf \n", nodeCnt, nXarr[ni], lYarr[li], 0.0 );
            nodeCnt++;
        } // end for ni

        // Node zone 4
        for ( int ni = 0; ni <= ne; ni++ ) {
            uniqueNodeId[5][li][ni][0] = nodeCnt;
            uniqueNodeId[7][li][0][ni] = nodeCnt;
            fprintf( fp, "Point%d = %lf, %lf, %lf \n", nodeCnt, -nXarr[ni], lYarr[li], 0.0 );
            nodeCnt++;
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
                double zPos = rMzArr[mi] + ( ( b - rMzArr[mi] ) * static_cast<double>( ni ) /
                                             static_cast<double>( ne ) );
                fprintf( fp, "Point%d = %lf, %lf, %lf \n", nodeCnt, xPos, lYarr[li], zPos );
                nodeCnt++;
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
                fprintf( fp, "Point%d = %lf, %lf, %lf \n", nodeCnt, xPos, lYarr[li], zPos );
                nodeCnt++;
            } // end for ni
        } // end for pi

        // Node zone 7
        for ( int mi = 1; mi <= me; mi++ ) {
            for ( int ni = 0; ni <= ne; ni++ ) {
                uniqueNodeId[2][li][mi][ni] = nodeCnt;
                if ( mi == me ) {
                    uniqueNodeId[3][li][ni][pe] = nodeCnt;
                }
                double xPos =
                    rMxArr[mi] + ( ( mXarr[mi] - rMxArr[mi] ) * static_cast<double>( ni ) /
                                   static_cast<double>( ne ) );
                double zPos = rMzArr[mi] + ( ( b - rMzArr[mi] ) * static_cast<double>( ni ) /
                                             static_cast<double>( ne ) );
                fprintf( fp, "Point%d = %lf, %lf, %lf \n", nodeCnt, xPos, lYarr[li], -zPos );
                nodeCnt++;
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
                fprintf( fp, "Point%d = %lf, %lf, %lf \n", nodeCnt, xPos, lYarr[li], -zPos );
                nodeCnt++;
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
                double zPos = rMzArr[mi] + ( ( b - rMzArr[mi] ) * static_cast<double>( ni ) /
                                             static_cast<double>( ne ) );
                fprintf( fp, "Point%d = %lf, %lf, %lf \n", nodeCnt, -xPos, lYarr[li], zPos );
                nodeCnt++;
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
                fprintf( fp, "Point%d = %lf, %lf, %lf \n", nodeCnt, -xPos, lYarr[li], zPos );
                nodeCnt++;
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
                double zPos = rMzArr[mi] + ( ( b - rMzArr[mi] ) * static_cast<double>( ni ) /
                                             static_cast<double>( ne ) );
                fprintf( fp, "Point%d = %lf, %lf, %lf \n", nodeCnt, -xPos, lYarr[li], -zPos );
                nodeCnt++;
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
                fprintf( fp, "Point%d = %lf, %lf, %lf \n", nodeCnt, -xPos, lYarr[li], -zPos );
                nodeCnt++;
            } // end for ni
        } // end for pi

    } // end for li

    fprintf( fp, "\n" );

    int numElem = 4 * le * ne * ( me + pe );
    fprintf( fp, "NumberOfElements = %d \n", numElem );

    int elemCnt = 0;
    for ( int li = 0; li < le; li++ ) {
        for ( int ei = 0; ei < 8; ei++ ) {
            for ( int zi = 0; zi < ze[ei]; zi++ ) {
                for ( int xi = 0; xi < xe[ei]; xi++ ) {
                    int p[8];
                    p[0] = uniqueNodeId[ei][li][zi][xi];
                    p[1] = uniqueNodeId[ei][li][zi][xi + 1];
                    p[2] = uniqueNodeId[ei][li + 1][zi][xi + 1];
                    p[3] = uniqueNodeId[ei][li + 1][zi][xi];
                    p[4] = uniqueNodeId[ei][li][zi + 1][xi];
                    p[5] = uniqueNodeId[ei][li][zi + 1][xi + 1];
                    p[6] = uniqueNodeId[ei][li + 1][zi + 1][xi + 1];
                    p[7] = uniqueNodeId[ei][li + 1][zi + 1][xi];
                    fprintf( fp,
                             "Elem%d = %d, %d, %d, %d, %d, %d, %d, %d \n",
                             elemCnt,
                             p[0],
                             p[1],
                             p[2],
                             p[3],
                             p[4],
                             p[5],
                             p[6],
                             p[7] );
                    elemCnt++;
                } // end for xi
            } // end for zi
        } // end for ei
    } // end for li

    fprintf( fp, "\n" );

    fprintf( fp, "NumberOfBoundaryNodeIds = 2 \n\n" );

    // Top
    fprintf( fp, "NumberOfBoundaryNodes1 = %d \n", ( ( ( 2 * me ) + 1 ) * ( le + 1 ) ) );
    fprintf( fp, "BoundaryNodeId1 = " );

    for ( int li = 0; li <= le; li++ ) {
        // Elem zone 1
        for ( int mi = 0; mi <= me; mi++ ) {
            fprintf( fp, "%d, ", uniqueNodeId[0][li][ne][mi] );
        } // end for mi

        // Elem zone 5
        for ( int mi = 1; mi < me; mi++ ) {
            fprintf( fp, "%d, ", uniqueNodeId[4][li][mi][ne] );
        } // end for mi

        // Elem zone 5: mi = me
        if ( li < le ) {
            fprintf( fp, "%d, ", uniqueNodeId[4][li][me][ne] );
        } else {
            fprintf( fp, "%d ", uniqueNodeId[4][li][me][ne] );
        }
    } // end for li

    fprintf( fp, "\n\n" );

    // Bottom
    fprintf( fp, "NumberOfBoundaryNodes2 = %d \n", ( ( ( 2 * me ) + 1 ) * ( le + 1 ) ) );
    fprintf( fp, "BoundaryNodeId2 = " );

    for ( int li = 0; li <= le; li++ ) {
        // Elem zone 3
        for ( int mi = 0; mi <= me; mi++ ) {
            fprintf( fp, "%d, ", uniqueNodeId[2][li][mi][ne] );
        } // end for mi

        // Elem zone 7
        for ( int mi = 1; mi < me; mi++ ) {
            fprintf( fp, "%d, ", uniqueNodeId[6][li][ne][mi] );
        } // end for mi

        // Elem zone 7: mi = me
        if ( li < le ) {
            fprintf( fp, "%d, ", uniqueNodeId[6][li][ne][me] );
        } else {
            fprintf( fp, "%d ", uniqueNodeId[6][li][ne][me] );
        }
    } // end for li

    fprintf( fp, "\n\n" );

    fprintf( fp, "NumberOfBoundarySideIds = 0 \n\n" );

    fprintf( fp, "} \n\n" );

    fprintf( fp, " a = %lf\n b = %lf\n c = %lf\n r = %lf\n", a, b, c, r );
    fprintf( fp, " le = %d\n me = %d\n ne = %d\n pe = %d\n", le, me, ne, pe );

    fprintf( fp, "\n\n" );

    fclose( fp );
}


/********************************************************
 * write a 2 element mesh                                *
 ********************************************************/
void write2elementMesh(
    double a, int ny, int nz, double Lx, double Ly, double Lz, const std::string &file )
{
    double hx = Lx / 2.0;
    double hy = Ly / ( static_cast<double>( ny - 1 ) );
    double hz = Lz / ( static_cast<double>( nz - 1 ) );

    auto NODE = [ny]( int xi, int yi, int zi ) {
        return ( ( xi ) + ( 3 * ( yi ) ) + ( 3 * ( zi ) * ( ny ) ) );
    };

    FILE *fp = fopen( file.data(), "w" );

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
    } // end for zi

    fprintf( fp, "\n" );

    int numElem = 2 * ( ny - 1 ) * ( nz - 1 );
    fprintf( fp, "NumberOfElements = %d \n", numElem );

    for ( int zi = 0, pi = 0; zi < ( nz - 1 ); zi++ ) {
        for ( int yi = 0; yi < ( ny - 1 ); yi++ ) {
            for ( int xi = 0; xi < 2; xi++, pi++ ) {
                int p0 = NODE( xi, yi, zi );
                int p1 = NODE( ( xi + 1 ), yi, zi );
                int p2 = NODE( ( xi + 1 ), ( yi + 1 ), zi );
                int p3 = NODE( xi, ( yi + 1 ), zi );
                int p4 = NODE( xi, yi, ( zi + 1 ) );
                int p5 = NODE( ( xi + 1 ), yi, ( zi + 1 ) );
                int p6 = NODE( ( xi + 1 ), ( yi + 1 ), ( zi + 1 ) );
                int p7 = NODE( xi, ( yi + 1 ), ( zi + 1 ) );
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
        } // end for yi
    } // end for zi

    fprintf( fp, "NumberOfBoundaryNodeIds = 4 \n\n" );

    // x = 0, z = 0 edge
    fprintf( fp, "NumberOfBoundaryNodes1 = %d \n", ny );
    fprintf( fp, "BoundaryNodeId1 = " );
    for ( int yi = 0; yi < ny; yi++ ) {
        fprintf( fp, "%d", NODE( 0, yi, 0 ) );
        if ( yi < ( ny - 1 ) ) {
            fprintf( fp, ", " );
        }
    } // end for yi
    fprintf( fp, "\n\n" );

    // x = 0, z = (nz - 1) edge
    fprintf( fp, "NumberOfBoundaryNodes2 = %d \n", ny );
    fprintf( fp, "BoundaryNodeId2 = " );
    for ( int yi = 0; yi < ny; yi++ ) {
        fprintf( fp, "%d", NODE( 0, yi, ( nz - 1 ) ) );
        if ( yi < ( ny - 1 ) ) {
            fprintf( fp, ", " );
        }
    } // end for yi
    fprintf( fp, "\n\n" );

    // x = 2, z = 0 edge
    fprintf( fp, "NumberOfBoundaryNodes3 = %d \n", ny );
    fprintf( fp, "BoundaryNodeId3 = " );
    for ( int yi = 0; yi < ny; yi++ ) {
        fprintf( fp, "%d", NODE( 2, yi, 0 ) );
        if ( yi < ( ny - 1 ) ) {
            fprintf( fp, ", " );
        }
    } // end for yi
    fprintf( fp, "\n\n" );

    // x = 2, z = (nz - 1) edge
    fprintf( fp, "NumberOfBoundaryNodes4 = %d \n", ny );
    fprintf( fp, "BoundaryNodeId4 = " );
    for ( int yi = 0; yi < ny; yi++ ) {
        fprintf( fp, "%d", NODE( 2, yi, ( nz - 1 ) ) );
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


/********************************************************
 * write a 7 element mesh                                *
 ********************************************************/
void write7elementMesh( int NumberOfBoundaryNodeIds, const std::string &file )
{
    FILE *fp = fopen( file.data(), "w" );

    fprintf( fp, "Mesh { \n" );

    fprintf( fp, "NumberOfNodes = 16 \n" );
    fprintf( fp, "Point0 = 0.0, 0.0, 0.0 \n" );
    fprintf( fp, "Point1 = 1.0, 0.0, 0.0 \n" );
    fprintf( fp, "Point2 = 1.0, 1.0, 0.0 \n" );
    fprintf( fp, "Point3 = 0.0, 1.0, 0.0 \n" );
    fprintf( fp, "Point4 = 0.0, 0.0, 1.0 \n" );
    fprintf( fp, "Point5 = 1.0, 0.0, 1.0 \n" );
    fprintf( fp, "Point6 = 1.0, 1.0, 1.0 \n" );
    fprintf( fp, "Point7 = 0.0, 1.0, 1.0 \n" );
    fprintf( fp, "Point8 = 0.249, 0.342, 0.192 \n" );
    fprintf( fp, "Point9 = 0.826, 0.288, 0.288 \n" );
    fprintf( fp, "Point10 = 0.850, 0.649, 0.263 \n" );
    fprintf( fp, "Point11 = 0.273, 0.750, 0.230 \n" );
    fprintf( fp, "Point12 = 0.320, 0.186, 0.643 \n" );
    fprintf( fp, "Point13 = 0.677, 0.305, 0.683 \n" );
    fprintf( fp, "Point14 = 0.788, 0.693, 0.644 \n" );
    fprintf( fp, "Point15 = 0.165, 0.745, 0.702 \n" );

    fprintf( fp, "\n" );

    fprintf( fp, "NumberOfElements = 7 \n" );

    fprintf( fp, "Elem0 = 0, 1, 2, 3, 8, 9, 10, 11 \n" );
    fprintf( fp, "Elem1 = 12, 13, 14, 15, 4, 5, 6, 7 \n" );
    fprintf( fp, "Elem2 = 0, 8, 11, 3, 4, 12, 15, 7 \n" );
    fprintf( fp, "Elem3 = 9, 1, 2, 10, 13, 5, 6, 14 \n" );
    fprintf( fp, "Elem4 = 0, 1, 9, 8, 4, 5, 13, 12 \n" );
    fprintf( fp, "Elem5 = 11, 10, 2, 3, 15, 14, 6, 7 \n" );
    fprintf( fp, "Elem6 = 8, 9, 10, 11, 12, 13, 14, 15 \n" );

    fprintf( fp, "\n" );

    if ( NumberOfBoundaryNodeIds == 2 ) {
        fprintf( fp, "NumberOfBoundaryNodeIds = 2 \n\n" );
        // z = 0 surface
        fprintf( fp, "NumberOfBoundaryNodes1 = 4 \n" );
        fprintf( fp, "BoundaryNodeId1 = 0, 1, 2, 3" );
        fprintf( fp, "\n\n" );
        // z = 1 surface
        fprintf( fp, "NumberOfBoundaryNodes2 = 4 \n" );
        fprintf( fp, "BoundaryNodeId2 = 4, 5, 6, 7" );
        fprintf( fp, "\n\n" );
    } else if ( NumberOfBoundaryNodeIds == 8 ) {
        fprintf( fp, "NumberOfBoundaryNodeIds = 8 \n\n" );
        for ( int i = 0; i < 8; i++ ) {
            fprintf( fp, "NumberOfBoundaryNodes%d = 1 \n", ( i + 1 ) );
            fprintf( fp, "BoundaryNodeId%d = %d", ( i + 1 ), i );
            fprintf( fp, "\n\n" );
        }
    } else {
        AMP_ERROR( "Not valid" );
    }
    fprintf( fp, "NumberOfBoundarySideIds = 0\n\n" );

    fprintf( fp, "}\n" );

    fclose( fp );
}


/********************************************************
 * write constrained mesh                                *
 ********************************************************/
void writeConstrainedMesh(
    int nx, int ny, int nz, double Lx, double Ly, double Lz, const std::string &file )
{
    auto NODE = [nx, ny]( int xi, int yi, int zi ) { return xi + yi * nx + zi * nx * ny; };

    double hx = Lx / ( static_cast<double>( nx - 1 ) );
    double hy = Ly / ( static_cast<double>( ny - 1 ) );
    double hz = Lz / ( static_cast<double>( nz - 1 ) );

    FILE *fp = fopen( file.data(), "w" );

    fprintf( fp, "Mesh { \n" );

    fprintf( fp, "NumberOfNodes = %d \n", ( nx * ny * nz ) );

    for ( int zi = 0, pi = 0; zi < nz; zi++ ) {
        for ( int yi = 0; yi < ny; yi++ ) {
            for ( int xi = 0; xi < nx; xi++, pi++ ) {
                double p[3];
                p[0] = ( static_cast<double>( xi ) ) * hx;
                p[1] = ( static_cast<double>( yi ) ) * hy;
                p[2] = ( static_cast<double>( zi ) ) * hz;
                fprintf( fp, "Point%d = %lf, %lf, %lf \n", pi, p[0], p[1], p[2] );
            } // end for xi
        } // end for yi
    } // end for zi

    fprintf( fp, "\n" );

    int numElem = ( nx - 1 ) * ( ny - 1 ) * ( nz - 1 );
    fprintf( fp, "NumberOfElements = %d \n", numElem );

    for ( int zi = 0, pi = 0; zi < ( nz - 1 ); zi++ ) {
        for ( int yi = 0; yi < ( ny - 1 ); yi++ ) {
            for ( int xi = 0; xi < ( nx - 1 ); xi++, pi++ ) {
                int p[8];
                p[0] = NODE( xi, yi, zi );
                p[1] = NODE( ( xi + 1 ), yi, zi );
                p[2] = NODE( ( xi + 1 ), ( yi + 1 ), zi );
                p[3] = NODE( xi, ( yi + 1 ), zi );
                p[4] = NODE( xi, yi, ( zi + 1 ) );
                p[5] = NODE( ( xi + 1 ), yi, ( zi + 1 ) );
                p[6] = NODE( ( xi + 1 ), ( yi + 1 ), ( zi + 1 ) );
                p[7] = NODE( xi, ( yi + 1 ), ( zi + 1 ) );
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
        } // end for yi
    } // end for zi

    fprintf( fp, "\n" );

    fprintf( fp, "NumberOfBoundaryNodeIds = 4 \n\n" );

    // x = Lx surface
    fprintf( fp, "NumberOfBoundaryNodes1 = %d \n", ( ny * nz ) );
    fprintf( fp, "BoundaryNodeId1 = " );

    for ( int zi = 0, cnt = 0; zi < nz; zi++ ) {
        for ( int yi = 0; yi < ny; yi++, cnt++ ) {
            int num = NODE( ( nx - 1 ), yi, zi );
            fprintf( fp, "%d", num );
            if ( cnt < ( ( ny * nz ) - 1 ) ) {
                fprintf( fp, ", " );
            }
        } // end for yi
    } // end for zi

    fprintf( fp, "\n\n" );

    // Center of top surface
    fprintf( fp, "NumberOfBoundaryNodes2 = 1 \n" );
    fprintf( fp,
             "BoundaryNodeId2 = %d \n\n",
             NODE( ( ( nx - 1 ) / 2 ), ( ( ny - 1 ) / 2 ), ( nz - 1 ) ) );

    // Center of bottom surface
    fprintf( fp, "NumberOfBoundaryNodes3 = 1 \n" );
    fprintf( fp, "BoundaryNodeId3 = %d \n\n", NODE( ( ( nx - 1 ) / 2 ), ( ( ny - 1 ) / 2 ), 0 ) );

    // Center of bottom front edge
    fprintf( fp, "NumberOfBoundaryNodes4 = 1 \n" );
    fprintf( fp, "BoundaryNodeId4 = %d \n\n", NODE( ( ( nx - 1 ) / 2 ), 0, 0 ) );

    fprintf( fp, "NumberOfBoundarySideIds = 0 \n\n" );

    fprintf( fp, "} \n" );

    fprintf( fp, " Lx = %lf\n Ly = %lf\n Lz = %lf\n", Lx, Ly, Lz );
    fprintf( fp, " nx = %d\n ny = %d\n nz = %d\n", nx, ny, nz );

    fprintf( fp, "\n\n" );

    fclose( fp );
}


/********************************************************
 * write cook mesh (mechanics)                           *
 ********************************************************/
void writeCookMesh( int nx, int ny, int nz, const std::string &file )
{
    auto NODE = [nx, ny]( int xi, int yi, int zi ) { return xi + yi * nx + zi * nx * ny; };

    double hx = 48.0 / static_cast<double>( nx - 1 );
    double hy = 1.0 / static_cast<double>( ny - 1 );

    FILE *fp = fopen( file.data(), "w" );

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
        } // end for yi
    } // end for zi

    fprintf( fp, "\n" );

    int numElem = ( nx - 1 ) * ( ny - 1 ) * ( nz - 1 );
    fprintf( fp, "NumberOfElements = %d \n", numElem );

    for ( int zi = 0, pi = 0; zi < ( nz - 1 ); zi++ ) {
        for ( int yi = 0; yi < ( ny - 1 ); yi++ ) {
            for ( int xi = 0; xi < ( nx - 1 ); xi++, pi++ ) {
                int p0 = NODE( xi, yi, zi );
                int p1 = NODE( ( xi + 1 ), yi, zi );
                int p2 = NODE( ( xi + 1 ), ( yi + 1 ), zi );
                int p3 = NODE( xi, ( yi + 1 ), zi );
                int p4 = NODE( xi, yi, ( zi + 1 ) );
                int p5 = NODE( ( xi + 1 ), yi, ( zi + 1 ) );
                int p6 = NODE( ( xi + 1 ), ( yi + 1 ), ( zi + 1 ) );
                int p7 = NODE( xi, ( yi + 1 ), ( zi + 1 ) );
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
        } // end for yi
    } // end for zi

    fprintf( fp, "\n" );

    fprintf( fp, "NumberOfBoundaryNodeIds = 2 \n\n" );

    // x = 0 surface
    fprintf( fp, "NumberOfBoundaryNodes1 = %d \n", ( ny * nz ) );
    fprintf( fp, "BoundaryNodeId1 = " );
    for ( int zi = 0, cnt = 0; zi < nz; zi++ ) {
        for ( int yi = 0; yi < ny; yi++, cnt++ ) {
            fprintf( fp, "%d", NODE( 0, yi, zi ) );
            if ( cnt < ( ( ny * nz ) - 1 ) ) {
                fprintf( fp, ", " );
            }
        } // end for yi
    } // end for zi
    fprintf( fp, "\n\n" );

    // x = Lx surface
    fprintf( fp, "NumberOfBoundaryNodes2 = %d \n", ( ny * nz ) );
    fprintf( fp, "BoundaryNodeId2 = " );
    for ( int zi = 0, cnt = 0; zi < nz; zi++ ) {
        for ( int yi = 0; yi < ny; yi++, cnt++ ) {
            fprintf( fp, "%d", NODE( ( nx - 1 ), yi, zi ) );
            if ( cnt < ( ( ny * nz ) - 1 ) ) {
                fprintf( fp, ", " );
            }
        } // end for yi
    } // end for zi
    fprintf( fp, "\n\n" );
    fprintf( fp, "NumberOfBoundarySideIds = 0 \n\n" );
    fprintf( fp, "} \n" );
    fclose( fp );
}


/********************************************************
 * write AMG test mesh                                   *
 ********************************************************/
void writeAMGMesh(
    int nx, int ny, int nz, double Lx, double Ly, double Lz, const std::string &filename )
{
    auto NODE = [nx, ny]( int xi, int yi, int zi ) { return xi + ( yi * nx ) + ( zi * nx * ny ); };

    double hx = Lx / ( static_cast<double>( nx - 1 ) );
    double hy = Ly / ( static_cast<double>( ny - 1 ) );
    double hz = Lz / ( static_cast<double>( nz - 1 ) );

    auto fp = fopen( filename.data(), "w" );

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
        } // end for yi
    } // end for zi

    fprintf( fp, "\n" );

    int numElem = ( nx - 1 ) * ( ny - 1 ) * ( nz - 1 );
    fprintf( fp, "NumberOfElements = %d \n", numElem );

    for ( int zi = 0, pi = 0; zi < ( nz - 1 ); zi++ ) {
        for ( int yi = 0; yi < ( ny - 1 ); yi++ ) {
            for ( int xi = 0; xi < ( nx - 1 ); xi++, pi++ ) {
                int p[8];
                p[0] = NODE( xi, yi, zi );
                p[1] = NODE( ( xi + 1 ), yi, zi );
                p[2] = NODE( ( xi + 1 ), ( yi + 1 ), zi );
                p[3] = NODE( xi, ( yi + 1 ), zi );
                p[4] = NODE( xi, yi, ( zi + 1 ) );
                p[5] = NODE( ( xi + 1 ), yi, ( zi + 1 ) );
                p[6] = NODE( ( xi + 1 ), ( yi + 1 ), ( zi + 1 ) );
                p[7] = NODE( xi, ( yi + 1 ), ( zi + 1 ) );
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
        } // end for yi
    } // end for zi

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
            int num = NODE( 0, yi, zi );
            fprintf( fp, "%d", num );
            if ( cnt < ( ( ny * nz ) - 1 ) ) {
                fprintf( fp, ", " );
            }
        } // end for yi
    } // end for zi

    fprintf( fp, "\n\n" );

    // x = Lx surface
    {
        int num = ( ny * nz );
        fprintf( fp, "NumberOfBoundaryNodes2 = %d \n", num );
        fprintf( fp, "BoundaryNodeId2 = " );
    }

    for ( int zi = 0, cnt = 0; zi < nz; zi++ ) {
        for ( int yi = 0; yi < ny; yi++, cnt++ ) {
            int num = NODE( ( nx - 1 ), yi, zi );
            fprintf( fp, "%d", num );
            if ( cnt < ( ( ny * nz ) - 1 ) ) {
                fprintf( fp, ", " );
            }
        } // end for yi
    } // end for zi

    fprintf( fp, "\n\n" );

    // y = 0 surface
    {
        int num = ( nx * nz );
        fprintf( fp, "NumberOfBoundaryNodes3 = %d \n", num );
        fprintf( fp, "BoundaryNodeId3 = " );
    }

    for ( int zi = 0, cnt = 0; zi < nz; zi++ ) {
        for ( int xi = 0; xi < nx; xi++, cnt++ ) {
            int num = NODE( xi, 0, zi );
            fprintf( fp, "%d", num );
            if ( cnt < ( ( nx * nz ) - 1 ) ) {
                fprintf( fp, ", " );
            }
        } // end for yi
    } // end for zi

    fprintf( fp, "\n\n" );

    // y = Ly surface
    {
        int num = ( nx * nz );
        fprintf( fp, "NumberOfBoundaryNodes4 = %d \n", num );
        fprintf( fp, "BoundaryNodeId4 = " );
    }

    for ( int zi = 0, cnt = 0; zi < nz; zi++ ) {
        for ( int xi = 0; xi < nx; xi++, cnt++ ) {
            int num = NODE( xi, ( ny - 1 ), zi );
            fprintf( fp, "%d", num );
            if ( cnt < ( ( nx * nz ) - 1 ) ) {
                fprintf( fp, ", " );
            }
        } // end for yi
    } // end for zi

    fprintf( fp, "\n\n" );

    // z = 0 surface
    {
        int num = ( ny * nx );
        fprintf( fp, "NumberOfBoundaryNodes5 = %d \n", num );
        fprintf( fp, "BoundaryNodeId5 = " );
    }

    for ( int yi = 0, cnt = 0; yi < ny; yi++ ) {
        for ( int xi = 0; xi < nx; xi++, cnt++ ) {
            int num = NODE( xi, yi, 0 );
            fprintf( fp, "%d", num );
            if ( cnt < ( ( ny * nx ) - 1 ) ) {
                fprintf( fp, ", " );
            }
        } // end for xi
    } // end for yi

    fprintf( fp, "\n\n" );

    // z = Lz surface
    {
        int num = ( ny * nx );
        fprintf( fp, "NumberOfBoundaryNodes6 = %d \n", num );
        fprintf( fp, "BoundaryNodeId6 = " );
    }

    for ( int yi = 0, cnt = 0; yi < ny; yi++ ) {
        for ( int xi = 0; xi < nx; xi++, cnt++ ) {
            int num = NODE( xi, yi, ( nz - 1 ) );
            fprintf( fp, "%d", num );
            if ( cnt < ( ( ny * nx ) - 1 ) ) {
                fprintf( fp, ", " );
            }
        } // end for xi
    } // end for yi

    fprintf( fp, "\n\n" );

    {
        int num = 0;
        fprintf( fp, "NumberOfBoundarySideIds = %d \n\n", num );
    }

    fprintf( fp, "} \n" );

    fclose( fp );
}


} // namespace AMP::Mesh::MeshWriters
