
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <vector>

#define __PI__ 3.1415926535897932

int main( int argc, char **argv )
{

    if ( argc < 10 ) {
        std::cout << "Usage: exe le me ne pe a b c r output_file" << std::endl;
        return 1;
    }

    int const le = atoi( argv[1] );
    int const me = atoi( argv[2] );
    int const ne = atoi( argv[3] );
    int const pe = atoi( argv[4] );

    int const ze[] = { ne, pe, me, ne, me, ne, ne, pe };
    int const xe[] = { me, ne, ne, pe, ne, pe, me, ne };

    double const a = atof( argv[5] ); // X-dimension
    double const b = atof( argv[6] ); // Z-dimension
    double const c = atof( argv[7] ); // Y-dimension
    double const r = atof( argv[8] );

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
        mXarr[mi]  = static_cast<double>( mi ) * a / static_cast<double>( me );
        double th  = ( __PI__ / 2.0 ) - ( static_cast<double>( mi ) * ( __PI__ ) /
                                         ( 4.0 * static_cast<double>( me ) ) );
        rMxArr[mi] = r * cos( th );
        rMzArr[mi] = r * sin( th );
    } // end for mi

    for ( int ni = 0; ni <= ne; ni++ ) {
        nXarr[ni] = r + ( static_cast<double>( ni ) * ( a - r ) / static_cast<double>( ne ) );
        nZarr[ni] = r + ( static_cast<double>( ni ) * ( b - r ) / static_cast<double>( ne ) );
    } // end for ni

    for ( int pi = 0; pi <= pe; pi++ ) {
        pZarr[pi]  = static_cast<double>( pi ) * b / static_cast<double>( pe );
        double th  = static_cast<double>( pi ) * ( __PI__ ) / ( 4.0 * static_cast<double>( pe ) );
        rPxArr[pi] = r * cos( th );
        rPzArr[pi] = r * sin( th );
    } // end for pi

    FILE *fp = fopen( argv[9], "w" );

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
        }     // end for li
    }         // end for ei

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
        }     // end for mi

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
        }     // end for pi

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
        }     // end for mi

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
        }     // end for pi

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
        }     // end for mi

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
        }     // end for pi

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
        }     // end for mi

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
        }     // end for pi

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
            }     // end for zi
        }         // end for ei
    }             // end for li

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
