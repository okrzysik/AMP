#include "AMP/utils/NearestPairSearch.h"
#include "AMP/utils/Utilities.h"
#include <cmath>


// subroutine to find the minimum distance between any two points
static std::pair<int, int>
find_min_dist2( const int ndim, const int N, const double *const *x, int split, double &dist )
{
    const int N_brute = 10; // Number of points below which we will switch to a brute force search
    if ( N < N_brute ) {
        // Use a brute force method
        double dist1 = 1e100;
        std::pair<int, int> index;
        for ( int i = 0; i < N; i++ ) {
            for ( int j = i + 1; j < N; j++ ) {
                double dist2 = 0.0;
                for ( int d = 0; d < ndim; d++ )
                    dist2 += ( x[d][i] - x[d][j] ) * ( x[d][i] - x[d][j] );
                if ( dist2 < dist1 ) {
                    dist1        = dist2;
                    index.first  = i;
                    index.second = j;
                }
            }
        }
        dist = sqrt( dist1 );
        return index;
    }
    // First get the points sorted along the desired direction
    double *x2[5] = { nullptr, nullptr, nullptr, nullptr, nullptr };
    for ( int d = 0; d < ndim; d++ )
        x2[d] = new double[N];
    for ( int i = 0; i < N; i++ )
        x2[split][i] = x[split][i];
    auto I = new int[N];
    for ( int i = 0; i < N; i++ )
        I[i] = i;
    AMP::Utilities::quicksort( N, x2[split], I );
    for ( int d = 0; d < ndim; d++ ) {
        if ( d != split ) {
            for ( int i = 0; i < N; i++ )
                x2[d][i] = x[d][I[i]];
        }
    }
    // Divide the points in half and recursively find the closest pair in each half
    double median = x2[split][N / 2];
    int k         = N / 2 + 1;
    while ( x2[split][k] <= median ) {
        k++;
        if ( k == N ) {
            break;
        }
    }
    if ( N - k < N_brute ) {
        // We could not evenly divide the points, switch to a brute force method
        for ( int d = 0; d < ndim; d++ )
            delete[] x2[d];
        delete[] I;
        double dist1 = 1e100;
        std::pair<int, int> index;
        for ( int i = 0; i < N; i++ ) {
            for ( int j = i + 1; j < N; j++ ) {
                double dist2 = 0.0;
                for ( int d = 0; d < ndim; d++ )
                    dist2 += ( x[d][i] - x[d][j] ) * ( x[d][i] - x[d][j] );
                if ( dist2 < dist1 ) {
                    dist1        = dist2;
                    index.first  = i;
                    index.second = j;
                }
            }
        }
        dist = sqrt( dist1 );
        return index;
    } else {
        median =
            0.5 * ( x2[split][k - 1] + x2[split][k] ); // Move the median to the middle of the split
    }
    double dx     = x2[split][k] - x2[split][k - 1];
    double *x3[5] = { nullptr, nullptr, nullptr, nullptr, nullptr };
    for ( int d = 0; d < ndim; d++ )
        x3[d] = &x2[d][k];
    int split2   = ( split + 1 ) % ndim;
    double dist1 = 0, dist2 = 0;
    std::pair<int, int> pair1 = find_min_dist2( ndim, k, x2, split2, dist1 );
    std::pair<int, int> pair2 = find_min_dist2( ndim, N - k, x3, split2, dist2 );
    pair2.first += k;
    pair2.second += k;
    std::pair<int, int> index;
    if ( dist1 <= dist2 ) {
        dist  = dist1;
        index = pair1;
    } else {
        dist  = dist2;
        index = pair2;
    }
    // Check for pairs of points that are closer and span the median
    auto i1 = static_cast<int>( AMP::Utilities::findfirst(
        N, x2[split], median - dist + dx / 2 ) ); // We only care about points > median-dist+dx/2
    auto i2 = static_cast<int>( AMP::Utilities::findfirst(
        N, x2[split], median + dist - dx / 2 ) ); // We only care about points < median+dist-dx/2
    if ( i2 - i1 <= 1 ) {
        // No pairs can be closer
    } else if ( i2 - i1 > 0.8 * N || i2 - i1 < N_brute ) {
        // We need to check most of the remaining points, use a brute force method
        for ( int i = i1; i < i2; i++ ) {
            for ( int j = i + 1; j < i2; j++ ) {
                double dist3 = 0.0;
                for ( int d = 0; d < ndim; d++ )
                    dist3 += ( x2[d][i] - x2[d][j] ) * ( x2[d][i] - x2[d][j] );
                if ( dist3 < dist * dist ) {
                    dist         = sqrt( dist3 );
                    index.first  = i;
                    index.second = j;
                }
            }
        }
    } else {
        // We can search a small subset of points recursively
        for ( int d = 0; d < ndim; d++ )
            x3[d] = &x2[d][i1];
        double dist3               = 1e100;
        std::pair<int, int> index3 = find_min_dist2( ndim, i2 - i1, x3, 0, dist3 );
        if ( dist3 < dist ) {
            dist         = dist3;
            index.first  = index3.first + i1;
            index.second = index3.second + i1;
        }
    }
    // Convert the index to the original x positions
    index.first  = I[index.first];
    index.second = I[index.second];
    // Free memory and return
    for ( int d = 0; d < ndim; d++ )
        delete[] x2[d];
    delete[] I;
    return index;
}


std::pair<int, int> find_min_dist( const int ndim, const int N, const double *x )
{
    if ( N < 2 )
        return std::pair<int, int>( 0, 0 );
    if ( ndim == 1 ) {
        // Special case for 1D
        auto I  = new int[N];
        auto x2 = new double[N];
        for ( int i = 0; i < N; i++ )
            I[i] = i;
        for ( int i = 0; i < N; i++ )
            x2[i] = x[i];
        AMP::Utilities::quicksort( N, x2, I );
        int index   = 1;
        double dist = x2[1] - x2[0];
        for ( int i = 2; i < N; i++ ) {
            if ( ( x2[i] - x2[i - 1] ) < dist ) {
                index = i;
                dist  = x2[i] - x2[i - 1];
            }
        }
        std::pair<int, int> tmp( I[index - 1], I[index] );
        delete[] I;
        delete[] x2;
        return tmp;
    }
    double *x2[5] = { nullptr, nullptr, nullptr, nullptr, nullptr };
    for ( int d = 0; d < ndim; d++ ) {
        x2[d] = new double[N];
        for ( int i = 0; i < N; i++ )
            x2[d][i] = x[d + i * ndim];
    }
    double dist               = 0;
    std::pair<int, int> index = find_min_dist2( ndim, N, x2, 0, dist );
    for ( int d = 0; d < ndim; d++ )
        delete[] x2[d];
    return index;
}
