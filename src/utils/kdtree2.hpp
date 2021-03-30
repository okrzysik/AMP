#ifndef included_AMP_kdtree2_hpp
#define included_AMP_kdtree2_hpp

#include "AMP/utils/Utilities.h"
#include "AMP/utils/kdtree2.h"

#include "ProfilerApp.h"

#include <cstring>
#include <memory>


namespace AMP {


/********************************************************
 *  External instantiations                              *
 ********************************************************/
extern template class kdtree2<1, int>;
extern template class kdtree2<2, int>;
extern template class kdtree2<3, int>;


/********************************************************
 * Constructor                                           *
 ********************************************************/
template<uint8_t NDIM, class TYPE>
kdtree2<NDIM, TYPE>::kdtree2( size_t N, const std::array<double, NDIM> *x, const TYPE *data )
    : d_N( N ),
      d_split_dim( 0 ),
      d_split( 0 ),
      d_left( nullptr ),
      d_right( nullptr ),
      d_data( nullptr )
{
    // Update the box
    d_lb.fill( 1e100 );
    d_ub.fill( -1e100 );
    for ( size_t i = 0; i < N; i++ ) {
        for ( int d = 0; d < NDIM; d++ ) {
            d_lb[d] = std::min( d_lb[d], x[i][d] );
            d_ub[d] = std::max( d_ub[d], x[i][d] );
        }
    }
    // If we have more than the threshold split the tree
    constexpr uint64_t threshold = 40; // Optimize for performance
    if ( N > threshold ) {
        // Split the tree and recurse
        splitData( N, x, data );
    } else {
        // Store the data
        d_data       = new data_struct;
        d_data->x    = std::vector<Point>( x, x + N );
        d_data->data = std::vector<TYPE>( data, data + N );
    }
}
template<uint8_t NDIM, class TYPE>
void kdtree2<NDIM, TYPE>::splitData( size_t N, const Point *x, const TYPE *data )
{
    // Choose the splitting direction (and tolerance)
    int dir    = 0;
    double tmp = d_ub[0] - d_lb[0];
    for ( int d = 0; d < NDIM; d++ ) {
        if ( ( d_ub[d] - d_lb[d] ) > tmp ) {
            dir = d;
            tmp = d_ub[d] - d_lb[d];
        }
    }
    // Resort the data along the splitting direction
    auto x2 = new double[N];
    for ( size_t i = 0; i < N; i++ )
        x2[i] = x[i][dir];
    auto index = new uint64_t[N];
    for ( size_t i = 0; i < N; i++ )
        index[i] = i;
    AMP::Utilities::quicksort( N, x2, index );
    auto t1 = new Point[N];
    auto t2 = new TYPE[N];
    for ( size_t i = 0; i < N; i++ ) {
        t1[i] = x[index[i]];
        t2[i] = data[index[i]];
    }
    delete[] index;
    // Find the ideal point to split
    size_t k = find_split( N, x2 );
    AMP_ASSERT( k > 0 );
    // Recursively split
    d_split_dim = dir;
    d_split     = 0.5 * ( x2[k - 1] + x2[k] );
    delete[] x2;
    d_left  = new kdtree2( k, t1, t2 );
    d_right = new kdtree2( N - k, &t1[k], &t2[k] );
    delete[] t1;
    delete[] t2;
}
template<uint8_t NDIM, class TYPE>
kdtree2<NDIM, TYPE>::~kdtree2()
{
    delete d_left;
    delete d_right;
    delete d_data;
}


/********************************************************
 * Return the domain box                                 *
 ********************************************************/
template<uint8_t NDIM, class TYPE>
void kdtree2<NDIM, TYPE>::add( const Point &p, const TYPE &data )
{
    d_N++;
    // Update the bounding box
    for ( int d = 0; d < NDIM; d++ ) {
        d_lb[d] = std::min( d_lb[d], p[d] );
        d_ub[d] = std::max( d_ub[d], p[d] );
    }
    if ( d_left ) {
        // Figure out which half we belong to and split
        if ( p[d_split_dim] <= d_split )
            d_left->add( p, data );
        else
            d_right->add( p, data );
    } else {
        // Add the point to the current leaf node
        d_data->x.push_back( p );
        d_data->data.push_back( data );
        // Split the leaf node if needed
        constexpr uint64_t threshold = 40; // Optimize for performance
        if ( d_N > threshold ) {
            splitData( d_N, d_data->x.data(), d_data->data.data() );
            delete d_data;
            d_data = nullptr;
        }
    }
}


/********************************************************
 * Return the domain box                                 *
 ********************************************************/
template<uint8_t NDIM, class TYPE>
std::array<double, 2 * NDIM> kdtree2<NDIM, TYPE>::box() const
{
    std::array<double, 2 * NDIM> b;
    for ( int d = 0; d < NDIM; d++ ) {
        b[2 * d + 0] = d_lb[d];
        b[2 * d + 1] = d_ub[d];
    }
    return b;
}


/********************************************************
 * Find the ideal point to split such that we divide     *
 *   both the space and points as much as possible       *
 ********************************************************/
template<uint8_t NDIM, class TYPE>
size_t kdtree2<NDIM, TYPE>::find_split( size_t N, const double *x )
{
    // Find the largest gap such that we also divide the points and space
    double lb = x[0];
    double ub = x[N - 1];
    int k     = 0;
    double q  = 0;
    for ( size_t i = 1; i < N; i++ ) {
        // Compute the quality of the split at the current location
        double q2 = ( x[i] - x[i - 1] ) * ( x[i] - lb ) * ( ub - x[i - 1] ) * i * ( N - i - 0 );
        if ( q2 > q ) {
            q = q2;
            k = i;
        }
    }
    return k;
}


/********************************************************
 * Nearest neighbor search                               *
 ********************************************************/
template<uint8_t NDIM, class TYPE>
constexpr double kdtree2<NDIM, TYPE>::norm( const Point &x, const Point &y )
{
    if constexpr ( NDIM == 1 ) {
        return ( x[0] - y[0] ) * ( x[0] - y[0] );
    } else if constexpr ( NDIM == 2 ) {
        return ( x[0] - y[0] ) * ( x[0] - y[0] ) + ( x[1] - y[1] ) * ( x[1] - y[1] );
    } else if constexpr ( NDIM == 3 ) {
        return ( x[0] - y[0] ) * ( x[0] - y[0] ) + ( x[1] - y[1] ) * ( x[1] - y[1] ) +
               ( x[2] - y[2] ) * ( x[2] - y[2] );
    } else {
        double dist = 0;
        for ( int d = 0; d < NDIM; d++ )
            dist += ( x[d] - y[d] ) * ( x[d] - y[d] );
        return dist;
    }
}
template<uint8_t NDIM, class TYPE>
std::tuple<std::array<double, NDIM>, TYPE>
kdtree2<NDIM, TYPE>::findNearest( const kdtree2::Point &x ) const
{
    std::tuple<Point, TYPE> nearest;
    // First, find dive into the structure to find where the position would be stored
    if ( d_left != nullptr ) {
        // Drill down the tree to find the node that should contain the point
        // As we travel back check the neighboring trees for any points that might be closer
        if ( x[d_split_dim] <= d_split ) {
            nearest = d_left->findNearest( x );
            d_right->checkNearest( x, nearest );
        } else {
            nearest = d_right->findNearest( x );
            d_left->checkNearest( x, nearest );
        }
    } else {
        // We are at the final node, find the closest value using the naive approach
        size_t k     = 0;
        double dist1 = norm( x, d_data->x[0] );
        for ( size_t i = 1; i < d_N; i++ ) {
            double dist2 = norm( x, d_data->x[i] );
            if ( dist2 < dist1 ) {
                k     = i;
                dist1 = dist2;
            }
        }
        nearest = std::tie( d_data->x[k], d_data->data[k] );
    }
    return nearest;
}
template<uint8_t NDIM, class TYPE>
void kdtree2<NDIM, TYPE>::checkNearest( const kdtree2::Point &x,
                                        std::tuple<kdtree2::Point, TYPE> &nearest ) const
{
    // Check if the point (and its radius) intersects with the current box
    double dist1 = norm( x, std::get<0>( nearest ) );
    double dist2 = 0.0;
    for ( int k = 0; k < NDIM; k++ ) {
        double d = std::max( d_lb[k] - x[k], x[k] - d_ub[k] );
        if ( d > 0.0 )
            dist2 += d * d;
    }
    if ( dist2 > dist1 )
        return;
    // Recursively search the subtrees
    if ( d_left != nullptr ) {
        d_left->checkNearest( x, nearest );
        d_right->checkNearest( x, nearest );
        return;
    }
    // We are at a base node, check the points for any that might be closer
    int64_t k = -1;
    for ( size_t i = 0; i < d_N; i++ ) {
        dist2 = norm( x, d_data->x[i] );
        if ( dist2 < dist1 ) {
            k     = i;
            dist1 = dist2;
        }
    }
    if ( k != -1 )
        nearest = std::tie( d_data->x[k], d_data->data[k] );
}


} // namespace AMP

#endif
