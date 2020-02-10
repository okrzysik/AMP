#include "AMP/utils/kdtree.h"
#include "AMP/utils/Utilities.h"

#if USE_AMP_MESH
#include "AMP/ampmesh/Mesh.h"
#endif

#include "ProfilerApp.h"

#include <cmath>
#include <iostream>
#include <limits>

#define ERROR_MSG AMP_ERROR


namespace AMP {


/********************************************************
 * Constructor                                           *
 ********************************************************/
kdtree::kdtree( const int N_dim, const size_t N, const double *const *x ) : kdtree()
{
    initialize( N_dim, N, x );
}
void kdtree::initialize( const int N_dim, const size_t N, const double *const *x )
{
    if ( N_dim < 1 || N_dim > 64 )
        ERROR_MSG( "Invalid dimension" );
    if ( N >
         (size_t)
             std::numeric_limits<int>::max() ) // Internal data structure still uses int for index
        ERROR_MSG( "Maximum number of points exceeeded" );
    PROFILE_START( "constructor" );
    // Create the head node
    d_dim            = N_dim;
    d_N              = N;
    d_tree.N_dim     = d_dim;
    d_tree.split_dim = 0;
    d_tree.x_start   = new double[d_dim];
    d_tree.x_end     = new double[d_dim];
    for ( unsigned char k = 0; k < d_dim; k++ ) {
        d_tree.x_start[k] = 1e100;
        d_tree.x_end[k]   = -1e100;
        for ( size_t i = 0; i < N; i++ ) {
            if ( x[k][i] < d_tree.x_start[k] )
                d_tree.x_start[k] = x[k][i];
            if ( x[k][i] > d_tree.x_end[k] )
                d_tree.x_end[k] = x[k][i];
        }
    }
    // Add the points to the current level and sort them
    d_tree.N     = static_cast<unsigned int>( N );
    d_tree.index = new int[N];
    d_tree.x     = new double[N * N_dim];
    auto arr     = new double[N];
    auto *brr    = new int64_t[N];
    for ( size_t i = 0; i < N; i++ ) {
        arr[i] = x[0][i];
        brr[i] = i;
    }
    AMP::Utilities::quicksort( N, arr, brr );
    for ( size_t i = 0; i < N; i++ ) {
        auto j          = (int) brr[i];
        d_tree.index[i] = j;
        for ( int k = 0; k < N_dim; k++ )
            d_tree.x[N_dim * i + k] = x[k][j];
    }
    delete[] arr;
    delete[] brr;
    // Divide the current tree
    kdtree::split_tree( &d_tree );
    PROFILE_STOP( "constructor" );
}
#if USE_AMP_MESH
kdtree::kdtree( const std::vector<AMP::Mesh::MeshPoint<double>> &x ) : kdtree()
{
    if ( x.empty() )
        return;
    int ndim = x[0].ndim();
    double *x2[5];
    for ( int d = 0; d < ndim; d++ )
        x2[d] = new double[x.size()];
    for ( size_t i = 0; i < x.size(); i++ ) {
        for ( int d = 0; d < ndim; d++ )
            x2[d][i] = x[i][d];
    }
    initialize( ndim, x.size(), x2 );
    for ( int d = 0; d < ndim; d++ )
        delete[] x2[d];
}
#endif


/********************************************************
 * Destructor                                            *
 ********************************************************/
kdtree::~kdtree() = default;


/********************************************************
 * Specialized constructors                              *
 ********************************************************/
std::shared_ptr<kdtree> kdtree::create2d( const size_t N, const double *x, const double *y )
{
    const double *x_tmp[2];
    x_tmp[0] = x;
    x_tmp[1] = y;
    return std::make_shared<kdtree>( 2, N, x_tmp );
}


/********************************************************
 * Specialized constructor (3D)                          *
 ********************************************************/
std::shared_ptr<kdtree>
kdtree::create3d( const size_t N, const double *x, const double *y, const double *z )
{
    const double *x_tmp[3];
    x_tmp[0] = x;
    x_tmp[1] = y;
    x_tmp[2] = z;
    return std::make_shared<kdtree>( 3, N, x_tmp );
}


/********************************************************
 * Return the bounding box of the tree                   *
 ********************************************************/
std::vector<double> kdtree::box()
{
    std::vector<double> range( 2 * d_dim );
    for ( unsigned int d = 0; d < d_dim; d++ ) {
        range[2 * d + 0] = d_tree.x_start[d];
        range[2 * d + 1] = d_tree.x_end[d];
    }
    return range;
}


/********************************************************
 * Return the memory_usage                               *
 ********************************************************/
size_t kdtree::memory_usage() const
{
    return sizeof( kdtree ) - sizeof( kdtree_struct ) + d_tree.memory_usage();
}


/********************************************************
 * Recursively splt the tree                             *
 ********************************************************/
void kdtree::split_tree( kdtree_struct *tree )
{
    // Don't split the tree in 1D (we have a sorted array to use)
    if ( tree->N_dim == 1 )
        return;
    // Check if the tree needs to be split
    unsigned char N_dim      = tree->N_dim;
    unsigned char split_dim1 = tree->split_dim;
    unsigned char split_dim2 = ( tree->split_dim + 1 ) % tree->N_dim;
    tree->left               = nullptr;
    tree->right              = nullptr;
    if ( tree->N < 40 )
        return; // Brute force search is more efficient
    // Check if all the points have the same x[k] value
    double eps = 1e-12 *
                 ( std::fabs( tree->x_start[split_dim1] ) + std::abs( tree->x_end[split_dim1] ) ) /
                 2.0;
    if ( std::fabs( tree->x_end[split_dim1] - tree->x_start[split_dim1] ) < eps ) {
        // The points share the same value, we should try splitting along a different direction,
        // return for now
        return;
    }
    // Determine the number of points in each half
    tree->x_split = std::fabs( tree->x_end[split_dim1] + tree->x_start[split_dim1] ) / 2;
    for ( size_t i = 0; i < tree->N / 2; i++ ) {
        size_t j1 = ( tree->N / 2 ) - i - 1;
        size_t j2 = ( tree->N / 2 ) + i;
        if ( std::fabs( tree->x[N_dim * j1 + split_dim1] - tree->x[N_dim * j2 + split_dim1] ) >
             eps ) {
            tree->x_split =
                ( tree->x[N_dim * j1 + split_dim1] + tree->x[N_dim * j2 + split_dim1] ) / 2.0;
            break;
        }
    }
    size_t N1 = 0;
    size_t N2 = 0;
    for ( size_t i = 0; i < tree->N; i++ ) {
        if ( tree->x[N_dim * i + split_dim1] <= tree->x_split )
            N1++;
        else
            break;
    }
    N2 = tree->N - N1;
    if ( N1 == 0 || N2 == 0 ) {
        // This should not occur
        return;
    }
    // Create the left and right nodes
    tree->left             = new kdtree_struct;
    tree->right            = new kdtree_struct;
    tree->left->N_dim      = tree->N_dim;
    tree->right->N_dim     = tree->N_dim;
    tree->left->x          = tree->x;
    tree->right->x         = tree->x;
    tree->left->split_dim  = split_dim2;
    tree->right->split_dim = split_dim2;
    // Divide the points in half
    tree->left->N        = static_cast<unsigned int>( N1 );
    tree->right->N       = static_cast<unsigned int>( N2 );
    tree->left->index    = new int[N1];
    tree->right->index   = new int[N2];
    tree->left->x        = new double[N_dim * N1];
    tree->right->x       = new double[N_dim * N2];
    tree->left->x_start  = new double[N_dim];
    tree->right->x_start = new double[N_dim];
    tree->left->x_end    = new double[N_dim];
    tree->right->x_end   = new double[N_dim];
    // Sort the indices in the left-half of the tree
    auto arr  = new double[tree->left->N];
    auto *brr = new int64_t[tree->left->N];
    for ( size_t i = 0; i < tree->left->N; i++ ) {
        arr[i] = tree->x[N_dim * i + split_dim2];
        brr[i] = i;
    }
    AMP::Utilities::quicksort( tree->left->N, arr, brr );
    for ( size_t i = 0; i < tree->left->N; i++ ) {
        auto j               = (size_t) brr[i];
        tree->left->index[i] = tree->index[j];
        for ( size_t k = 0; k < N_dim; k++ )
            tree->left->x[N_dim * i + k] = tree->x[N_dim * j + k];
    }
    delete[] arr;
    delete[] brr;
    size_t m            = N1 / 2;
    tree->left->x_split = tree->left->x[N_dim * m + split_dim2];
    for ( size_t k = 0; k < N_dim; k++ ) {
        tree->left->x_start[k] = 1e100;
        tree->left->x_end[k]   = -1e100;
        for ( size_t i = 0; i < tree->left->N; i++ ) {
            if ( tree->left->x[N_dim * i + k] < tree->left->x_start[k] )
                tree->left->x_start[k] = tree->left->x[N_dim * i + k];
            if ( tree->left->x[N_dim * i + k] > tree->left->x_end[k] )
                tree->left->x_end[k] = tree->left->x[N_dim * i + k];
        }
    }
    // Sort the indices in the left-half of the tree
    arr = new double[tree->right->N];
    brr = new int64_t[tree->right->N];
    for ( size_t i = 0; i < tree->right->N; i++ ) {
        arr[i] = tree->x[N_dim * ( i + N1 ) + split_dim2];
        brr[i] = i + N1;
    }
    AMP::Utilities::quicksort( tree->right->N, arr, brr );
    for ( size_t i = 0; i < tree->right->N; i++ ) {
        auto j                = (size_t) brr[i];
        tree->right->index[i] = tree->index[j];
        for ( size_t k = 0; k < N_dim; k++ )
            tree->right->x[N_dim * i + k] = tree->x[N_dim * j + k];
    }
    delete[] arr;
    delete[] brr;
    m                    = N2 / 2;
    tree->right->x_split = tree->right->x[N_dim * m + split_dim2];
    for ( size_t k = 0; k < N_dim; k++ ) {
        tree->right->x_start[k] = 1e100;
        tree->right->x_end[k]   = -1e100;
        for ( size_t i = 0; i < tree->right->N; i++ ) {
            if ( tree->right->x[N_dim * i + k] < tree->right->x_start[k] )
                tree->right->x_start[k] = tree->right->x[N_dim * i + k];
            if ( tree->right->x[N_dim * i + k] > tree->right->x_end[k] )
                tree->right->x_end[k] = tree->right->x[N_dim * i + k];
        }
    }
    // We no longer need local copies of the points, free the memory
    delete[] tree->x;
    delete[] tree->index;
    tree->x     = nullptr;
    tree->index = nullptr;
    // Further split each half of the tree
    kdtree::split_tree( tree->left );
    kdtree::split_tree( tree->right );
}


/********************************************************
 * Nearest neighbor search                               *
 ********************************************************/
size_t kdtree::find_nearest( const double *x, double *dist, double *pos ) const
{
    PROFILE_START( "find_nearest single", 5 );
    double dist2 = 1e100;
    double pos2[64];
    size_t index = find_nearest_tree( &d_tree, x, dist2, pos2 );
    if ( dist != nullptr )
        *dist = dist2;
    if ( pos != nullptr ) {
        for ( unsigned int d = 0; d < d_dim; d++ )
            pos[d] = pos2[d];
    }
    PROFILE_STOP( "find_nearest single", 5 );
    return index;
}
void kdtree::find_nearest( int N, const double *x, size_t *index, double *dist, double *pos ) const
{
    if ( N == 0 )
        return;
    if ( N < 0 )
        ERROR_MSG( "N must be >= 0" );
    if ( index == nullptr )
        ERROR_MSG( "index may not be null" );
    PROFILE_START( "find_nearest multiple", 3 );
    for ( int i = 0; i < N; i++ ) {
        double dist2 = 1e100;
        double pos2[64];
        index[i] = find_nearest_tree( &d_tree, &x[d_dim * i], dist2, pos2 );
        if ( dist != nullptr )
            dist[i] = dist2;
        if ( pos != nullptr ) {
            for ( unsigned int d = 0; d < d_dim; d++ )
                pos[d_dim * i + d] = pos2[d];
        }
    }
    PROFILE_STOP( "find_nearest multiple", 3 );
}
size_t kdtree::find_nearest2d( const double x, const double y ) const
{
    PROFILE_START( "find_nearest 2d", 5 );
    double xy[2];
    xy[0] = x;
    xy[1] = y;
    double dist2;
    double pos2[2];
    size_t index = find_nearest_tree( &d_tree, xy, dist2, pos2 );
    PROFILE_STOP( "find_nearest 2d", 5 );
    return index;
}

size_t kdtree::find_nearest3d( const double x, const double y, const double z ) const
{
    PROFILE_START( "find_nearest 3d", 5 );
    double xy[3];
    xy[0] = x;
    xy[1] = y;
    xy[2] = z;
    double dist2;
    double pos2[3];
    size_t index = find_nearest_tree( &d_tree, xy, dist2, pos2 );
    PROFILE_STOP( "find_nearest 3d", 5 );
    return index;
}

size_t
kdtree::find_nearest_tree( const kdtree_struct *tree, const double *x, double &dist, double *pos )
{
    unsigned char N_dim = tree->N_dim;
    size_t index        = tree->N;
    if ( N_dim == 1 && tree->x != nullptr ) {
        // Perform a hash search on 1D arrays
        index = AMP::Utilities::findfirst( tree->N, tree->x, x[0] );
        if ( index == tree->N )
            index = tree->N - 1;
        dist = std::fabs( x[0] - tree->x[index] );
        *pos = tree->x[index];
        return tree->index[index];
    }
    // First, find dive into the structure to find where the position would be stored
    if ( tree->left != nullptr ) {
        // Drill down the tree to find the node that should contain the point
        if ( x[tree->split_dim] <= tree->x_split )
            index = find_nearest_tree( tree->left, x, dist, pos );
        else
            index = find_nearest_tree( tree->right, x, dist, pos );
        // As we travel back check the neighboring trees for any points that might be closer
        if ( x[tree->split_dim] <= tree->x_split )
            check_neighbor( tree->right, x, &index, dist, pos );
        else
            check_neighbor( tree->left, x, &index, dist, pos );
    } else {
        // We are at the final node, find the closest value using the naive approach
        size_t j        = tree->N;
        double min_dist = 1e100;
        double dist2;
        if ( N_dim == 2 ) {
            for ( size_t i = 0; i < tree->N; i++ ) {
                dist2 = ( tree->x[2 * i + 0] - x[0] ) * ( tree->x[2 * i + 0] - x[0] );
                dist2 += ( tree->x[2 * i + 1] - x[1] ) * ( tree->x[2 * i + 1] - x[1] );
                if ( dist2 < min_dist ) {
                    j        = i;
                    min_dist = dist2;
                }
            }
        } else {
            for ( size_t i = 0; i < tree->N; i++ ) {
                dist2 = 0.0;
                for ( unsigned char k = 0; k < tree->N_dim; k++ )
                    dist2 += ( tree->x[N_dim * i + k] - x[k] ) * ( tree->x[N_dim * i + k] - x[k] );
                if ( dist2 < min_dist ) {
                    j        = i;
                    min_dist = dist2;
                }
            }
        }
        index = tree->index[j];
        for ( unsigned char k = 0; k < N_dim; k++ )
            pos[k] = tree->x[N_dim * j + k];
        dist = sqrt( min_dist );
    }
    return index;
}
void kdtree::check_neighbor(
    const kdtree_struct *tree, const double *x, size_t *index, double &dist, double *pos )
{
    // Check if the point (and its radius) intersects with the current box
    unsigned char N_dim = tree->N_dim;
    double dist1        = dist * dist;
    double dist2        = 0.0;
    for ( unsigned char k = 0; k < N_dim; k++ ) {
        double d = std::max( tree->x_start[k] - x[k], x[k] - tree->x_end[k] );
        if ( d > 0.0 )
            dist2 += d * d;
    }
    if ( dist2 > dist1 )
        return;
    // Recursively search the subtrees
    if ( tree->left != nullptr ) {
        check_neighbor( tree->left, x, index, dist, pos );
        check_neighbor( tree->right, x, index, dist, pos );
        return;
    }
    // We are at a base node, check the points for any that might be closer
    if ( N_dim == 2 ) {
        for ( size_t i = 0; i < tree->N; i++ ) {
            dist2 = ( tree->x[2 * i + 0] - x[0] ) * ( tree->x[2 * i + 0] - x[0] );
            dist2 += ( tree->x[2 * i + 1] - x[1] ) * ( tree->x[2 * i + 1] - x[1] );
            if ( dist2 < dist1 ) {
                *index = tree->index[i];
                dist   = sqrt( dist2 );
                pos[0] = tree->x[2 * i + 0];
                pos[1] = tree->x[2 * i + 1];
                dist1  = dist2;
            }
        }
    } else {
        for ( size_t i = 0; i < tree->N; i++ ) {
            dist2 = 0.0;
            for ( unsigned char k = 0; k < tree->N_dim; k++ )
                dist2 += ( tree->x[N_dim * i + k] - x[k] ) * ( tree->x[N_dim * i + k] - x[k] );
            if ( dist2 < dist1 ) {
                *index = tree->index[i];
                dist   = sqrt( dist2 );
                for ( unsigned char k = 0; k < N_dim; k++ )
                    pos[k] = tree->x[N_dim * i + k];
                dist1 = dist2;
            }
        }
    }
}


/********************************************************
 * kdtree_struct functions                               *
 ********************************************************/
kdtree::kdtree_struct::kdtree_struct()
{
    N_dim     = 0;
    split_dim = -1;
    N         = 0;
    x_split   = 0.0;
    x_start   = nullptr;
    x_end     = nullptr;
    x         = nullptr;
    index     = nullptr;
    left      = nullptr;
    right     = nullptr;
}
kdtree::kdtree_struct::~kdtree_struct()
{
    delete left;
    delete right;
    delete[] x_start;
    delete[] x_end;
    delete[] index;
    delete[] x;
}
size_t kdtree::kdtree_struct::memory_usage() const
{
    size_t bytes = sizeof( kdtree_struct );
    if ( x_start != nullptr )
        bytes += N_dim * sizeof( double );
    if ( x_end != nullptr )
        bytes += N_dim * sizeof( double );
    if ( x != nullptr )
        bytes += N * N_dim * sizeof( int );
    if ( index != nullptr )
        bytes += N * sizeof( double );
    if ( left != nullptr )
        bytes += left->memory_usage();
    if ( right != nullptr )
        bytes += right->memory_usage();
    return bytes;
}


} // namespace AMP
