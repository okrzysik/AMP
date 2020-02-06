#include <algorithm>
#include <stdio.h>
#include <string>
#include <vector>

#include "ProfilerApp.h"


//#define R_INT 0x20000000
#define R_INT 1000000


// Compute the L2 norm of a vector ||x||
inline double L2norm( size_t N, const double *x, const bool *mask = nullptr )
{
    double norm = 0;
    if ( mask == nullptr ) {
        for ( size_t i = 0; i < N; i++ )
            norm += x[i] * x[i];
    } else {
        for ( size_t i = 0; i < N; i++ )
            norm += mask[i] ? x[i] * x[i] : 0;
    }
    return sqrt( norm );
}


// Compute the L2 error norm of a vector ||x1-x2||
inline double L2errNorm( size_t N, const double *x1, const double *x2, const bool *mask = nullptr )
{
    double norm = 0;
    if ( mask == nullptr ) {
        for ( size_t i = 0; i < N; i++ )
            norm += ( x1[i] - x2[i] ) * ( x1[i] - x2[i] );
    } else {
        for ( size_t i = 0; i < N; i++ )
            norm += mask[i] ? ( x1[i] - x2[i] ) * ( x1[i] - x2[i] ) : 0;
    }
    return sqrt( norm );
}


// Random number generator in the interval [0,1]
inline double rand_double()
{
    const double rmax = static_cast<double>( RAND_MAX );
    double x          = static_cast<double>( rand() ) / rmax;
    x += 1e-4 * static_cast<double>( rand() ) / rmax;
    x += 1e-8 * static_cast<double>( rand() ) / rmax;
    x += 1e-12 * static_cast<double>( rand() ) / rmax;
    x /= 1.0 + 1e-4 + 1e-8 + 1e-12;
    return x;
}


// Random number generator in the interval [0,R_INT-1]
// Note: this is not a perfect distribution but should be good enough for our purposes
inline int rand_int()
{
    unsigned int i1 =
        ( static_cast<unsigned int>( rand() ) * 0x9E3779B9 ) % 32767; // 2^32*0.5*(sqrt(5)-1)
    unsigned int i2 =
        ( static_cast<unsigned int>( rand() ) * 0x9E3779B9 ) % 32767; // 2^32*0.5*(sqrt(5)-1)
    return static_cast<int>( ( i1 + ( i2 << 15 ) ) % R_INT );
}


// Structure to store a point
template<int NDIM>
struct PointInt {
    int x[NDIM];
    PointInt() { memset( x, 0, sizeof( x ) ); }
    inline uint64_t R2() const
    {
        uint64_t R2 = 0;
        for ( auto &elem : x ) {
            uint64_t x2 = abs( elem );
            R2 += x2 * x2;
        }
        return R2;
    }
    inline bool operator==( const PointInt &rhs ) const
    {
        bool equal = true;
        for ( int d = 0; d < NDIM; d++ )
            equal = equal && x[d] == rhs.x[d];
        return equal;
    }
    inline bool operator!=( const PointInt &rhs ) const { return !operator==( rhs ); }
    inline bool operator>=( const PointInt &rhs ) const
    {
        for ( int d = 0; d < NDIM; d++ ) {
            if ( x[d] > rhs.x[d] )
                return true;
            if ( rhs.x[d] > x[d] )
                return false;
        }
        return true;
    }
    inline bool operator<=( const PointInt &rhs ) const
    {
        for ( int d = 0; d < NDIM; d++ ) {
            if ( x[d] < rhs.x[d] )
                return true;
            if ( rhs.x[d] < x[d] )
                return false;
        }
        return true;
    }
    inline bool operator>( const PointInt &rhs ) const
    {
        return operator>=( rhs ) && !operator==( rhs );
    }
    inline bool operator<( const PointInt &rhs ) const
    {
        return operator<=( rhs ) && !operator==( rhs );
    }
};


// Create a set of random points in a Nd-hypersphere
// Note: we use a combination of a fixed grid and random points to control the error
template<class TYPE>
std::vector<TYPE> createRandomPoints( int ndim, int N );
template<int NDIM>
std::vector<PointInt<NDIM>> createRandomPointsInt( int N )
{
    std::vector<PointInt<NDIM>> points;
    if ( N == NDIM + 1 ) {
        points.resize( N );
        for ( int d = 0; d < NDIM; d++ )
            points[d + 1].x[d] = R_INT;
        return points;
    }
    PROFILE_START( "createRandomPoints", 1 );
    const uint64_t R_INT2 = static_cast<uint64_t>( R_INT ) * static_cast<uint64_t>( R_INT );
    points.reserve( N + 10 );
    // Create a logical Nd-hypercube on [-1,1] and keep only the points within R<=R_INT
    if ( N > 10 ) {
        int Nd          = static_cast<int>( floor( pow( N, 1.0 / NDIM ) ) );
        Nd              = std::min( N / 2, Nd );
        Nd              = 2 * ( Nd / 2 ) + 1;
        const double dx = 2.0 / static_cast<double>( Nd - 1 );
        for ( size_t k = 0; k < pow( Nd, NDIM ); k++ ) {
            PointInt<NDIM> p;
            size_t j = k;
            for ( int d = 0; d < NDIM; d++ ) {
                double pos = ( j % Nd ) * dx - 1.0;
                p.x[d]     = static_cast<int>( round( static_cast<double>( R_INT ) * pos ) );
                j /= Nd;
            }
            if ( p.R2() <= R_INT2 )
                points.push_back( p );
        }
    }
    // Add random points
    while ( static_cast<int>( points.size() ) < N + 5 ) {
        PointInt<NDIM> p;
        for ( int d = 0; d < NDIM; d++ )
            p.x[d] = ( rand() % 2 == 1 ? 1 : -1 ) * rand_int();
        if ( p.R2() <= R_INT2 )
            points.push_back( p );
    }
    // Check for duplicates and remove them
    std::sort( points.begin(), points.end() );
    for ( size_t i = points.size() - 1; i > 0; i-- ) {
        if ( points[i] == points[i - 1] ) {
            points[i] = points.back();
            points.resize( points.size() - 1 );
        }
    }
    points.resize( std::min<size_t>( points.size(), N ) );
    // Resort the points in random order
    std::random_shuffle( points.begin(), points.end() );
    PROFILE_STOP( "createRandomPoints", 1 );
    return points;
}
std::vector<int> getPointListInt( int ndim, int N )
{
    std::vector<int> points;
    if ( ndim == 1 ) {
        std::vector<PointInt<1>> points2 = createRandomPointsInt<1>( N );
        points.resize( ndim * points2.size(), 0 );
        for ( size_t i = 0; i < points2.size(); i++ )
            points[i] = points2[i].x[0];
    } else if ( ndim == 2 ) {
        std::vector<PointInt<2>> points2 = createRandomPointsInt<2>( N );
        points.resize( ndim * points2.size(), 0 );
        for ( size_t i = 0; i < points2.size(); i++ ) {
            points[2 * i + 0] = points2[i].x[0];
            points[2 * i + 1] = points2[i].x[1];
        }
    } else if ( ndim == 3 ) {
        std::vector<PointInt<3>> points2 = createRandomPointsInt<3>( N );
        points.resize( ndim * points2.size(), 0 );
        for ( size_t i = 0; i < points2.size(); i++ ) {
            points[3 * i + 0] = points2[i].x[0];
            points[3 * i + 1] = points2[i].x[1];
            points[3 * i + 2] = points2[i].x[2];
        }
    } else if ( ndim == 4 ) {
        std::vector<PointInt<4>> points2 = createRandomPointsInt<4>( N );
        points.resize( ndim * points2.size(), 0 );
        for ( size_t i = 0; i < points2.size(); i++ ) {
            for ( int d = 0; d < ndim; d++ )
                points[i * ndim + d] = points2[i].x[d];
        }
    } else if ( ndim == 5 ) {
        std::vector<PointInt<5>> points2 = createRandomPointsInt<5>( N );
        points.resize( ndim * points2.size(), 0 );
        for ( size_t i = 0; i < points2.size(); i++ ) {
            for ( int d = 0; d < ndim; d++ )
                points[i * ndim + d] = points2[i].x[d];
        }
    }
    return points;
}
template<>
std::vector<int> createRandomPoints<int>( int ndim, int N )
{
    return getPointListInt( ndim, N );
}
/*template<>
std::vector<double> createRandomPoints<double>( int ndim, int N )
{
    std::vector<int> points = getPointListInt( ndim, N );
    std::vector<double> points2(points.size(),0);
    const double R_intd = R_INT;
    const double rand_max = RAND_MAX;
    for (size_t i=0; i<points.size(); i++) {
        points2[i] = static_cast<double>(points[i])/R_intd;
        points2[i] += 1.0/R_intd*((rand()/rand_max)-0.5);     // Wiggle the point slightly
    }
    return points2;
}*/
template<>
std::vector<double> createRandomPoints<double>( int ndim, int N )
{
    PROFILE_START( "createRandomPointsDouble", 1 );
    std::vector<double> points( N * ndim, 0.0 );
    int i = 0;
    // Create a Nd-hypercube on [-1,1] and keep only the points within R<=1
    int Nd = static_cast<int>( floor( pow( N, 1.0 / ndim ) ) );
    Nd     = std::min( N / 2, Nd );
    for ( size_t k = 0; k < pow( Nd, ndim ); k++ ) {
        double x[10] = { 0 };
        double R     = 0.0;
        size_t j     = k;
        for ( int d = 0; d < ndim; d++ ) {
            x[d] = 2.0 * ( ( j % Nd ) + 0.5 ) / static_cast<double>( Nd ) - 1.0;
            R += x[d] * x[d];
            j /= Nd;
        }
        if ( R <= 1.0 ) {
            for ( int d = 0; d < ndim; d++ )
                points[d + i * ndim] = x[d];
            i++;
        }
    }
    // Add random points
    while ( i < N ) {
        double x[10] = { 0 };
        double R     = 0.0;
        for ( int d = 0; d < ndim; d++ ) {
            x[d] = 2.0 * rand_double() - 1.0;
            R += x[d] * x[d];
        }
        if ( R <= 1.0 ) {
            for ( int d = 0; d < ndim; d++ )
                points[d + i * ndim] = x[d];
            i++;
        }
    }
    PROFILE_STOP( "createRandomPointsDouble", 1 );
    return points;
}
