#include "AMP/ampmesh/Geometry.h"
#include "AMP/ampmesh/LogicalGeometry.h"
#include "AMP/ampmesh/MeshUtilities.h"
#include "AMP/ampmesh/MultiGeometry.h"
#include "AMP/utils/UnitTest.h"

#include "ProfilerApp.h"

#include <algorithm>
#include <random>


namespace AMP {
namespace Geometry {


// Generate a random direction
static inline Point genRandDir( int ndim )
{
    static std::random_device rd;
    static std::mt19937 gen( rd() );
    static std::uniform_real_distribution<double> dis( -1, 1 );
    return normalize( Point( ndim, { dis( gen ), dis( gen ), dis( gen ) } ) );
}


// Run all geometry based tests
void testGeometry( const AMP::Geometry::Geometry &geom, AMP::UnitTest &ut )
{
    // Get the physical dimension
    int ndim  = geom.getDim();
    auto name = geom.getName();
    PROFILE_START( "testGeometry " + name );
    // First get the centroid and the range
    auto center = geom.centroid();
    auto box    = geom.box();
    bool pass   = center.ndim() == ndim;
    for ( int d = 0; d < ndim; d++ )
        pass = pass && center[d] >= box.first[d] && center[d] <= box.second[d];
    if ( !pass )
        ut.failure( "testGeometry centroid/box: " + name );
    // Use a series of rays projecting from the centroid to get points on the surface
    PROFILE_START( "testGeometry-surface " + name );
    size_t N = 10000;
    std::vector<Point> surfacePoints;
    surfacePoints.reserve( N );
    bool all_hit = true;
    while ( surfacePoints.size() < N ) {
        // Keep each surface the ray hits
        auto dir = genRandDir( ndim );
        double d = geom.distance( center, dir );
        AMP_ASSERT( d == d );
        all_hit   = all_hit && d <= 0;
        double d2 = 0;
        while ( fabs( d ) < 1e100 ) {
            d2 += fabs( d );
            surfacePoints.push_back( center + d2 * dir );
            d2 += 1e-6;
            d = geom.distance( center + d2 * dir, dir );
            if ( surfacePoints.size() > N )
                break;
        }
    }
    pass = pass && !surfacePoints.empty();
    if ( surfacePoints.empty() )
        ut.failure( "testGeometry unable to get surface: " + name );
    if ( geom.inside( center ) && !all_hit )
        ut.failure( "testGeometry failed all rays hit: " + name );
    PROFILE_STOP( "testGeometry-surface " + name );
    // Verify each surface point is "inside" the object
    PROFILE_START( "testGeometry-inside " + name );
    bool pass_inside = true;
    for ( const auto &tmp : surfacePoints ) {
        bool inside = geom.inside( tmp );
        if ( !inside ) {
            pass_inside = false;
            std::cout << "testGeometry-inside: " << tmp << std::endl;
        }
    }
    pass = pass && pass_inside;
    if ( !pass_inside )
        ut.failure( "testGeometry surface inside geometry: " + name );
    PROFILE_STOP( "testGeometry-inside " + name );
    // Project each surface point in a random direction and back propagate to get the same point
    PROFILE_START( "testGeometry-distance " + name );
    bool pass_projection = true;
    auto length          = box.second - box.first;
    const double d0      = 0.2 * std::max( { length.x(), length.y(), length.z() } );
    for ( const auto &tmp : surfacePoints ) {
        auto ang = genRandDir( ndim );
        auto pos = tmp - d0 * ang;
        double d = fabs( geom.distance( pos, ang ) );
        for ( int it = 0; it < 1000 && d < d0 - 1e-5; it++ ) {
            // We may have crossed multiple surfaces, find the original
            d += 1e-6;
            auto pos2 = pos + d * ang;
            d += fabs( geom.distance( pos2, ang ) );
        }
        if ( fabs( d - d0 ) > 1e-5 ) {
            printf( "testGeometry-distance: %f %f\n", d0, d );
            pass_projection = false;
        }
    }
    pass = pass && pass_projection;
    if ( !pass_projection )
        ut.failure( "testGeometry distances do not match: " + name );
    PROFILE_STOP( "testGeometry-distance " + name );
    // Get a set of interior points by randomly sampling the space
    PROFILE_START( "testGeometry-sample " + name );
    static std::random_device rd;
    static std::mt19937 gen( rd() );
    std::uniform_real_distribution<double> dist[3];
    for ( int d = 0; d < ndim; d++ )
        dist[d] = std::uniform_real_distribution<double>( box.first[d], box.second[d] );
    std::vector<Point> interiorPoints;
    for ( int i = 0; i < 10000; i++ ) {
        Point pos( ndim, { 0, 0, 0 } );
        for ( int d = 0; d < ndim; d++ )
            pos[d] = dist[d]( gen );
        if ( geom.inside( pos ) )
            interiorPoints.push_back( pos );
    }
    PROFILE_STOP( "testGeometry-sample " + name );
    // Test logical transformation (if valid)
    PROFILE_START( "testGeometry-logical " + name );
    auto geom2 = dynamic_cast<const AMP::Geometry::LogicalGeometry *>( &geom );
    if ( geom2 ) {
        bool pass2 = true;
        for ( const auto &p : surfacePoints ) {
            auto p2 = geom2->logical( p );
            auto p3 = geom2->physical( p2 );
            for ( int d = 0; d < ndim; d++ )
                pass2 = pass2 && fabs( p[d] - p3[d] ) < 1e-6;
        }
        for ( const auto &p : interiorPoints ) {
            auto p2 = geom2->logical( p );
            auto p3 = geom2->physical( p2 );
            for ( int d = 0; d < ndim; d++ )
                pass2 = pass2 && fabs( p[d] - p3[d] ) < 1e-6;
        }
        pass = pass && pass2;
        if ( !pass_inside )
            ut.failure( "testGeometry physical-logical-physical: " + geom2->getName() );
    }
    PROFILE_STOP( "testGeometry-logical " + name );
    // Test getting surface normals
    if ( geom2 ) {
        bool passNorm = true;
        for ( const auto &p : surfacePoints ) {
            auto norm = geom2->surfaceNorm( p );
            double n  = sqrt( norm.x() * norm.x() + norm.y() * norm.y() + norm.z() * norm.z() );
            // auto p1   = p - 1e-5 * norm;
            // auto p2   = p + 1e-5 * norm;
            passNorm = passNorm && fabs( n - 1.0 ) < 1e-6;
            // passNorm  = passNorm && geom.inside( p1 ) && !geom.inside( p2 );
        }
        pass = pass && passNorm;
        if ( !passNorm )
            ut.failure( "testGeometry surfaceNorm: " + geom2->getName() );
    }
    // Test getting the volume
    {
        PROFILE_START( "testGeometry-volume " + name );
        double volume    = geom.volume();
        double boxVolume = 1.0;
        for ( int d = 0; d < ndim; d++ )
            boxVolume *= box.second[d] - box.first[d];
        bool passVol = volume > 0;
        if ( ndim == static_cast<int>( geom.getGeomType() ) )
            passVol = passVol && volume <= boxVolume;
        pass = pass && passVol;
        if ( !passVol )
            ut.failure( "testGeometry volume: " + name );
        // Test mesh utilities volume overlap
        auto multigeom = dynamic_cast<const MultiGeometry *>( &geom );
        if ( ndim == static_cast<int>( geom.getGeomType() ) && multigeom == nullptr ) {
            auto tmp      = AMP::Mesh::volumeOverlap( geom, std::vector<int>( ndim, 35 ) );
            double vol2   = tmp.sum();
            bool passVol2 = fabs( vol2 - volume ) < 0.01 * volume;
            pass          = pass && passVol2;
            if ( !passVol2 )
                ut.failure( "testGeometry volumeOverlap: " + name );
        }
        PROFILE_STOP( "testGeometry-volume " + name );
    }
    // Finished with all tests
    if ( pass )
        ut.passes( "testGeometry: " + name );
    PROFILE_STOP( "testGeometry " + name );
}


} // namespace Geometry
} // namespace AMP
