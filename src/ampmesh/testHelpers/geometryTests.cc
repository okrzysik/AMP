#include "AMP/ampmesh/Geometry.h"
#include "AMP/ampmesh/LogicalGeometry.h"
#include "AMP/ampmesh/MeshUtilities.h"
#include "AMP/ampmesh/MultiGeometry.h"
#include "AMP/utils/UnitTest.h"

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
    int ndim = geom.getDim();
    // First get the centroid and the range
    auto center = geom.centroid();
    auto box    = geom.box();
    bool pass   = center.ndim() == ndim;
    for ( int d = 0; d < ndim; d++ )
        pass = pass && center[d] >= box.first[d] && center[d] <= box.second[d];
    if ( !pass )
        ut.failure( "testGeometry centroid/box: " + geom.getName() );
    // Use a series of rays projecting from the centroid to get points on the surface
    std::vector<Point> surfacePoints;
    bool all_hit = true;
    for ( int i = 0; i < 10000; i++ ) {
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
        }
    }
    pass = pass && !surfacePoints.empty();
    if ( surfacePoints.empty() )
        ut.failure( "testGeometry unable to get surface: " + geom.getName() );
    if ( geom.inside( center ) && !all_hit )
        ut.failure( "testGeometry failed all rays hit: " + geom.getName() );
    // Verify each surface point is "inside" the object
    bool pass_inside = true;
    for ( const auto &tmp : surfacePoints )
        pass_inside = pass_inside && geom.inside( tmp );
    pass = pass && pass_inside;
    if ( !pass_inside )
        ut.failure( "testGeometry surface inside geometry: " + geom.getName() );
    // Project each surface point in a random direction and back propagate to get the same point
    bool pass_projection = true;
    auto length          = box.second - box.first;
    const double d0      = 0.2 * std::max( { length.x(), length.y(), length.z() } );
    for ( const auto &tmp : surfacePoints ) {
        auto ang = genRandDir( ndim );
        auto pos = tmp - d0 * ang;
        double d = fabs( geom.distance( pos, ang ) );
        while ( d < d0 - 1e-5 ) {
            // We may have crossed multiple surfaces, find the original
            d += 1e-6;
            auto pos2 = pos + d * ang;
            d += fabs( geom.distance( pos2, ang ) );
        }
        if ( fabs( d - d0 ) > 1e-5 )
            pass_projection = false;
    }
    pass = pass && pass_projection;
    if ( !pass_projection )
        ut.failure( "testGeometry distances do not match: " + geom.getName() );
    // Get a set of interior points by randomly sampling the space
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
    // Test logical transformation (if valid)
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
        double volume    = geom.volume();
        double boxVolume = 1.0;
        for ( int d = 0; d < ndim; d++ )
            boxVolume *= box.second[d] - box.first[d];
        bool passVol = volume > 0;
        if ( ndim == static_cast<int>( geom.getGeomType() ) )
            passVol = passVol && volume <= boxVolume;
        pass = pass && passVol;
        if ( !passVol )
            ut.failure( "testGeometry volume: " + geom.getName() );
        // Test mesh utilities volume overlap
        auto multigeom = dynamic_cast<const MultiGeometry *>( &geom );
        if ( ndim == static_cast<int>( geom.getGeomType() ) && multigeom == nullptr ) {
            auto tmp      = AMP::Mesh::volumeOverlap( geom, { 35, 35, 35 } );
            double vol2   = tmp.sum();
            bool passVol2 = fabs( vol2 - volume ) < 0.01 * volume;
            pass          = pass && passVol2;
            if ( !passVol2 )
                ut.failure( "testGeometry volumeOverlap: " + geom.getName() );
        }
    }
    // Finished with all tests
    if ( pass )
        ut.passes( "testGeometry: " + geom.getName() );
}


} // namespace Geometry
} // namespace AMP
