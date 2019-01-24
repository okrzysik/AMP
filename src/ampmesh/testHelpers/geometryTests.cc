#include "AMP/ampmesh/Geometry.h"
#include "AMP/utils/UnitTest.h"

#include <random>


namespace AMP {
namespace Geometry {


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
    std::random_device rd;
    std::mt19937 gen( rd() );
    std::uniform_real_distribution<double> dis( -1, 1 );
    bool all_hit = true;
    for ( int i = 0; i < 10000; i++ ) {
        Point dir   = { dis( gen ), dis( gen ), dis( gen ) };
        double norm = sqrt( dir.x() * dir.x() + dir.y() * dir.y() + dir.z() * dir.z() );
        dir         = { dir.x() / norm, dir.y() / norm, dir.z() / norm };
        double dist = geom.distance( center, dir );
        AMP_ASSERT( dist == dist );
        all_hit = all_hit && dist <= 0;
        if ( dist < 1e100 )
            surfacePoints.push_back( center + dir * std::abs( dist ) );
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
    // Project each surface point beyond the object and back propagate to get the same point
    bool pass_projection = true;
    for ( const auto &tmp : surfacePoints ) {
        double d0       = 1e-2;
        auto ang        = normalize( center - tmp );
        auto pos        = tmp - d0 * ang;
        double d        = geom.distance( pos, ang );
        pass_projection = pass_inside && fabs( fabs( d ) - d0 ) < 1e-5;
    }
    pass = pass && pass_projection;
    if ( !pass_projection )
        ut.failure( "testGeometry distances do not match: " + geom.getName() );
    // Get a set of interior points by randomly sampling the space
    auto range = geom.box();
    std::uniform_real_distribution<double> dis_x( range.first.x(), range.second.x() );
    std::uniform_real_distribution<double> dis_y( range.first.y(), range.second.y() );
    std::uniform_real_distribution<double> dis_z( range.first.z(), range.second.z() );
    std::vector<Point> interiorPoints;
    for ( int i = 0; i < 10000; i++ ) {
        Point pos( ndim, { dis_x( gen ), dis_y( gen ), dis_z( gen ) } );
        if ( geom.inside( pos ) )
            interiorPoints.push_back( pos );
    }
    // Test logical transformation (if valid)
    if ( geom.isLogical() ) {
        bool pass2 = true;
        for ( const auto &p : surfacePoints ) {
            auto p2 = geom.logical( p );
            auto p3 = geom.physical( p2 );
            for ( int d = 0; d < ndim; d++ )
                pass2 = pass2 && fabs( p[d] - p3[d] ) < 1e-6;
        }
        for ( const auto &p : interiorPoints ) {
            auto p2 = geom.logical( p );
            auto p3 = geom.physical( p2 );
            for ( int d = 0; d < ndim; d++ )
                pass2 = pass2 && fabs( p[d] - p3[d] ) < 1e-6;
        }
        pass = pass && pass2;
        if ( !pass_inside )
            ut.failure( "testGeometry physical-logical-physical: " + geom.getName() );
    }
    if ( pass )
        ut.passes( "testGeometry: " + geom.getName() );
}


} // namespace Geometry
} // namespace AMP
