#include "AMP/geometry/MultiGeometry.h"
#include "AMP/IO/RestartManager.h"
#include "AMP/geometry/Geometry.h"
#include "AMP/utils/UtilityMacros.h"


namespace AMP::Geometry {


MultiGeometry::MultiGeometry( const std::vector<std::shared_ptr<Geometry>> &geom )
    : Geometry(), d_geom( geom )
{
    d_physicalDim = 0;
    for ( const auto &tmp : d_geom )
        d_physicalDim = std::max( d_physicalDim, tmp->getDim() );
}
Point MultiGeometry::nearest( const Point &pos ) const
{
    if ( d_geom.empty() )
        return Point();
    auto p = d_geom[0]->nearest( pos );
    auto d = ( p - pos ).norm();
    for ( size_t i = 1; i < d_geom.size(); i++ ) {
        auto p2 = d_geom[i]->nearest( pos );
        auto d2 = ( p2 - pos ).norm();
        if ( d2 < d ) {
            d = d2;
            p = p2;
        }
    }
    return p;
}
double MultiGeometry::distance( const Point &pos, const Point &dir ) const
{
    double dist = std::numeric_limits<double>::infinity();
    for ( const auto &geom : d_geom ) {
        double d  = geom->distance( pos, dir );
        double d1 = fabs( d );
        double d2 = fabs( dist );
        if ( fabs( d1 - d2 ) < 1e-6 * d1 ) {
            // The two distances are ~ equal, take the inside distance if it exists
            dist = std::min( d, dist );
        } else if ( d1 < d2 ) {
            // Keep the smallest distance
            dist = d;
        }
    }
    return dist;
}
bool MultiGeometry::inside( const Point &pos ) const
{
    bool inside = false;
    for ( const auto &geom : d_geom )
        inside = inside || geom->inside( pos );
    return inside;
}
int MultiGeometry::NSurface() const
{
    AMP_ERROR( "NSurface is not valid for MultiGeometry" );
    return 0;
}
int MultiGeometry::surface( const Point & ) const
{
    AMP_ERROR( "surface is not valid for MultiGeometry" );
    return 0;
}
Point MultiGeometry::surfaceNorm( const Point & ) const
{
    AMP_ERROR( "surfaceNorm is not valid for MultiGeometry" );
    return Point();
}
Point MultiGeometry::centroid() const
{
    Point p( d_physicalDim, { 0, 0, 0 } );
    double V = 0;
    for ( size_t i = 0; i < d_geom.size(); i++ ) {
        if ( d_geom[i]->getGeomType() < getGeomType() ) {
            // Surfaces do not contribute to the centroid
            continue;
        }
        double Vi = d_geom[i]->volume();
        V += Vi;
        p += Vi * d_geom[i]->centroid();
    }
    p *= 1.0 / V;
    return p;
}
std::pair<Point, Point> MultiGeometry::box() const
{
    Point lb = { 1e200, 1e200, 1e200 };
    Point ub = { -1e200, -1e200, -1e200 };
    for ( const auto &geom : d_geom ) {
        auto [lb2, ub2] = geom->box();
        lb.x()          = std::min( lb.x(), lb2.x() );
        lb.y()          = std::min( lb.y(), lb2.y() );
        lb.z()          = std::min( lb.z(), lb2.z() );
        ub.x()          = std::max( ub.x(), ub2.x() );
        ub.y()          = std::max( ub.y(), ub2.y() );
        ub.z()          = std::max( ub.z(), ub2.z() );
    }
    for ( int d = d_physicalDim; d < 3; d++ ) {
        lb[d] = 0;
        ub[d] = 0;
    }
    lb.setNdim( d_physicalDim );
    ub.setNdim( d_physicalDim );
    return std::make_pair( lb, ub );
}
double MultiGeometry::volume() const
{
    double V = 0;
    for ( const auto &geom : d_geom )
        V += geom->volume();
    return V;
}
void MultiGeometry::displace( const double *x )
{
    for ( const auto &geom : d_geom )
        geom->displace( x );
}
std::unique_ptr<AMP::Geometry::Geometry> MultiGeometry::clone() const
{
    std::vector<std::shared_ptr<Geometry>> geom2;
    for ( const auto &geom : d_geom )
        geom2.push_back( geom->clone() );
    return std::make_unique<MultiGeometry>( geom2 );
}


/********************************************************
 * Compare the geometry                                  *
 ********************************************************/
bool MultiGeometry::operator==( const Geometry &rhs ) const
{
    if ( &rhs == this )
        return true;
    auto geom = dynamic_cast<const MultiGeometry *>( &rhs );
    if ( !geom )
        return false;
    return d_geom == geom->d_geom;
}


/********************************************************
 * Write/Read restart data                               *
 ********************************************************/
void MultiGeometry::registerChildObjects( AMP::IO::RestartManager *manager ) const
{
    for ( auto geom : d_geom )
        manager->registerObject( geom );
}
void MultiGeometry::writeRestart( int64_t fid ) const
{
    Geometry::writeRestart( fid );
    std::vector<uint64_t> geomID( d_geom.size() );
    for ( size_t i = 0; i < d_geom.size(); i++ )
        geomID[i] = d_geom[i]->getID();
    writeHDF5( fid, "geomID", geomID );
}
MultiGeometry::MultiGeometry( int64_t fid, AMP::IO::RestartManager *manager ) : Geometry( fid )
{
    std::vector<uint64_t> geomID;
    readHDF5( fid, "geomID", geomID );
    d_geom.resize( geomID.size() );
    for ( size_t i = 0; i < d_geom.size(); i++ )
        d_geom[i] = manager->getData<Geometry>( geomID[i] );
}


} // namespace AMP::Geometry
