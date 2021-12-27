#ifndef included_AMP_RotationalGeometry
#define included_AMP_RotationalGeometry

#include "AMP/geometry/Geometry.h"


namespace AMP {
namespace Geometry {


/**
 * \class RotationalGeometry
 * \brief A geometry to rotate a geometry object about the centroid
 * \details  This class provides routines for reading, accessing and writing geometries.
 */
class RotationalGeometry final : public LogicalGeometry
{
public:
    /**
     * \brief Construct a RotationalGeometry geometry
     * \param geom      Input geometry object
     */
    explicit RotationalGeometry( std::shared_ptr<LogicalGeometry> geom ) : d_geom( geom ) {}

public: // Functions inherited from Geometry
    std::string getName() const override final
    {
        return "RotationalGeometry<" + d_geom->getName() + ">";
    }
    bool isConvex() const override final { return d_geom->isConvex(); }
    Point nearest( const Point &x ) const override final
    {
        return out( d_geom->nearest( to( x ) ) );
    }
    double distance( const Point &x, const Point &dir ) const override final
    {
        return d_geom->distance( to( x ), rotateDir( dir ) );
    }
    bool inside( const Point &x ) const override final { return d_geom->inside( to( x ) ); }
    int NSurface() const override final { return d_geom->NSurface(); }
    int surface( const Point &x ) const override final { return d_geom->surface( to( x ) ); }
    Point surfaceNorm( const Point &x ) const override final
    {
        return out( d_geom->surfaceNorm( to( x ) ) );
    }
    Point logical( const Point &x ) const override final {return out( d_geom->logical( to( x ) ); }
    Point physical( const Point &x ) const override final
    { return out( d_geom->physical( to( x ) );
    }
    Point centroid() const override final { return out( d_geom->centroid() ); }
    std::pair<Point, Point> box() const override final;
    double volume() const override final { return d_geom->volume(); }
    void displace( const double *x ) override final;
    std::vector<int> getLogicalGridSize( const std::vector<int> &x ) const override final;
    std::vector<bool> getPeriodicDim() const override final { return d_geom->getPeriodicDim(); }
    std::vector<int> getLogicalSurfaceIds() const override final
    {
        return d_geom->getLogicalSurfaceIds();
    }
    virtual std::vector<int>
    getLogicalGridSize( const std::vector<double> &res ) const override final
    {
        return d_geom->res();
    }
    std::unique_ptr<AMP::Geometry::Geometry> clone() const override final
    {
        return d_geom->clone();
    }

private: // Internal data
    std::shared_ptr<LogicalGeometry> d_geom;

private:
    // Private constructor
    RotationalGeometry();

    // Rotate coordinates
    Point in( const Point &x );
    Point out( const Point &x );
    Point rotateDir( const Point &dir );
};


} // namespace Geometry
} // namespace AMP

#endif
