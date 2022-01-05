#ifndef included_AMP_Geometry_SphereSurface
#define included_AMP_Geometry_SphereSurface

#include "AMP/geometry/LogicalGeometry.h"

#include <vector>


namespace AMP::Geometry {


/**
 * \class Geometry
 * \brief A class used to abstract away geometry information from an application or mesh.
 * \details  This class provides routines for reading, accessing and writing geometries.
 */
class SphereSurface : public LogicalGeometry
{
public:
    /**
     * \brief Construct a SphereSurface geometry
     * \param db        Input database
     */
    explicit SphereSurface( std::shared_ptr<const AMP::Database> db );

    /**
     * \brief Construct a SphereSurface geometry
     * \param R     The radius
     */
    explicit SphereSurface( double R );

    // Functions inherited from Geometry
    std::string getName() const override final { return "SphereSurface"; }
    AMP::Mesh::GeomType getGeomType() const override final { return AMP::Mesh::GeomType::Face; }
    bool isConvex() const override final { return true; }
    Point nearest( const Point &pos ) const override final;
    double distance( const Point &pos, const Point &dir ) const override final;
    bool inside( const Point &pos ) const override final;
    int NSurface() const override final { return 1; }
    int surface( const Point & ) const override final { return 0; }
    Point surfaceNorm( const Point &x ) const override final;
    Point logical( const Point &x ) const override final;
    Point physical( const Point &x ) const override final;
    Point centroid() const override final;
    std::pair<Point, Point> box() const override final;
    double volume() const override final;
    void displace( const double *x ) override final;
    std::vector<int> getLogicalGridSize( const std::vector<int> &x ) const override final;
    virtual std::vector<int>
    getLogicalGridSize( const std::vector<double> &res ) const override final;
    std::vector<bool> getPeriodicDim() const override final;
    std::vector<int> getLogicalSurfaceIds() const override final;
    std::unique_ptr<AMP::Geometry::Geometry> clone() const override final;
    bool operator==( const Geometry &rhs ) const override final;

protected:
    // Internal data
    double d_r;
    double d_offset[3];

private:
    // Private constructor
    SphereSurface();
};


} // namespace AMP::Geometry

#endif
