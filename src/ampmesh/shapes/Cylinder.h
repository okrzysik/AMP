#ifndef included_AMP_Geometry_Cylinder
#define included_AMP_Geometry_Cylinder

#include "AMP/ampmesh/LogicalGeometry.h"

#include <vector>


namespace AMP {
namespace Geometry {


/**
 * \class Geometry
 * \brief A class used to abstract away geometry information from an application or mesh.
 * \details  This class provides routines for reading, accessing and writing geometries.
 */
class Cylinder : public LogicalGeometry
{
public:
    /**
     * \brief Construct a Cylinder geometry
     * \param db        Input database
     */
    explicit Cylinder( std::shared_ptr<AMP::Database> db );

    /**
     * \brief Construct a Cylinder geometry
     * \param r         The radius of the cylinder
     * \param z_min     The lower z-coordinate
     * \param z_max     The upper z-coordinate
     */
    explicit Cylinder( double r, double z_min, double z_max );

    // Functions inherited from Geometry
    std::string getName() const override final { return "Cylinder"; }
    bool isConvex() const override final { return true; }
    Point nearest( const Point &pos ) const override final;
    double distance( const Point &pos, const Point &dir ) const override final;
    bool inside( const Point &pos ) const override final;
    int NSurface() const override final { return 3; }
    int surface( const Point &x ) const override final;
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
    double d_z_min;
    double d_z_max;
    double d_offset[3];

private:
    // Private constuctor
    Cylinder();
};


} // namespace Geometry
} // namespace AMP

#endif
