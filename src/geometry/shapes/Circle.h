#ifndef included_AMP_Geometry_Circle
#define included_AMP_Geometry_Circle

#include "AMP/geometry/LogicalGeometry.h"

#include <vector>


namespace AMP::Geometry {


/**
 * \class Geometry
 * \brief A class used to abstract away geometry information from an application or mesh.
 * \details  This class provides routines for reading, accessing and writing geometries.
 */
class Circle : public LogicalGeometry
{
public:
    /**
     * \brief Construct a Circle geometry
     * \param db        Input database
     */
    explicit Circle( std::shared_ptr<const AMP::Database> db );

    /**
     * \brief Construct a Circle geometry
     * \param R         The circle radius
     */
    explicit Circle( double R );

    // Functions inherited from Geometry
    std::string getName() const override final { return "Circle"; }
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
    double d_R;
    double d_offset[2];

private:
    // Private constructor
    Circle();
};


} // namespace AMP::Geometry

#endif
