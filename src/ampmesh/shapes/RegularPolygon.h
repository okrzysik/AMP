#ifndef included_AMP_Geometry_RegularPolygon
#define included_AMP_Geometry_RegularPolygon

#include "AMP/ampmesh/LogicalGeometry.h"

#include <vector>


namespace AMP {
namespace Geometry {


/**
 * \class RegularPolygon
 * \brief A class used to abstract away geometry information from an application or mesh.
 * \details  This class provides routines for reading, accessing and writing regular polygons.
 */
class RegularPolygon : public LogicalGeometry
{
public:
    /**
     * \brief Construct a regular polygon geometry
     * \param db        Input database
     */
    explicit RegularPolygon( std::shared_ptr<AMP::Database> db );

    /**
     * \brief Construct a regular polygon geometry
     * \param N         The number of sides
     * \param R         The circumcircle radius
     */
    explicit RegularPolygon( int N, double R );

    // Functions inherited from Geometry
    std::string getName() const override { return "RegularPolygon"; }
    bool isConvex() const override final { return true; }
    Point nearest( const Point &pos ) const override final;
    double distance( const Point &pos, const Point &dir ) const override final;
    bool inside( const Point &pos ) const override final;
    int NSurface() const override final;
    int surface( const Point & ) const override final;
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
    int d_N;
    double d_R;
    double d_offset[2];
    std::vector<std::array<double, 2>> d_verticies;
    std::vector<Point> d_norm;

protected:
    // Internal functions
    std::tuple<Point, double, int> nearest2( const Point &pos ) const;
    void computeNorms();


private:
    // Private constuctor
    RegularPolygon();
};


} // namespace Geometry
} // namespace AMP

#endif
