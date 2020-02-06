#ifndef included_AMP_Geometry_Circle
#define included_AMP_Geometry_Circle

#include "AMP/ampmesh/LogicalGeometry.h"

#include <vector>


namespace AMP {
namespace Geometry {


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
    explicit Circle( std::shared_ptr<AMP::Database> db );

    /**
     * \brief Construct a Circle geometry
     * \param range     The range of the Circle [xmin, xmax, ymin, ymax, zmin, zmax, ...]
     */
    explicit Circle( double R );

    // Functions inherited from Geometry
    virtual std::string getName() const override final { return "Circle"; }
    virtual bool isConvex() const override final { return true; }
    virtual double distance( const Point &pos, const Point &dir ) const override final;
    virtual bool inside( const Point &pos ) const override final;
    virtual int NSurface() const override final { return 1; }
    virtual int surface( const Point & ) const override final { return 0; }
    virtual Point surfaceNorm( const Point &x ) const override final;
    virtual Point logical( const Point &x ) const override final;
    virtual Point physical( const Point &x ) const override final;
    virtual Point centroid() const override final;
    virtual std::pair<Point, Point> box() const override final;
    virtual void displace( const double *x ) override final;
    virtual std::vector<int> getLogicalGridSize( const std::vector<int> &x ) const override final;
    virtual std::vector<bool> getPeriodicDim() const override final;
    virtual std::vector<int> getLogicalSurfaceIds() const override final;
    virtual std::unique_ptr<AMP::Geometry::Geometry> clone() const override final;

protected:
    // Internal data
    double d_R;
    double d_offset[2];

private:
    // Private constuctor
    Circle();
};


} // namespace Geometry
} // namespace AMP

#endif
