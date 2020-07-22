#ifndef included_AMP_Geometry_SquareFrustum
#define included_AMP_Geometry_SquareFrustum

#include "AMP/ampmesh/LogicalGeometry.h"

#include <vector>


namespace AMP {
namespace Geometry {


/**
 * \class SquareFrustum
 * \brief A geometry for a square frustum
 * \details  This class provides routines for reading, accessing and writing geometries.
 */
class SquareFrustum : public LogicalGeometry
{
public:
    /**
     * \brief Construct a SquareFrustum geometry
     * \param db        Input database
     */
    explicit SquareFrustum( std::shared_ptr<AMP::Database> db );

    /**
     * \brief Construct a SquareFrustum geometry
     * \param range     The range of the SquareFrustum [xmin, xmax, ymin, ymax, zmin, zmax, ...]
     * \param dir       The direction of the pyramid { -x, x, -y, y, -z, z }
     * \param height    The height of the pyramid
     */
    explicit SquareFrustum( const std::vector<double> &range, int dir, double height );

public: // Functions inherited from Geometry
    virtual std::string getName() const override final { return "SquareFrustum"; }
    virtual bool isConvex() const override final { return true; }
    virtual Point nearest( const Point &pos ) const override final;
    virtual double distance( const Point &pos, const Point &dir ) const override final;
    virtual bool inside( const Point &pos ) const override final;
    virtual int NSurface() const override final { return 6; }
    virtual int surface( const Point &x ) const override final;
    virtual Point surfaceNorm( const Point &x ) const override final;
    virtual Point logical( const Point &x ) const override final;
    virtual Point physical( const Point &x ) const override final;
    virtual Point centroid() const override final;
    virtual std::pair<Point, Point> box() const override final;
    virtual double volume() const override final;
    virtual void displace( const double *x ) override final;
    virtual std::vector<int> getLogicalGridSize( const std::vector<int> &x ) const override final;
    virtual std::vector<int>
    getLogicalGridSize( const std::vector<double> &res ) const override final;
    virtual std::vector<bool> getPeriodicDim() const override final;
    virtual std::vector<int> getLogicalSurfaceIds() const override final;
    virtual std::unique_ptr<AMP::Geometry::Geometry> clone() const override final;

protected:
    // Internal data
    uint8_t d_dir;
    double d_range[6];        // The bounding box size
    double d_pyramid_size[3]; // The underlying rotated pyramid size
    double d_scale_height;    // Ratio of frustum to pyramid height
    double d_volume;
    Point d_face[6][4]; // Points forming each face
    Point d_normal[6];  // Normal to each face

private:
    // Private constuctor
    SquareFrustum();
    // Initialize the data
    void initialize( const std::vector<double> &range, int dir, double height );
};


} // namespace Geometry
} // namespace AMP

#endif
