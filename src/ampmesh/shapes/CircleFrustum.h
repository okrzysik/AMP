#ifndef included_AMP_Geometry_CircleFrustum
#define included_AMP_Geometry_CircleFrustum

#include "AMP/ampmesh/LogicalGeometry.h"

#include <array>
#include <vector>


namespace AMP {
namespace Geometry {


/**
 * \class CircleFrustum
 * \brief A geometry for a square frustum
 * \details  This class provides routines for reading, accessing and writing geometries.
 */
class CircleFrustum : public LogicalGeometry
{
public:
    /**
     * \brief Construct a CircleFrustum geometry
     * \param db        Input database
     */
    explicit CircleFrustum( std::shared_ptr<const AMP::Database> db );

    /**
     * \brief Construct a CircleFrustum geometry
     * \param r         The the radii of the frustrum (base/top)
     * \param dir       The direction of the pyramid { -x, x, -y, y, -z, z }
     * \param height    The height of the frustrum
     */
    explicit CircleFrustum( const std::array<double, 2> &r, int dir, double height );

public: // Functions inherited from Geometry
    std::string getName() const override final { return "CircleFrustum"; }
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
    std::vector<bool> getPeriodicDim() const override final;
    std::vector<int> getLogicalSurfaceIds() const override final;
    virtual std::vector<int>
    getLogicalGridSize( const std::vector<double> &res ) const override final;
    std::unique_ptr<AMP::Geometry::Geometry> clone() const override final;
    bool operator==( const Geometry &rhs ) const override final;

protected:              // Internal data
    uint8_t d_dir;      // The direction of the center axis
    double d_r[2];      // The two radii
    double d_h;         // The height of the frustrum
    double d_offset[3]; // The offset
    Point d_C;          // Apex of cone
    double d_theta;     // Apex angle

private:
    // Private constructor
    CircleFrustum();
};


} // namespace Geometry
} // namespace AMP

#endif
