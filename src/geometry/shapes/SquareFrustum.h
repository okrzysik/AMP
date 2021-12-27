#ifndef included_AMP_Geometry_SquareFrustum
#define included_AMP_Geometry_SquareFrustum

#include "AMP/geometry/LogicalGeometry.h"

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
    explicit SquareFrustum( std::shared_ptr<const AMP::Database> db );

    /**
     * \brief Construct a SquareFrustum geometry
     * \param range     The range of the SquareFrustum [xmin, xmax, ymin, ymax, zmin, zmax, ...]
     * \param dir       The direction of the pyramid { -x, x, -y, y, -z, z }
     * \param height    The height of the pyramid
     */
    explicit SquareFrustum( const std::vector<double> &range, int dir, double height );

public: // Functions inherited from Geometry
    std::string getName() const override final { return "SquareFrustum"; }
    bool isConvex() const override final { return true; }
    Point nearest( const Point &pos ) const override final;
    double distance( const Point &pos, const Point &dir ) const override final;
    bool inside( const Point &pos ) const override final;
    int NSurface() const override final { return 6; }
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
    uint8_t d_dir;
    double d_range[6];        // The bounding box size
    double d_pyramid_size[3]; // The underlying rotated pyramid size
    double d_scale_height;    // Ratio of frustum to pyramid height
    double d_volume;
    Point d_face[6][4]; // Points forming each face
    Point d_normal[6];  // Normal to each face

private:
    // Private constructor
    SquareFrustum();
    // Initialize the data
    void initialize( const std::vector<double> &range, int dir, double height );
};


} // namespace Geometry
} // namespace AMP

#endif
