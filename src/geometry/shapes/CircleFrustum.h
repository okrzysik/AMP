#ifndef included_AMP_Geometry_CircleFrustum
#define included_AMP_Geometry_CircleFrustum

#include "AMP/geometry/LogicalGeometry.h"

#include <array>
#include <vector>


namespace AMP::Geometry {


/**
 * \class CircleFrustum
 * \brief A geometry for a circular frustum
 * \details  This class provides routines for reading, accessing and writing geometries.
 */
class CircleFrustum final : public LogicalGeometry
{
public:
    /**
     * \brief Construct a CircleFrustum geometry
     * \param db        Input database
     */
    explicit CircleFrustum( std::shared_ptr<const AMP::Database> db );

    /**
     * \brief Construct a CircleFrustum geometry
     * \param r         The the radii of the frustum (base/top)
     * \param dir       The direction of the pyramid { -x, x, -y, y, -z, z }
     * \param height    The height of the frustum
     */
    explicit CircleFrustum( const std::array<double, 2> &r, int dir, double height );

    //! Construct from restart
    CircleFrustum( int64_t );

    //! Copy contructor
    CircleFrustum( const CircleFrustum & ) = default;

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
    ArraySize getLogicalGridSize( const ArraySize &x ) const override final;
    ArraySize getLogicalGridSize( const std::vector<double> &res ) const override final;
    std::unique_ptr<AMP::Geometry::Geometry> clone() const override final;
    bool operator==( const Geometry &rhs ) const override final;
    void writeRestart( int64_t ) const override;

protected:                                        // Internal data
    uint8_t d_dir                  = 0;           // The direction of the center axis
    double d_h                     = 0;           // The height of the frustrum
    std::array<double, 2> d_r      = { 0, 0 };    // The two radii
    std::array<double, 3> d_offset = { 0, 0, 0 }; // The offset
    double d_theta                 = 0;           // Apex angle

private:
    // Private constructor
    CircleFrustum() = delete;
    // Initialize the data
    void initialize( int dir, const std::array<double, 2> &r, double h );
    // Convert a point/angle to the reference frame (does not subtract offset)
    Point convertToReference( const Point & ) const;
    // Convert a point/angle from the reference frame (does not add offset)
    Point convertFromReference( const Point & ) const;
};


} // namespace AMP::Geometry

#endif
