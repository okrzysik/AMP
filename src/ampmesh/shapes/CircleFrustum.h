#ifndef included_AMP_Geometry_CircleFrustum
#define included_AMP_Geometry_CircleFrustum

#include "AMP/ampmesh/Geometry.h"

#include <array>
#include <vector>


namespace AMP {
namespace Geometry {


/**
 * \class CircleFrustum
 * \brief A geometry for a square frustum
 * \details  This class provides routines for reading, accessing and writing geometries.
 */
class CircleFrustum : public Geometry
{
public:
    /**
     * \brief Construct a CircleFrustum
     * \param r         The the radii of the frustrum (base/top)
     * \param dir       The direction of the pyramid { -x, x, -y, y, -z, z }
     * \param height    The height of the frustrum
     */
    explicit CircleFrustum( const std::array<double, 2> &r, int dir, double height );

public: // Functions inherited from Geometry
    virtual std::string getName() const override final { return "CircleFrustum"; }
    virtual double distance( const Point &pos, const Point &dir ) const override final;
    virtual bool inside( const Point &pos ) const override final;
    virtual int NSurface() const override final { return 3; }
    /**
     * \brief    Get the surface id
     * \details     This function will return the surface id closest to the point
     * \param[in] x     Current position
     * @return          Returns the surface id [-,+,cylinder]
     */
    virtual int surface( const Point &x ) const override final;
    virtual Point surfaceNorm( const Point &x ) const override final;
    virtual Point logical( const Point &x ) const override final;
    virtual Point physical( const Point &x ) const override final;
    virtual Point centroid() const override final;
    virtual std::pair<Point, Point> box() const override final;
    virtual void displaceMesh( const double *x ) override final;

protected:              // Internal data
    uint8_t d_dir;      // The direction of the center axis
    double d_r[2];      // The two radii
    double d_h;         // The height of the frustrum
    double d_offset[3]; // The offset
    Point d_C;          // Apex of cone
    double d_theta;     // Apex angle

private:
    // Private constuctor
    CircleFrustum();
};


} // namespace Geometry
} // namespace AMP

#endif
