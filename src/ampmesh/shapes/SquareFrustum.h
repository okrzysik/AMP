#ifndef included_AMP_Geometry_SquareFrustum
#define included_AMP_Geometry_SquareFrustum

#include "AMP/ampmesh/Geometry.h"

#include <vector>


namespace AMP {
namespace Geometry {


/**
 * \class SquareFrustum
 * \brief A geometry for a square frustum
 * \details  This class provides routines for reading, accessing and writing geometries.
 */
class SquareFrustum : public Geometry
{
public:
    /**
     * \brief Construct a SquareFrustum
     * \param range     The range of the SquareFrustum [xmin, xmax, ymin, ymax, zmin, zmax, ...]
     * \param dir       The direction of the pyramid { -x, x, -y, y, -z, z }
     * \param height    The height of the pyramid
     */
    explicit SquareFrustum( const std::vector<double> &range, int dir, double height );

public: // Functions inherited from Geometry
    virtual uint8_t getDim() const override final { return 3; }
    virtual double distance( const Point &pos, const Point &dir ) const override final;
    virtual bool inside( const Point &pos ) const override final;
    virtual int surface( const Point &x ) const override final;
    virtual Point surfaceNorm( const Point &x ) const override final;
    virtual Point logical( const Point &x ) const override final;
    virtual Point physical( const Point &x ) const override final;
    virtual void displaceMesh( const double *x ) override final;

protected:
    // Internal data
    uint8_t d_dir;
    double d_range[6];        // The bounding box size
    double d_pyramid_size[3]; // The underlying rotated pyramid size
    double d_scale_height;

private:
    // Private constuctor
    SquareFrustum();
};


} // namespace Geometry
} // namespace AMP

#endif
