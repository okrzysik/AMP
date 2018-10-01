#ifndef included_AMP_Geometry_SphereSurface
#define included_AMP_Geometry_SphereSurface

#include "AMP/ampmesh/Geometry.h"

#include <vector>


namespace AMP {
namespace Geometry {


/**
 * \class Geometry
 * \brief A class used to abstract away geometry information from an application or mesh.
 * \details  This class provides routines for reading, accessing and writing geometries.
 */
class SphereSurface : public Geometry
{
public:
    /**
     * \brief Construct a SphereSurface
     * \param range     The range of the SphereSurface [xmin, xmax, ymin, ymax, zmin, zmax, ...]
     */
    explicit SphereSurface( double R );

    // Functions inherited from Geometry
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
    double d_r;
    double d_offset[3];

private:
    // Private constuctor
    SphereSurface();
};


} // namespace Geometry
} // namespace AMP

#endif
