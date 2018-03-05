#ifndef included_AMP_Geometry_Tube
#define included_AMP_Geometry_Tube

#include "AMP/ampmesh/Geometry.h"

#include <vector>


namespace AMP {
namespace Geometry {


/**
 * \class Geometry
 * \brief A class used to abstract away geometry information from an application or mesh.
 * \details  This class provides routines for reading, accessing and writing geometries.
 */
class Tube : public Geometry
{
public:
    /**
     * \brief Construct a SphereSurface
     * \param r_min     The minimum radius
     * \param r_max     The maximum radius
     * \param z_min     The minimum z coordinate
     * \param z_max     The maximum z coordinate
     */
    explicit Tube( double r_min, double r_max, double z_min, double z_max );

    // Functions inherited from Geometry
    virtual uint8_t getDim() const override final { return 3; }
    virtual double distance( const Point<double> &pos,
                             const Point<double> &dir ) const override final;
    virtual bool inside( const Point<double> &pos ) const override final;
    virtual int surface( const Point<double> &x ) const override final;
    virtual Point<double> surfaceNorm( const Point<double> &x ) const override final;
    virtual Point<double> logical( const Point<double> &x ) const override final;
    virtual Point<double> physical( const Point<double> &x ) const override final;
    virtual void displaceMesh( const double *x ) override final;

protected:
    // Internal data
    double d_r_min, d_r_max, d_z_min, d_z_max;
    double d_offset[3];

private:
    // Private constuctor
    Tube();
};


} // namespace Geometry
} // namespace AMP

#endif
