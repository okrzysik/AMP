#ifndef included_AMP_Geometry_Shell
#define included_AMP_Geometry_Shell

#include "AMP/ampmesh/Geometry.h"

#include <vector>


namespace AMP {
namespace Geometry {


/**
 * \class Geometry
 * \brief A class used to abstract away geometry information from an application or mesh.
 * \details  This class provides routines for reading, accessing and writing geometries.
 */
class Shell : public Geometry
{
public:
    /**
     * \brief Construct a Shell
     * \param r_min     The minimum radius of the shell
     * \param r_max     The maximum radius of the shell
     */
    explicit Shell( double r_min, double r_max );

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
    double d_r_min, d_r_max;
    std::array<double, 3> d_offset;

private:
    // Private constuctor
    Shell();
};

} // namespace Geometry
} // namespace AMP

#endif
