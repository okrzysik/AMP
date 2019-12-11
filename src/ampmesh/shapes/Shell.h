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
     * \brief Construct a Shell geometry
     * \param db        Input database
     */
    explicit Shell( std::shared_ptr<AMP::Database> db );

    /**
     * \brief Construct a Shell geometry
     * \param r_min     The minimum radius of the shell
     * \param r_max     The maximum radius of the shell
     */
    explicit Shell( double r_min, double r_max );

    // Functions inherited from Geometry
    virtual std::string getName() const override final { return "Shell"; }
    virtual double distance( const Point &pos, const Point &dir ) const override final;
    virtual bool inside( const Point &pos ) const override final;
    virtual int NSurface() const override final { return 2; }
    virtual int surface( const Point &x ) const override final;
    virtual Point surfaceNorm( const Point &x ) const override final;
    virtual Point logical( const Point &x ) const override final;
    virtual Point physical( const Point &x ) const override final;
    virtual Point centroid() const override final;
    virtual std::pair<Point, Point> box() const override final;
    virtual void displaceMesh( const double *x ) override final;
    virtual std::vector<int> getLogicalGridSize( const std::vector<int> &x ) const override final;
    virtual std::vector<bool> getPeriodicDim() const override final;
    virtual std::vector<int> getLogicalSurfaceIds() const override final;
    virtual std::shared_ptr<AMP::Geometry::Geometry> clone() const override final;

protected:
    // Internal data
    double d_r_min, d_r_max;
    double d_offset[3];

private:
    // Private constuctor
    Shell();
};


} // namespace Geometry
} // namespace AMP

#endif
