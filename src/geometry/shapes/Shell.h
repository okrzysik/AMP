#ifndef included_AMP_Geometry_Shell
#define included_AMP_Geometry_Shell

#include "AMP/geometry/LogicalGeometry.h"

#include <vector>


namespace AMP::Geometry {


/**
 * \class Shell
 * \brief A class used to abstract away geometry information from an application or mesh.
 * \details  This class provides routines for reading, accessing and writing a spherical shell.
 */
class Shell : public LogicalGeometry
{
public:
    /**
     * \brief Construct a Shell geometry
     * \param db        Input database
     */
    explicit Shell( std::shared_ptr<const AMP::Database> db );

    /**
     * \brief Construct a Shell geometry
     * \param r_min     The minimum radius of the shell
     * \param r_max     The maximum radius of the shell
     */
    explicit Shell( double r_min, double r_max );

    // Functions inherited from Geometry
    std::string getName() const override final { return "Shell"; }
    bool isConvex() const override final { return true; }
    Point nearest( const Point &pos ) const override final;
    double distance( const Point &pos, const Point &dir ) const override final;
    bool inside( const Point &pos ) const override final;
    int NSurface() const override final { return 2; }
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
    std::unique_ptr<AMP::Geometry::Geometry> clone() const override final;
    bool operator==( const Geometry &rhs ) const override final;

protected:
    // Internal data
    double d_r_min, d_r_max;
    double d_offset[3];

private:
    // Private constructor
    Shell();
};


} // namespace AMP::Geometry

#endif
