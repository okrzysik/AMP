#ifndef included_AMP_Geometry_Box
#define included_AMP_Geometry_Box

#include "ampmesh/Geometry.h"

#include <vector>


namespace AMP {
namespace Geometry {


/**
 * \class Geometry
 * \brief A class used to abstract away geometry information from an application or mesh.
 * \details  This class provides routines for reading, accessing and writing geometries.
 */
class Box final: public Geometry
{
public:

    /**
     * \brief Construct a box
     * \param range     The range of the box [xmin, xmax, ymin, ymax, zmin, zmax, ...]
     */
    Box( const std::vector<double>& range ): d_range(range) {}

    /**
     * \brief    Calculate the distance to the object given a ray
     * \details  This function computes the distance to the object given a ray.
     *     If the ray is inside the object, this distance is negitive.  If the 
     *     ray will never intersect the object, this distance is inf.  
     * \param pos   Current position of ray
     * \param dir   Direction of ray (should be normalized for most uses)
     * @return      Returns the distance to the nearest surface (intersection = pos + dir*distance)
     */
    virtual double distance( const std::initializer_list<double>& pos, const std::initializer_list<double>& dir ) const override;

    /**
     * \brief    Is the point in the geometry
     * \details  This function checks if the ray is in the geometry.  If it is on the surface,
     *     it will return true.
     * \param pos   Current position of ray
     * @return      Returns true if the point is inside the geometry (or on the surface)
     */
    virtual bool inside( const std::initializer_list<double>& pos ) const override;

    /**
     * \brief    Displace the entire mesh
     * \details  This function will displace the entire mesh by a scalar value.
     *   This function is a blocking call for the mesh communicator, and requires
     *   the same value on all processors.  The displacement vector should be the
     *   size of the physical dimension.
     * \param x  Displacement vector
     */
    virtual void displaceMesh( const std::vector<double> &x ) override;

    /**
     * \brief Return the box
     * @return      The range of the box [xmin, xmax, ymin, ymax, zmin, zmax, ...]
     */
    inline std::vector<double> range() const { return d_range; }

private:

    // Private constuctor
    Box();

    std::vector<double> d_range;

};

} // Mesh namespace
} // Geometry namespace

#endif
