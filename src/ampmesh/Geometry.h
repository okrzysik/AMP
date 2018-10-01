#ifndef included_AMP_Geometry
#define included_AMP_Geometry

#include "AMP/ampmesh/MeshPoint.h"
#include "AMP/utils/shared_ptr.h"

#include <vector>


namespace AMP {
namespace Geometry {


using Point = AMP::Mesh::MeshPoint<double>;


/**
 * \class Geometry
 * \brief A class used to abstract away geometry information from an application or mesh.
 * \details  This class provides routines for reading, accessing and writing geometries.
 */
class Geometry
{
public:
    //! Destructor
    virtual ~Geometry() {}

    /**
     *\typedef shared_ptr
     *\brief  Name for the shared pointer.
     *\details  Use this typedef for a reference counted pointer to a geometry object.
     */
    typedef AMP::shared_ptr<AMP::Geometry::Geometry> shared_ptr;

    /**
     *\typedef const_shared_ptr
     *\brief  Name for the const shared pointer.
     *\details  Use this typedef for a reference counted pointer to a geometry object.
     */
    typedef AMP::shared_ptr<const AMP::Geometry::Geometry> const_shared_ptr;


    /**
     * \brief    Get the number of dimensions for the object
     * \details  This function returns the number of physical dimensions for the geometry
     * @return      Returns the distance to the nearest surface (intersection = pos + dir*distance)
     */
    virtual uint8_t getDim() const = 0;

    /**
     * \brief    Calculate the distance to the object given a ray
     * \details  This function computes the distance to the object given a ray.
     *     If the ray is inside the object, this distance is negitive.  If the
     *     ray will never intersect the object, this distance is inf.
     * \param[in] pos   Current position of ray
     * \param[in] dir   Direction of ray (should be normalized for most uses)
     * @return          Returns the distance to the nearest surface
     *                  (intersection = pos + dir*distance)
     */
    virtual double distance( const Point &pos, const Point &dir ) const = 0;

    /**
     * \brief    Is the point in the geometry
     * \details  This function checks if the ray is in the geometry.  If it is on the surface,
     *     it will return true.
     * \param[in] pos   Current position
     * @return          Returns true if the point is inside the geometry (or on the surface)
     */
    virtual bool inside( const Point &pos ) const = 0;

    /**
     * \brief    Get the surface id
     * \details     This function will return the surface id closest to the point
     * \param[in] x     Current position
     * @return          Returns the surface id
     */
    virtual int surface( const Point &x ) const = 0;

    /**
     * \brief    Return the outward normal to a surface
     * \details  This function will return the surface id and outward normal to the surface at the
     given point
     * \param[in] x     Current position
     * @return          Returns the surface normal

     */
    virtual Point surfaceNorm( const Point &x ) const = 0;

    /**
     * \brief    Return the logical position
     * \details  This function will return logical coordinate given the physical coordinates
     * \param[in] x     Physical coordinate of the point
     * @return          Returns the logical coordinates
     */
    virtual Point logical( const Point &x ) const = 0;

    /**
     * \brief    Return the physical position
     * \details  This function will return physical coordinate given the logical coordinates
     * \param[in] x     Logical coordinate of the point
     * @return          Returns the physical coordinates
     */
    virtual Point physical( const Point &x ) const = 0;

    /**
     * \brief    Displace the entire mesh
     * \details  This function will displace the entire mesh by a scalar value.
     *   This function is a blocking call for the mesh communicator, and requires
     *   the same value on all processors.  The displacement vector should be the
     *   size of the physical dimension.
     * \param[int] x    Displacement vector
     */
    virtual void displaceMesh( const double *x ) = 0;


protected:
    //!  Empty constructor for the base class
    Geometry() {}
};


} // namespace Geometry
} // namespace AMP

#endif
