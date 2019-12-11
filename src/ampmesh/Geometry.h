#ifndef included_AMP_Geometry
#define included_AMP_Geometry

#include "AMP/ampmesh/MeshPoint.h"
#include <memory>

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

    //! Get the name of the geometry
    virtual std::string getName() const = 0;

    /**
     *\typedef shared_ptr
     *\brief  Name for the shared pointer.
     *\details  Use this typedef for a reference counted pointer to a geometry object.
     */
    typedef std::shared_ptr<AMP::Geometry::Geometry> shared_ptr;

    /**
     *\typedef const_shared_ptr
     *\brief  Name for the const shared pointer.
     *\details  Use this typedef for a reference counted pointer to a geometry object.
     */
    typedef std::shared_ptr<const AMP::Geometry::Geometry> const_shared_ptr;


    /**
     * \brief    Get the number of dimensions for the object
     * \details  This function returns the number of physical dimensions for the geometry
     * @return      Returns the number of physical dimensions
     */
    inline uint8_t getDim() const { return d_physicalDim; }

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
     * \brief    Get the number of surfaces
     * \details     This function will return the number of unique surfaces
     * @return          Returns the number of unique surfaces
     */
    virtual int NSurface() const = 0;

    /**
     * \brief    Get the surface id
     * \details     This function will return the surface id closest to the point
     * \param[in] x     Current position
     * @return          Returns the surface id (0:NSurface-1)
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
     * \brief    Return the centroid
     * \details  This function will return centroid of the object
     * @return          Returns the physical coordinates
     */
    virtual Point centroid() const = 0;

    /**
     * \brief    Return the bounding box
     * \details  This function will return the bounding box of the object
     * @return          Returns the bounding box [lb,ub]
     */
    virtual std::pair<Point, Point> box() const = 0;

    /**
     * \brief    Displace the entire geometry
     * \details  This function will displace the entire geometry by a scalar value.
     *   The displacement vector should be the size of the physical dimension.
     * \param[int] x    Displacement vector
     */
    virtual void displaceMesh( const double *x ) = 0;

    /**
     * \brief    Is the geometry logically rectangular
     * \details  This function will return true if the underlying geometry is logically
     *    rectangular.  This would mean logical() and physical() are valid operations.
     */
    inline bool isLogical() const { return d_logicalDim != 0; }

    /**
     * \brief    Return the number of logical dimensions
     * \details  This function will return the number of logical dimensions
     *    of the underlying geometry.  If the geometry is not logically rectangular
     *    this function should return 0.
     */
    inline uint8_t getLogicalDim() const { return d_logicalDim; }

    /**
     * \brief    Return the logical grid size
     * \details  This function will return the dimensions of a logical grid
     *    given a size that makes sense for the object.
     *    If the coordinates cannot map to a logical grid, this function should throw
     *    a runtime exception.
     * \param[int] x    Input size
     * @return          Return the logical boundary ids (2*logicalDim)
     */
    virtual std::vector<int> getLogicalGridSize( const std::vector<int> &x ) const = 0;

    /**
     * \brief    Return the logical grid periodic dimensions
     * \details  This function will return a vector indicating which logical grid
     *    dimensions are periodic.  If the coordinates cannot map to a logical grid,
     *    this function should throw a runtime exception.
     * @return          Return the periodic dimensions
     */
    virtual std::vector<bool> getPeriodicDim() const = 0;

    /**
     * \brief    Return the surface ids for the logical boundaries
     * \details  This function will return the surface ids for each logical boundary.
     *    If a logical boundary does not map to a surface, it will return -1.
     *    If the coordinates cannot map to a logical grid, this function should
     *    throw a runtime exception.
     * @return          Return the logical boundary ids (2*logicalDim)
     */
    virtual std::vector<int> getLogicalSurfaceIds() const = 0;

    //! Clone the object
    virtual std::shared_ptr<AMP::Geometry::Geometry> clone() const = 0;

public:
    /**
     * \brief   Create a geometry
     * \details  This function will create a geometry based on
     *   the input database.
     * \param params Parameters for constructing a geometry from an input database
     */
    static std::shared_ptr<AMP::Geometry::Geometry>
    buildGeometry( std::shared_ptr<AMP::Database> db );


protected:
    //!  Empty constructor for the base class
    Geometry() : d_physicalDim( 0 ), d_logicalDim( 0 ) {}

    // Delete copy constructors
    Geometry( Geometry && )      = delete;
    Geometry( const Geometry & ) = default;
    Geometry &operator=( Geometry && ) = delete;
    Geometry &operator=( const Geometry & ) = delete;

protected: // Helper functions
    // Compute the normal to the plane formed by 3 points
    static inline Point normal( const Point &a, const Point &b, const Point &c )
    {
        return normalize( cross( b - a, c - a ) );
    }

protected: // Internal data
    uint8_t d_physicalDim;
    uint8_t d_logicalDim;
};


} // namespace Geometry
} // namespace AMP

#endif
