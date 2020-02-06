#ifndef included_AMP_MeshGeometry
#define included_AMP_MeshGeometry

#include "AMP/ampmesh/Geometry.h"

#include <memory>
#include <vector>


namespace AMP::Mesh {
class Mesh;        // Forward declare Mesh
class MeshElement; // Forward declare Mesh
} // namespace AMP::Mesh


namespace AMP::Geometry {


/**
 * \class MeshGeometry
 * \brief A class used to abstract away geometry information from an application or mesh.
 * \details  This class provides a geometry implimentation based on a surface mesh
 */
class MeshGeometry final : Geometry
{
public:
    //! Default constructor
    MeshGeometry( std::unique_ptr<AMP::Mesh::Mesh> mesh );

    //! Destructor
    virtual ~MeshGeometry() = default;

    //! Get the name of the geometry
    virtual std::string getName() const override { return "MeshGeometry"; }

    /**
     * \brief    Is the object convex
     * \details  Check if the geometric object is convex
     * @return      Returns true if the object is convex
     */
    virtual bool isConvex() const override;

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
    virtual double distance( const Point &pos, const Point &dir ) const override;

    /**
     * \brief    Is the point in the geometry
     * \details  This function checks if the ray is in the geometry.  If it is on the surface,
     *     it will return true.
     * \param[in] pos   Current position
     * @return          Returns true if the point is inside the geometry (or on the surface)
     */
    virtual bool inside( const Point &pos ) const override;

    /**
     * \brief    Get the number of surfaces
     * \details     This function will return the number of unique surfaces
     * @return          Returns the number of unique surfaces
     */
    virtual int NSurface() const override;

    /**
     * \brief    Get the surface id
     * \details     This function will return the surface id closest to the point
     * \param[in] x     Current position
     * @return          Returns the surface id (0:NSurface-1)
     */
    virtual int surface( const Point &x ) const override;

    /**
     * \brief    Return the outward normal to a surface
     * \details  This function will return the surface id and outward normal to the surface at the
     given point
     * \param[in] x     Current position
     * @return          Returns the surface normal

     */
    virtual Point surfaceNorm( const Point &x ) const override;

    /**
     * \brief    Return the centroid
     * \details  This function will return centroid of the object
     * @return          Returns the physical coordinates
     */
    virtual Point centroid() const override;

    /**
     * \brief    Return the bounding box
     * \details  This function will return the bounding box of the object
     * @return          Returns the bounding box [lb,ub]
     */
    virtual std::pair<Point, Point> box() const override;

    /**
     * \brief    Displace the entire geometry
     * \details  This function will displace the entire geometry by a scalar value.
     *   The displacement vector should be the size of the physical dimension.
     * \param[int] x    Displacement vector
     */
    virtual void displace( const double *x ) override;

    //! Clone the object
    virtual std::unique_ptr<AMP::Geometry::Geometry> clone() const override;


private: // Internal functions
    // Initialize the internal data
    void initialize();

    // Get the nearest element to a point
    AMP::Mesh::MeshElement getNearest( const Point &x ) const;


private: // Internal data
    std::unique_ptr<AMP::Mesh::Mesh> d_mesh;
    std::vector<int> d_surfaceIds;
    Point d_centroid;
};


} // namespace AMP::Geometry


#endif
