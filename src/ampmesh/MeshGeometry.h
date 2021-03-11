#ifndef included_AMP_MeshGeometry
#define included_AMP_MeshGeometry

#include "AMP/ampmesh/Geometry.h"
#include "AMP/ampmesh/Mesh.h"
#include "AMP/ampmesh/MeshUtilities.h"

#include <memory>
#include <vector>


namespace AMP::Geometry {


/**
 * \class MeshGeometry
 * \brief A class used to abstract away geometry information from an application or mesh.
 * \details  This class provides a geometry implimentation based on a surface mesh
 */
class MeshGeometry final : public Geometry
{
public:
    //! Default constructor
    MeshGeometry( std::shared_ptr<AMP::Mesh::Mesh> mesh );

    //! Destructor
    virtual ~MeshGeometry() = default;

    //! Copy constructor
    MeshGeometry( const MeshGeometry & ) = delete;

    //! Assignment operator
    MeshGeometry &operator=( const MeshGeometry & ) = delete;

    //! Get the name of the geometry
    std::string getName() const override { return "MeshGeometry"; }

    /**
     * \brief    Is the object convex
     * \details  Check if the geometric object is convex
     * @return      Returns true if the object is convex
     */
    bool isConvex() const override;

    /**
     * \brief    Calculate the nearest point on the surface
     * \details  This function computes the nearest point on the surface
     * \param[in] pos   Current position of ray
     * @return          Returns the nearest surface point
     */
    Point nearest( const Point &pos ) const override;

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
    double distance( const Point &pos, const Point &dir ) const override;

    /**
     * \brief    Is the point in the geometry
     * \details  This function checks if the ray is in the geometry.  If it is on the surface,
     *     it will return true.
     * \param[in] pos   Current position
     * @return          Returns true if the point is inside the geometry (or on the surface)
     */
    bool inside( const Point &pos ) const override;

    /**
     * \brief    Get the number of surfaces
     * \details     This function will return the number of unique surfaces
     * @return          Returns the number of unique surfaces
     */
    int NSurface() const override;

    /**
     * \brief    Get the surface id
     * \details     This function will return the surface id closest to the point
     * \param[in] x     Current position
     * @return          Returns the surface id (0:NSurface-1)
     */
    int surface( const Point &x ) const override;

    /**
     * \brief    Return the outward normal to a surface
     * \details  This function will return the surface id and outward normal to the surface at the
     given point
     * \param[in] x     Current position
     * @return          Returns the surface normal

     */
    Point surfaceNorm( const Point &x ) const override;

    /**
     * \brief    Return the centroid
     * \details  This function will return centroid of the object
     * @return          Returns the physical coordinates
     */
    Point centroid() const override;

    /**
     * \brief    Return the bounding box
     * \details  This function will return the bounding box of the object
     * @return          Returns the bounding box [lb,ub]
     */
    std::pair<Point, Point> box() const override;

    /**
     * \brief    Return the volume
     * \details  This function will return the interior volume of the object
     * @return          Returns the volume
     */
    double volume() const override;

    /**
     * \brief    Displace the entire geometry
     * \details  This function will displace the entire geometry by a scalar value.
     *   The displacement vector should be the size of the physical dimension.
     * \param[in] x     Displacement vector
     */
    void displace( const double *x ) override;

    //! Clone the object
    std::unique_ptr<AMP::Geometry::Geometry> clone() const override;

    //! Get the mesh
    const AMP::Mesh::Mesh &getMesh() const { return *d_mesh; }

    //! Check if two geometries are equal
    bool operator==( const Geometry &rhs ) const override;


private: // Internal functions
    // Get the nearest element to a point
    std::vector<AMP::Mesh::MeshElement> getNearestElements( const Point &x ) const;
    std::pair<AMP::Mesh::MeshElement, Point> getNearestPoint( const Point &x ) const;


private: // Internal data
    std::shared_ptr<AMP::Mesh::Mesh> d_mesh;
    std::vector<int> d_surfaceIds;
    AMP::Mesh::ElementFinder d_find;
};


} // namespace AMP::Geometry


#endif
