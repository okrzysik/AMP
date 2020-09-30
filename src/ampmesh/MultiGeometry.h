#ifndef included_AMP_MultiGeometry
#define included_AMP_MultiGeometry

#include "AMP/ampmesh/LogicalGeometry.h"


namespace AMP {
namespace Geometry {


/**
 * \class MultiGeometry
 * \brief A concrete geometry class for a multi-mesh
 *
 * \details  This class provides routines for creating and accessing geometry functions for a
 * multi-mesh.
 */
class MultiGeometry : public Geometry
{
public:
    /**
     * \brief Contructor to create a MultiMesh from existing meshes
     * \details  This constructor takes a list of meshes and a communicator
     *    and generates the appropriate multimesh
     * \param geom      Sub-geometries
     */
    explicit MultiGeometry( const std::vector<Geometry::shared_ptr> &geom );


    //! Destructor
    virtual ~MultiGeometry() = default;


    //! Access the geometries
    inline const auto &getGeometries() const { return d_geom; }

public: // Functions inherited from Geometry
    std::string getName() const override final { return "MultiGeometry"; }
    bool isConvex() const override final { return false; }
    Point nearest( const Point &pos ) const override final;
    double distance( const Point &pos, const Point &dir ) const override final;
    bool inside( const Point &pos ) const override final;
    int NSurface() const override final;
    int surface( const Point &x ) const override final;
    Point surfaceNorm( const Point &x ) const override final;
    Point centroid() const override final;
    std::pair<Point, Point> box() const override final;
    double volume() const override final;
    void displace( const double *x ) override final;
    std::unique_ptr<AMP::Geometry::Geometry> clone() const override final;

private:
    //! Empty constructor for a mesh
    MultiGeometry() = delete;

    std::vector<Geometry::shared_ptr> d_geom;
};

} // namespace Geometry
} // namespace AMP

#endif
