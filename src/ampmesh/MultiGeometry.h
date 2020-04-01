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
     * \param name      Name of the new mesh
     * \param comm      Desired communicator for the multimesh
     * \param meshes    Meshes to be used as part of the multimesh
     */
    explicit MultiGeometry( const std::vector<Geometry::shared_ptr> &geom );


    //! Destructor
    virtual ~MultiGeometry() = default;


    //! Access the geometries
    inline const auto &getGeometries() const { return d_geom; }

public: // Functions inherited from Geometry
    virtual std::string getName() const override final { return "MultiGeometry"; }
    virtual bool isConvex() const override final { return false; }
    virtual Point nearest( const Point &pos ) const override final;
    virtual double distance( const Point &pos, const Point &dir ) const override final;
    virtual bool inside( const Point &pos ) const override final;
    virtual int NSurface() const override final;
    virtual int surface( const Point &x ) const override final;
    virtual Point surfaceNorm( const Point &x ) const override final;
    virtual Point centroid() const override final;
    virtual std::pair<Point, Point> box() const override final;
    virtual void displace( const double *x ) override final;
    virtual std::unique_ptr<AMP::Geometry::Geometry> clone() const override final;

private:
    //! Empty constructor for a mesh
    MultiGeometry() = delete;

    std::vector<Geometry::shared_ptr> d_geom;
};

} // namespace Geometry
} // namespace AMP

#endif
