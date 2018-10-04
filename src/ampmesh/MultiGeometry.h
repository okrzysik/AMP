#ifndef included_AMP_MultiGeometry
#define included_AMP_MultiGeometry

#include "AMP/ampmesh/Geometry.h"


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
    MultiGeometry( const std::vector<Geometry::shared_ptr> &geom );


    //! Destructor
    virtual ~MultiGeometry() = default;


public: // Functions inherited from Geometry
    virtual uint8_t getDim() const override final;
    virtual double distance( const Point &pos, const Point &dir ) const override final;
    virtual bool inside( const Point &pos ) const override final;
    virtual int surface( const Point &x ) const override final;
    virtual Point surfaceNorm( const Point &x ) const override final;
    virtual Point logical( const Point &x ) const override final;
    virtual Point physical( const Point &x ) const override final;
    virtual void displaceMesh( const double *x ) override final;

private:
    //! Empty constructor for a mesh
    MultiGeometry() = delete;

    std::vector<Geometry::shared_ptr> d_geom;
};

} // namespace Geometry
} // namespace AMP

#endif
