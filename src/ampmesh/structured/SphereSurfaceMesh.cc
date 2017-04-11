#include "ampmesh/structured/SphereSurfaceMesh.h"
#include "ampmesh/structured/BoxMesh.h"
#include "ampmesh/structured/BoxMeshHelpers.h"

#include "ampmesh/MultiIterator.h"
#include "ampmesh/structured/structuredMeshElement.h"
#include "ampmesh/structured/structuredMeshIterator.h"
#include "ampmesh/shapes/Box.h"

#ifdef USE_AMP_VECTORS
#include "vectors/Variable.h"
#include "vectors/Vector.h"
#include "vectors/VectorBuilder.h"
#endif
#ifdef USE_AMP_DISCRETIZATION
#include "discretization/DOF_Manager.h"
#include "discretization/simpleDOF_Manager.h"
#endif

namespace AMP {
namespace Mesh {


/****************************************************************
* Constructors                                                  *
****************************************************************/
SphereSurfaceMesh::SphereSurfaceMesh( MeshParameters::shared_ptr params ):
    BoxMesh( params )
{
    // Input options from the database
    PhysicalDim = 3;
    GeomDim     = GeomType::Face;
    auto size   = d_db->getIntegerArray( "Size" );
    auto range  = d_db->getDoubleArray( "Range" );
    d_max_gcw   = d_db->getIntegerWithDefault( "GCW", 2 );
    AMP_INSIST( size.size() == 1u, "Size must be an array of length 1" );
    AMP_INSIST( range.size() == 1u, "Range must be an array of length 1" );
    AMP_INSIST( (int) PhysicalDim == 3, "dim must be 3" );
    AMP_ASSERT( range[0] >= 0 );
    d_r = range[0];
    d_isPeriodic[0] = true;
    d_isPeriodic[1] = false;
    d_isPeriodic[2] = false;
    d_globalSize[0] = size[0];
    d_globalSize[1] = size[0]/2;
    d_globalSize[2] = 1;
    d_offset.fill( 0 );
    // Change the surface ids to match the standard ids
    // 0,1,2,3 - 4: Outer surface
    // 4 - 2: Bottom surface
    // 5 - 1: Top surface
    d_surfaceId[0] = -1;
    d_surfaceId[1] = -1;
    d_surfaceId[2] = -1;
    d_surfaceId[3] = -1;
    d_surfaceId[4] = -1;
    d_surfaceId[5] = -1;
    d_onSurface[0] = false;
    d_onSurface[1] = false;
    d_onSurface[2] = false;
    d_onSurface[3] = false;
    d_onSurface[4] = false;
    d_onSurface[5] = false;
    // Initialize the logical mesh
    BoxMesh::initialize();
    // Set the geometry
    //d_geometry.reset( new Geometry::Box( range ) );
    // Finalize the logical mesh
    BoxMesh::finalize();
}


/****************************************************************
* Estimate the mesh size                                        *
****************************************************************/
std::vector<size_t> SphereSurfaceMesh::estimateLogicalMeshSize( const MeshParameters::shared_ptr &params )
{
    auto db = params->getDatabase();
    std::vector<int> size = db->getIntegerArray( "Size" );
    AMP_ASSERT(size.size()==1u);
    std::vector<size_t> size2(2);
    size2[0] = size[0];
    size2[1] = size[0] / 2;
    return size2;
}


/****************************************************************
* Functions to displace the mesh                                *
****************************************************************/
int SphereSurfaceMesh::isMeshMovable( ) const
{
    return 1;
}
void SphereSurfaceMesh::displaceMesh( const std::vector<double> &x )
{
    AMP_ASSERT( x.size() == PhysicalDim );
    for ( int i = 0; i < PhysicalDim; i++ ) {
        d_offset[i] += x[i];
        d_box[2 * i + 0] += x[i];
        d_box[2 * i + 1] += x[i];
        d_box_local[2 * i + 0] += x[i];
        d_box_local[2 * i + 1] += x[i];
    }
    if ( d_geometry != nullptr )
        d_geometry->displaceMesh( x );
}
#ifdef USE_AMP_VECTORS
void SphereSurfaceMesh::displaceMesh( const AMP::LinearAlgebra::Vector::const_shared_ptr )
{
    AMP_ERROR( "displaceMesh (vector) violates SphereSurfaceMesh properties" );
}
#endif


/****************************************************************
* Copy the mesh                                                 *
****************************************************************/
AMP::shared_ptr<Mesh> SphereSurfaceMesh::copy() const
{
    return AMP::shared_ptr<SphereSurfaceMesh>( new SphereSurfaceMesh(*this) );
}


/****************************************************************
* Return the coordinate                                         *
****************************************************************/
void SphereSurfaceMesh::coord( const MeshElementIndex &index, double *pos ) const
{
    int i = index.index(0);
    int j = index.index(1);
    double x = static_cast<double>(i) / static_cast<double>(d_globalSize[0]);
    double y = static_cast<double>(j) / static_cast<double>(d_globalSize[1]);
    auto point = BoxMeshHelpers::map_logical_sphere_surface( d_r, x, y );
    pos[0] = std::get<0>( point ) + d_offset[0];
    pos[1] = std::get<1>( point ) + d_offset[1];
    pos[2] = std::get<2>( point ) + d_offset[2];
}


/****************************************************************
* Return the logical coordinates                                *
****************************************************************/
std::array<double,3> SphereSurfaceMesh::physicalToLogical( const double *x ) const
{
    const double x0 = x[0] - d_offset[0];
    const double y0 = x[1] - d_offset[1];
    const double z0 = x[2] - d_offset[2];
    const double r = sqrt( x0*x0 + y0*y0 + z0*z0 );
    auto point = BoxMeshHelpers::map_sphere_surface_logical( d_r, x0, y0, z0 );
    std::array<double,3> pos = { point.first, point.second, r/d_r-1 };
    return pos;
}


} // Mesh namespace
} // AMP namespace



