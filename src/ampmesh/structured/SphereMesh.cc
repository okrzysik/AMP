#include "AMP/ampmesh/structured/SphereMesh.h"

#include "AMP/ampmesh/structured/BoxMesh.h"
#include "AMP/ampmesh/structured/BoxMeshHelpers.h"

#include "AMP/ampmesh/MultiIterator.h"
#include "AMP/ampmesh/shapes/Box.h"
#include "AMP/ampmesh/structured/structuredMeshElement.h"
#include "AMP/ampmesh/structured/structuredMeshIterator.h"

#ifdef USE_AMP_VECTORS
#include "AMP/vectors/Variable.h"
#include "AMP/vectors/Vector.h"
#include "AMP/vectors/VectorBuilder.h"
#endif
#ifdef USE_AMP_DISCRETIZATION
#include "AMP/discretization/DOF_Manager.h"
#include "AMP/discretization/simpleDOF_Manager.h"
#endif

namespace AMP {
namespace Mesh {


/****************************************************************
 * Constructors                                                  *
 ****************************************************************/
SphereMesh::SphereMesh( MeshParameters::shared_ptr params ) : BoxMesh( params )
{
    // Input options from the database
    PhysicalDim = d_db->getInteger( "dim" );
    GeomDim     = (GeomType) PhysicalDim;
    auto size   = d_db->getIntegerArray( "Size" );
    auto range  = d_db->getDoubleArray( "Range" );
    d_max_gcw   = d_db->getIntegerWithDefault( "GCW", 2 );
    AMP_INSIST( size.size() == 1u, "Size must be an array of length 1" );
    AMP_INSIST( range.size() == 1u, "Range must be an array of length 3" );
    AMP_INSIST( (int) PhysicalDim == 3, "dim must be 3" );
    d_r             = range[0];
    d_globalSize[0] = 2 * size[0];
    d_globalSize[1] = 2 * size[0];
    d_globalSize[2] = 2 * size[0];
    d_offset.fill( 0 );
    // Change the surface ids to match the standard ids
    // 0,1,2,3 - 4: Outer surface
    // 4 - 2: Bottom surface
    // 5 - 1: Top surface
    d_surfaceId[0] = 4;
    d_surfaceId[1] = 4;
    d_surfaceId[2] = 4;
    d_surfaceId[3] = 4;
    d_surfaceId[4] = 2;
    d_surfaceId[5] = 1;
    d_onSurface[0] = true;
    d_onSurface[1] = true;
    d_onSurface[2] = true;
    d_onSurface[3] = true;
    d_onSurface[4] = true;
    d_onSurface[5] = true;
    // Initialize the logical mesh
    BoxMesh::initialize();
    // Set the geometry
    // d_geometry.reset( new Geometry::Box( range ) );
    // Finalize the logical mesh
    BoxMesh::finalize();
}


/****************************************************************
 * Estimate the mesh size                                        *
 ****************************************************************/
std::vector<size_t> SphereMesh::estimateLogicalMeshSize( const MeshParameters::shared_ptr &params )
{
    auto db               = params->getDatabase();
    std::vector<int> size = db->getIntegerArray( "Size" );
    AMP_ASSERT( size.size() == 1u );
    return std::vector<size_t>( 3, 2 * size[0] );
}


/****************************************************************
 * Functions to displace the mesh                                *
 ****************************************************************/
int SphereMesh::isMeshMovable() const { return 1; }
void SphereMesh::displaceMesh( const std::vector<double> &x )
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
void SphereMesh::displaceMesh( const AMP::LinearAlgebra::Vector::const_shared_ptr )
{
    AMP_ERROR( "displaceMesh (vector) violates SphereMesh properties" );
}
#endif


/****************************************************************
 * Copy the mesh                                                 *
 ****************************************************************/
AMP::shared_ptr<Mesh> SphereMesh::copy() const { return AMP::make_shared<SphereMesh>( *this ); }


/****************************************************************
 * Return the coordinate                                         *
 ****************************************************************/
void SphereMesh::coord( const MeshElementIndex &index, double *pos ) const
{
    int i      = index.index( 0 );
    int j      = index.index( 1 );
    int k      = index.index( 2 );
    double x   = static_cast<double>( i ) / static_cast<double>( d_globalSize[0] );
    double y   = static_cast<double>( j ) / static_cast<double>( d_globalSize[1] );
    double z   = static_cast<double>( k ) / static_cast<double>( d_globalSize[2] );
    auto point = BoxMeshHelpers::map_logical_sphere( d_r, x, y, z );
    pos[0]     = std::get<0>( point ) + d_offset[0];
    pos[1]     = std::get<1>( point ) + d_offset[1];
    pos[2]     = std::get<2>( point ) + d_offset[2];
}


/****************************************************************
 * Return the logical coordinates                                *
 ****************************************************************/
std::array<double, 3> SphereMesh::physicalToLogical( const double * ) const
{
    AMP_ERROR( "physicalToLogical is not supported in SphereMesh" );
    return std::array<double, 3>();
}


} // namespace Mesh
} // namespace AMP
