#include "ampmesh/structured/ShellMesh.h"
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
ShellMesh::ShellMesh( MeshParameters::shared_ptr params ):
    BoxMesh( params )
{
    // Input options from the database
    PhysicalDim = d_db->getInteger( "dim" );
    GeomDim     = (GeomType) PhysicalDim;
    auto size   = d_db->getIntegerArray( "Size" );
    auto range  = d_db->getDoubleArray( "Range" );
    d_max_gcw   = d_db->getIntegerWithDefault( "GCW", 2 );
    AMP_INSIST( size.size() == 2u, "Size must be an array of length 2" );
    AMP_INSIST( range.size() == 2u, "Range must be an array of length 2" );
    AMP_INSIST( (int) PhysicalDim == 3, "dim must be 3" );
    AMP_ASSERT( range[0] >= 0 && range[1] > 0 && (range[1]-range[0]) > 0 );
    d_r_min = range[0];
    d_r_max = range[1];
    d_isPeriodic[0] = true;
    d_isPeriodic[1] = false;
    d_isPeriodic[2] = false;
    d_globalSize[0] = size[1];
    d_globalSize[1] = size[1] / 2;
    d_globalSize[2] = size[0];
    d_offset.fill( 0 );
    // Change the surface ids to match the standard ids
    // 0,1,2,3 - 4: Outer surface
    // 4 - 2: Bottom surface
    // 5 - 1: Top surface
    d_surfaceId[0] = -1;
    d_surfaceId[1] = -1;
    d_surfaceId[2] = 1;
    d_surfaceId[3] = 2;
    d_surfaceId[4] = 3;
    d_surfaceId[5] = 4;
    d_onSurface[0] = false;
    d_onSurface[1] = false;
    d_onSurface[2] = true;
    d_onSurface[3] = true;
    d_onSurface[4] = true;
    d_onSurface[5] = true;
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
std::vector<size_t> ShellMesh::estimateLogicalMeshSize( const MeshParameters::shared_ptr &params )
{
    auto db = params->getDatabase();
    std::vector<int> size = db->getIntegerArray( "Size" );
    AMP_ASSERT(size.size()==2u);
    std::vector<size_t> size2(3);
    size2[0] = size[1];
    size2[1] = size[1] / 2;
    size2[2] = size[0];
    return size2;
}


/****************************************************************
* Functions to displace the mesh                                *
****************************************************************/
int ShellMesh::isMeshMovable( ) const
{
    return 1;
}
void ShellMesh::displaceMesh( const std::vector<double> &x )
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
void ShellMesh::displaceMesh( const AMP::LinearAlgebra::Vector::const_shared_ptr )
{
    AMP_ERROR( "displaceMesh (vector) violates ShellMesh properties" );
}
#endif


/****************************************************************
* Copy the mesh                                                 *
****************************************************************/
AMP::shared_ptr<Mesh> ShellMesh::copy() const
{
    return AMP::shared_ptr<ShellMesh>( new ShellMesh(*this) );
}


/****************************************************************
* Return the coordinate                                         *
****************************************************************/
void ShellMesh::coord( const MeshElementIndex &index, double *pos ) const
{
    int i = index.index(0);
    int j = index.index(1);
    int k = index.index(2);
    double x = static_cast<double>(i) / static_cast<double>(d_globalSize[0]);
    double y = static_cast<double>(j) / static_cast<double>(d_globalSize[1]);
    double z = static_cast<double>(k) / static_cast<double>(d_globalSize[2]);
    auto point = BoxMeshHelpers::map_logical_shell( d_r_min, d_r_max, x, y, z );
    pos[0] = std::get<0>( point ) + d_offset[0];
    pos[1] = std::get<1>( point ) + d_offset[1];
    pos[2] = std::get<2>( point ) + d_offset[2];
}


} // Mesh namespace
} // AMP namespace



