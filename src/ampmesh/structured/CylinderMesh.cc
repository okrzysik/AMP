#include "ampmesh/structured/CylinderMesh.h"

#include "ampmesh/structured/BoxMesh.h"
#include "ampmesh/structured/BoxMeshHelpers.h"

#include "ampmesh/MultiIterator.h"
#include "ampmesh/shapes/Box.h"
#include "ampmesh/structured/structuredMeshElement.h"
#include "ampmesh/structured/structuredMeshIterator.h"

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
CylinderMesh::CylinderMesh( MeshParameters::shared_ptr params ) : BoxMesh( params )
{
    // Input options from the database
    PhysicalDim = d_db->getInteger( "dim" );
    GeomDim     = (GeomType) PhysicalDim;
    auto size   = d_db->getIntegerArray( "Size" );
    auto range  = d_db->getDoubleArray( "Range" );
    d_max_gcw   = d_db->getIntegerWithDefault( "GCW", 2 );
    std::vector<unsigned char> per( 1, false );
    if ( d_db->keyExists( "Periodic" ) )
        per = d_db->getBoolArray( "Periodic" );
    AMP_INSIST( per.size() == 1u, "Periodic must be an array of length 1" );
    AMP_INSIST( size.size() == 2u, "Size must be an array of length 2" );
    AMP_INSIST( range.size() == 3u, "Range must be an array of length 3" );
    AMP_INSIST( (int) PhysicalDim == 3, "dim must be 3" );
    d_range[0]      = range[0];
    d_range[1]      = range[1];
    d_range[2]      = range[2];
    d_globalSize[0] = 2 * size[0];
    d_globalSize[1] = 2 * size[0];
    d_globalSize[2] = size[1];
    d_offset.fill( 0 );
    AMP_ASSERT( d_range[0] > 0 );
    AMP_ASSERT( d_range[2] > d_range[1] );
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
std::vector<size_t>
CylinderMesh::estimateLogicalMeshSize( const MeshParameters::shared_ptr &params )
{
    auto db               = params->getDatabase();
    std::vector<int> size = db->getIntegerArray( "Size" );
    AMP_ASSERT( size.size() == 2u );
    std::vector<size_t> size2( 3 );
    size2[0] = 2 * size[0];
    size2[1] = 2 * size[0];
    size2[2] = size[1];
    return size2;
}


/****************************************************************
 * Functions to displace the mesh                                *
 ****************************************************************/
int CylinderMesh::isMeshMovable() const { return 1; }
void CylinderMesh::displaceMesh( const std::vector<double> &x )
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
void CylinderMesh::displaceMesh( const AMP::LinearAlgebra::Vector::const_shared_ptr )
{
    AMP_ERROR( "displaceMesh (vector) violates CylinderMesh properties" );
}
#endif


/****************************************************************
 * Copy the mesh                                                 *
 ****************************************************************/
AMP::shared_ptr<Mesh> CylinderMesh::copy() const { return AMP::make_shared<CylinderMesh>( *this ); }


/****************************************************************
 * Return the coordinate                                         *
 ****************************************************************/
void CylinderMesh::coord( const MeshElementIndex &index, double *pos ) const
{
    int i      = index.index( 0 );
    int j      = index.index( 1 );
    int k      = index.index( 2 );
    double x   = static_cast<double>( i ) / static_cast<double>( d_globalSize[0] );
    double y   = static_cast<double>( j ) / static_cast<double>( d_globalSize[1] );
    double z   = static_cast<double>( k ) / static_cast<double>( d_globalSize[2] );
    auto point = BoxMeshHelpers::map_logical_circle( d_range[0], 2, x, y );
    pos[0]     = point.first + d_offset[0];
    pos[1]     = point.second + d_offset[1];
    pos[2]     = d_range[1] + z * ( d_range[2] - d_range[1] ) + d_offset[2];
}


/****************************************************************
 * Return the logical coordinates                                *
 ****************************************************************/
std::array<double, 3> CylinderMesh::physicalToLogical( const double * ) const
{
    AMP_ERROR( "physicalToLogical is not supported in CylinderMesh" );
    return std::array<double, 3>();
}


} // namespace Mesh
} // namespace AMP
