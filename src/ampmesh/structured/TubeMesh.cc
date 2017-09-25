#include "ampmesh/structured/TubeMesh.h"
#include "ampmesh/MultiIterator.h"
#include "ampmesh/shapes/Box.h"
#include "ampmesh/structured/BoxMesh.h"
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
TubeMesh::TubeMesh( MeshParameters::shared_ptr params ) : BoxMesh( params )
{
    for ( int d = 0; d < 3; d++ ) {
        d_globalSize[d] = 1;
        d_isPeriodic[d] = false;
        d_numBlocks[d]  = 0;
    }
    // Check for valid inputs
    AMP_INSIST( params.get(), "Params must not be null" );
    AMP_INSIST( d_comm != AMP_MPI( AMP_COMM_NULL ), "Communicator must be set" );
    AMP_INSIST( d_db.get(), "Database must exist" );
    // Input options from the database
    PhysicalDim = d_db->getInteger( "dim" );
    GeomDim     = (GeomType) PhysicalDim;
    auto size   = d_db->getIntegerArray( "Size" );
    auto range  = d_db->getDoubleArray( "Range" );
    d_max_gcw   = d_db->getIntegerWithDefault( "GCW", 2 );
    std::vector<unsigned char> per( 1, false );
    if ( d_db->keyExists( "Periodic" ) )
        per = d_db->getBoolArray( "Periodic" );
    AMP_INSIST( size.size() == 3u, "Size must be an array of length 3" );
    AMP_INSIST( per.size() == 1u, "Periodic must be an array of length 1" );
    AMP_INSIST( range.size() == 4u, "Range must be an array of length 4" );
    AMP_INSIST( PhysicalDim == 3, "dim must be 3" );
    for ( int i = 0; i < 3; i++ )
        d_globalSize[i] = size[i];
    for ( int i = 0; i < 4; i++ )
        d_range[i] = range[i];
    d_offset.fill( 0 );
    AMP_ASSERT( d_range[0] > 0 );
    AMP_ASSERT( d_range[1] > d_range[0] );
    AMP_ASSERT( d_range[3] > d_range[2] );
    // Change the surface ids to match the standard ids
    // 0 - 8: Inner surface
    // 1 - 4: Outer surface
    // 4 - 2: Bottom surface
    // 5 - 1: Top surface
    d_surfaceId[0] = 8;
    d_surfaceId[1] = 4;
    d_surfaceId[2] = -1;
    d_surfaceId[3] = -1;
    d_surfaceId[4] = 2;
    d_surfaceId[5] = 1;
    d_onSurface[0] = true;
    d_onSurface[1] = true;
    d_onSurface[2] = false;
    d_onSurface[3] = false;
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
std::vector<size_t> TubeMesh::estimateLogicalMeshSize( const MeshParameters::shared_ptr &params )
{
    auto db               = params->getDatabase();
    std::vector<int> size = db->getIntegerArray( "Size" );
    AMP_ASSERT( size.size() == 3u );
    std::vector<size_t> size2( size.size() );
    for ( size_t d = 0; d < size.size(); d++ )
        size2[d] = size[d];
    return size2;
}


/****************************************************************
 * Functions to displace the mesh                                *
 ****************************************************************/
int TubeMesh::isMeshMovable() const { return 1; }
void TubeMesh::displaceMesh( const std::vector<double> &x )
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
void TubeMesh::displaceMesh( const AMP::LinearAlgebra::Vector::const_shared_ptr )
{
    AMP_ERROR( "displaceMesh (vector) violates TubeMesh properties" );
}
#endif


/****************************************************************
 * Copy the mesh                                                 *
 ****************************************************************/
AMP::shared_ptr<Mesh> TubeMesh::copy() const
{
    return AMP::shared_ptr<TubeMesh>( new TubeMesh( *this ) );
}


/****************************************************************
 * Return the coordinate                                         *
 ****************************************************************/
void TubeMesh::coord( const MeshElementIndex &index, double *pos ) const
{
    const double pi = 3.141592653589793116;
    int i           = index.index( 0 );
    int j           = index.index( 1 );
    int k           = index.index( 2 );
    double x        = static_cast<double>( i ) / static_cast<double>( d_globalSize[0] );
    double y        = static_cast<double>( j ) / static_cast<double>( d_globalSize[1] );
    double z        = static_cast<double>( k ) / static_cast<double>( d_globalSize[2] );
    double r        = d_range[0] + x * ( d_range[1] - d_range[0] );
    double theta    = 2.0 * pi * y;
    pos[0]          = r * cos( theta ) + d_offset[0];
    pos[1]          = r * sin( theta ) + d_offset[1];
    pos[2]          = d_range[2] + z * ( d_range[3] - d_range[2] ) + d_offset[2];
}


/****************************************************************
 * Return the logical coordinates                                *
 ****************************************************************/
std::array<double, 3> TubeMesh::physicalToLogical( const double * ) const
{
    AMP_ERROR( "physicalToLogical is not supported in TubeMesh" );
    return std::array<double, 3>();
}


} // namespace Mesh
} // namespace AMP
