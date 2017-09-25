#include "ampmesh/structured/CircleMesh.h"
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
CircleMesh::CircleMesh( MeshParameters::shared_ptr params ) : BoxMesh( params )
{
    // Check for valid inputs
    AMP_INSIST( params.get(), "Params must not be null" );
    AMP_INSIST( d_comm != AMP_MPI( AMP_COMM_NULL ), "Communicator must be set" );
    AMP_INSIST( d_db.get(), "Database must exist" );
    // Input options from the database
    PhysicalDim           = d_db->getInteger( "dim" );
    GeomDim               = (GeomType) PhysicalDim;
    std::vector<int> size = d_db->getIntegerArray( "Size" );
    auto range            = d_db->getDoubleArray( "Range" );
    d_max_gcw             = d_db->getIntegerWithDefault( "GCW", 2 );
    AMP_INSIST( (int) PhysicalDim == 2, "dim must be size 2" );
    AMP_INSIST( size.size() == 1u, "Size must be an array of length 1" );
    AMP_INSIST( range.size() == 1u, "Size must be an array of length 1" );
    d_globalSize[0] = 2 * size[0];
    d_globalSize[1] = 2 * size[0];
    d_isPeriodic[0] = false;
    d_isPeriodic[1] = false;
    d_R             = range[0];
    d_offset.fill( 0 );
    // Change the surface ids to match the standard ids
    d_surfaceId[0] = 1;
    d_surfaceId[1] = 1;
    d_surfaceId[2] = 1;
    d_surfaceId[3] = 1;
    d_onSurface[0] = true;
    d_onSurface[1] = true;
    d_onSurface[2] = true;
    d_onSurface[3] = true;
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
std::vector<size_t> CircleMesh::estimateLogicalMeshSize( const MeshParameters::shared_ptr &params )
{
    auto db               = params->getDatabase();
    std::vector<int> size = db->getIntegerArray( "Size" );
    AMP_ASSERT( size.size() == 1u );
    std::vector<size_t> size2( 1, 2 * size[0] );
    return size2;
}


/****************************************************************
* Functions to displace the mesh                                *
****************************************************************/
int CircleMesh::isMeshMovable() const { return 1; }
void CircleMesh::displaceMesh( const std::vector<double> &x )
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
void CircleMesh::displaceMesh( const AMP::LinearAlgebra::Vector::const_shared_ptr )
{
    AMP_ERROR( "displaceMesh (vector) violates CircleMesh properties" );
}
#endif


/****************************************************************
* Copy the mesh                                                 *
****************************************************************/
AMP::shared_ptr<Mesh> CircleMesh::copy() const
{
    return AMP::shared_ptr<CircleMesh>( new CircleMesh( *this ) );
}


/****************************************************************
* Return the coordinate                                         *
****************************************************************/
void CircleMesh::coord( const MeshElementIndex &index, double *pos ) const
{
    // This maps from a a logically rectangular 2D mesh to a circular mesh using the mapping by:
    // Dona Calhoun, Christiane Helzel, Randall LeVeque, "Logically Rectangular Grids and Finite
    // GeomType::Volume
    //    Methods for PDEs in Circular and Spherical Domains", SIAM REVIEW, Vol. 50, No. 4, pp.
    //    723–752 (2008)
    int i      = index.index( 0 );
    int j      = index.index( 1 );
    double x   = static_cast<double>( i ) / static_cast<double>( d_globalSize[0] );
    double y   = static_cast<double>( j ) / static_cast<double>( d_globalSize[1] );
    auto point = BoxMeshHelpers::map_logical_circle( d_R, 2, x, y );
    pos[0]     = point.first + d_offset[0];
    pos[1]     = point.second + d_offset[1];
}


/****************************************************************
* Return the logical coordinates                                *
****************************************************************/
std::array<double, 3> CircleMesh::physicalToLogical( const double *x ) const
{
    auto point = BoxMeshHelpers::map_circle_logical( d_R, 2, x[0], x[1] );
    std::array<double, 3> y = { point.first, point.second, 0 };
    return y;
}


} // Mesh namespace
} // AMP namespace
