#include "ampmesh/structured/CubeMesh.h"
#include "ampmesh/structured/BoxMesh.h"
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
CubeMesh::CubeMesh( MeshParameters::shared_ptr params ):
    BoxMesh( params )
{
    // Input options from the database
    PhysicalDim           = d_db->getInteger( "dim" );
    GeomDim               = (GeomType) PhysicalDim;
    std::vector<int> size = d_db->getIntegerArray( "Size" );
    d_max_gcw             = d_db->getIntegerWithDefault( "GCW", 2 );
    std::vector<unsigned char> per(PhysicalDim,false);
    if ( d_db->keyExists( "Periodic" ) )
        per = d_db->getBoolArray( "Periodic" );
    AMP_INSIST( size.size() == PhysicalDim, "Size of field 'Size' must match dim" );
    for ( int d=0; d<PhysicalDim; d++ ) {
        AMP_INSIST( size[d] > 0, "All dimensions must have a size > 0" );
        d_globalSize[d] = size[d];
        d_isPeriodic[d] = per[d];
    }
    // Initialize the logical mesh
    BoxMesh::initialize();
    // Fill the coordinates
    std::vector<double> range( 2*PhysicalDim, 0 );
    if ( d_db->keyExists( "Range" ) ) {
        range = d_db->getDoubleArray( "Range" );
        AMP_INSIST( range.size() == 2 * PhysicalDim, "Range must be 2*dim for cube generator" );
        for (int d=0; d<PhysicalDim; d++) {
            d_coord[d].resize( d_globalSize[d]+1 );
            const double x0 = range[2*d];
            const double dx = ( range[2*d+1] - range[2*d+0] ) / d_globalSize[d];
            for (int i=0; i<=d_globalSize[d]; i++)
                d_coord[d][i] = x0 + i*dx;
        }
    } else if ( d_db->keyExists( "x_grid" ) ) {
        for ( int d = 0; d < PhysicalDim; d++ ) {
            if ( d == 0 ) {
                d_coord[d] = d_db->getDoubleArray( "x_grid" );
            } else if ( d == 1 ) {
                AMP_INSIST( d_db->keyExists( "y_grid" ),
                            "Field 'y_grid' must exist in database'" );
                d_coord[d] = d_db->getDoubleArray( "y_grid" );
            } else if ( d == 2 ) {
                AMP_INSIST( d_db->keyExists( "z_grid" ),
                            "Field 'z_grid' must exist in database'" );
                d_coord[d] = d_db->getDoubleArray( "z_grid" );
            } else {
                AMP_ERROR( "Physical Dimensions > 3 are not supported yet" );
            }
            AMP_ASSERT((int)d_coord[d].size()==d_globalSize[d]+1);
            range[2*d+0] = d_coord[d].front();
            range[2*d+1] = d_coord[d].back();
        }
    }
    // Set the geometry
    d_geometry.reset( new Geometry::Box( range ) );
    // Finalize the logical mesh
    BoxMesh::finalize();
}


/****************************************************************
* Estimate the mesh size                                        *
****************************************************************/
std::vector<size_t> CubeMesh::estimateLogicalMeshSize( const MeshParameters::shared_ptr &params )
{
    auto db = params->getDatabase();
    int dim = db->getInteger( "dim" );
    std::vector<int> size = db->getIntegerArray( "Size" );
    AMP_ASSERT((int)size.size()==dim);
    std::vector<size_t> size2(size.size());
    for ( size_t d=0; d<size.size(); d++)
        size2[d] = size[d];
    return size2;
}


/****************************************************************
* Functions to displace the mesh                                *
****************************************************************/
int CubeMesh::isMeshMovable( ) const
{
    return 1;
}void CubeMesh::displaceMesh( const std::vector<double> &x )
{
    AMP_ASSERT( x.size() == PhysicalDim );
    for ( int i = 0; i < PhysicalDim; i++ ) {
        for ( size_t j = 0; j < d_coord[i].size(); j++ )
            d_coord[i][j] += x[i];
        d_box[2 * i + 0] += x[i];
        d_box[2 * i + 1] += x[i];
        d_box_local[2 * i + 0] += x[i];
        d_box_local[2 * i + 1] += x[i];
    }
    if ( d_geometry != nullptr )
        d_geometry->displaceMesh( x );
}
#ifdef USE_AMP_VECTORS
void CubeMesh::displaceMesh( const AMP::LinearAlgebra::Vector::const_shared_ptr )
{
    AMP_ERROR( "displaceMesh (vector) violates CubeMesh properties" );
}
#endif


/****************************************************************
* Copy the mesh                                                 *
****************************************************************/
AMP::shared_ptr<Mesh> CubeMesh::copy() const
{
    return AMP::shared_ptr<CubeMesh>( new CubeMesh(*this) );
}


/****************************************************************
* Return the coordinate                                         *
****************************************************************/
void CubeMesh::coord( const MeshElementIndex &index, double *pos ) const
{
    AMP_ASSERT( index.type() == AMP::Mesh::Vertex );
    for ( int d  = 0; d < PhysicalDim; d++ ) {
        int i = index.index(d);
        if ( i>=0 && i <=d_globalSize[d] ) {
            pos[d] = d_coord[d][i];
        } else {
            int shift = 0;
            if ( i < 0 ) {
                i += d_globalSize[d];
                shift = -1;
            } else if ( i > d_globalSize[d] ) {
                i -= d_globalSize[d];
                shift = 1;
            }
            pos[d] = d_coord[d][i] + shift*(d_box[2*d+1]-d_box[2*d]);
        }
    }
}


} // Mesh namespace
} // AMP namespace



