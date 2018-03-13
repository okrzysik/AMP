#include "AMP/ampmesh/structured/CubeMesh.h"

#include "AMP/ampmesh/MultiIterator.h"
#include "AMP/ampmesh/shapes/Box.h"
#include "AMP/ampmesh/structured/BoxMesh.h"
#include "AMP/ampmesh/structured/structuredMeshElement.h"
#include "AMP/ampmesh/structured/structuredMeshIterator.h"

namespace AMP {
namespace Mesh {


/****************************************************************
 * Constructors                                                  *
 ****************************************************************/
CubeMesh::CubeMesh( MeshParameters::shared_ptr params ) : StructuredGeometryMesh( params )
{
    // Input options from the database
    PhysicalDim           = d_db->getInteger( "dim" );
    GeomDim               = (GeomType) PhysicalDim;
    std::vector<int> size = d_db->getIntegerArray( "Size" );
    d_max_gcw             = d_db->getIntegerWithDefault( "GCW", 2 );
    std::vector<unsigned char> per( PhysicalDim, false );
    if ( d_db->keyExists( "Periodic" ) )
        per = d_db->getBoolArray( "Periodic" );
    AMP_INSIST( size.size() == PhysicalDim, "Size of field 'Size' must match dim" );
    for ( int d = 0; d < PhysicalDim; d++ ) {
        AMP_INSIST( size[d] > 0, "All dimensions must have a size > 0" );
        d_globalSize[d] = size[d];
        d_isPeriodic[d] = per[d];
    }
    // Initialize the logical mesh
    BoxMesh::initialize();
    // Fill the coordinates and set the geometry
    if ( d_db->keyExists( "Range" ) ) {
        auto range = d_db->getDoubleArray( "Range" );
        AMP_INSIST( range.size() == 2 * PhysicalDim, "Range must be 2*dim for cube generator" );
        if ( PhysicalDim == 1 )
            d_geometry.reset( new Geometry::Box<1>( range ) );
        if ( PhysicalDim == 2 )
            d_geometry.reset( new Geometry::Box<2>( range ) );
        if ( PhysicalDim == 3 )
            d_geometry.reset( new Geometry::Box<3>( range ) );
    } else if ( d_db->keyExists( "x_grid" ) ) {
        std::vector<std::vector<double>> coord( PhysicalDim );
        for ( int d = 0; d < PhysicalDim; d++ ) {
            if ( d == 0 ) {
                coord[d] = d_db->getDoubleArray( "x_grid" );
            } else if ( d == 1 ) {
                AMP_INSIST( d_db->keyExists( "y_grid" ), "Field 'y_grid' must exist in database'" );
                coord[d] = d_db->getDoubleArray( "y_grid" );
            } else if ( d == 2 ) {
                AMP_INSIST( d_db->keyExists( "z_grid" ), "Field 'z_grid' must exist in database'" );
                coord[d] = d_db->getDoubleArray( "z_grid" );
            } else {
                AMP_ERROR( "Physical Dimensions > 3 are not supported yet" );
            }
            AMP_ASSERT( (int) coord[d].size() == d_globalSize[d] + 1 );
        }
        if ( PhysicalDim == 1 )
            d_geometry.reset( new Geometry::Grid<1>( coord ) );
        if ( PhysicalDim == 2 )
            d_geometry.reset( new Geometry::Grid<2>( coord ) );
        if ( PhysicalDim == 3 )
            d_geometry.reset( new Geometry::Grid<3>( coord ) );
    }
    // Finalize the logical mesh
    BoxMesh::finalize();
}


/****************************************************************
 * Estimate the mesh size                                        *
 ****************************************************************/
std::vector<size_t> CubeMesh::estimateLogicalMeshSize( const MeshParameters::shared_ptr &params )
{
    auto db               = params->getDatabase();
    int dim               = db->getInteger( "dim" );
    std::vector<int> size = db->getIntegerArray( "Size" );
    AMP_ASSERT( (int) size.size() == dim );
    std::vector<size_t> size2( size.size() );
    for ( size_t d = 0; d < size.size(); d++ )
        size2[d] = size[d];
    return size2;
}


/****************************************************************
 * Copy the mesh                                                 *
 ****************************************************************/
AMP::shared_ptr<Mesh> CubeMesh::clone() const { return AMP::make_shared<CubeMesh>( *this ); }


} // namespace Mesh
} // namespace AMP
