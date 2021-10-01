#include "AMP/ampmesh/structured/PureLogicalMesh.h"
#include "AMP/ampmesh/MeshParameters.h"


namespace AMP {
namespace Mesh {


/****************************************************************
 * Constructors                                                 *
 ****************************************************************/
PureLogicalMesh::PureLogicalMesh( std::shared_ptr<const MeshParameters> params ) : BoxMesh( params )
{
    // Basic defaults
    d_globalSize.fill( 1 );
    d_isPeriodic.fill( false );
    d_numBlocks.fill( 1 );
    // Check for valid inputs
    AMP_INSIST( params.get(), "Params must not be null" );
    AMP_INSIST( d_db.get(), "Database must exist" );
    if ( d_db->keyExists( "commSize" ) ) {
        d_rank = d_db->getScalar<int>( "commRank" );
        d_size = d_db->getScalar<int>( "commSize" );
        d_comm = AMP_COMM_NULL;
    } else {
        AMP_INSIST( !d_comm.isNull(), "Communicator must be set" );
    }
    // Fill basic mesh information
    auto size = d_db->getVector<int>( "Size" );
    AMP_INSIST( size.size() >= 1u && size.size() <= 3u, "bad value for size" );
    PhysicalDim = size.size();
    GeomDim     = static_cast<AMP::Mesh::GeomType>( size.size() );
    d_max_gcw   = d_db->getWithDefault<int>( "GCW", 2 );
    AMP_ASSERT( PhysicalDim == d_db->getWithDefault<int>( "dim", PhysicalDim ) );
    for ( size_t d = 0; d < size.size(); d++ ) {
        d_globalSize[d] = size[d];
        d_isPeriodic[d] = false;
    }
    for ( size_t d = 0; d < 2 * size.size(); d++ )
        d_surfaceId[d] = d;
    // Initialize the logical mesh
    BoxMesh::initialize();
    BoxMesh::finalize();
}
PureLogicalMesh::PureLogicalMesh( const PureLogicalMesh &mesh ) : BoxMesh( mesh )
{
    PhysicalDim  = mesh.PhysicalDim;
    GeomDim      = mesh.GeomDim;
    d_max_gcw    = mesh.d_max_gcw;
    d_comm       = mesh.d_comm;
    d_name       = mesh.d_name;
    d_box        = mesh.d_box;
    d_box_local  = mesh.d_box_local;
    d_isPeriodic = mesh.d_isPeriodic;
    d_globalSize = mesh.d_globalSize;
    d_blockSize  = mesh.d_blockSize;
    d_numBlocks  = mesh.d_numBlocks;
    d_surfaceId  = mesh.d_surfaceId;
    for ( int d = 0; d < 4; d++ ) {
        for ( int i = 0; i < 6; i++ )
            d_globalSurfaceList[i][d] = mesh.d_globalSurfaceList[i][d];
    }
}


/********************************************************
 * Create domain info                                    *
 ********************************************************/
void PureLogicalMesh::createBoundingBox()
{
    // Fill the bounding box
    d_box       = std::vector<double>( 2 * PhysicalDim );
    d_box_local = std::vector<double>( 2 * PhysicalDim );
    auto local  = getLocalBlock( d_comm.getRank() );
    for ( int d = 0; d < PhysicalDim; d++ ) {
        d_box_local[2 * d + 0] = local[2 * d + 0];
        d_box_local[2 * d + 1] = local[2 * d + 1];
        d_box[2 * d + 0]       = 0;
        d_box[2 * d + 1]       = d_globalSize[d];
    }
}


/********************************************************
 * Return the class name                                 *
 ********************************************************/
std::string PureLogicalMesh::meshClass() const { return "PureLogicalMesh"; }


/****************************************************************
 * Basic functions                                               *
 ****************************************************************/
Mesh::Movable PureLogicalMesh::isMeshMovable() const { return Mesh::Movable::Fixed; }
uint64_t PureLogicalMesh::positionHash() const { return 0; }
void PureLogicalMesh::displaceMesh( const std::vector<double> & )
{
    AMP_ERROR( "displaceMesh is not supported for PureLogicalMesh" );
}
void PureLogicalMesh::displaceMesh( std::shared_ptr<const AMP::LinearAlgebra::Vector> )
{
    AMP_ERROR( "displaceMesh is not supported for PureLogicalMesh" );
}
AMP::Geometry::Point PureLogicalMesh::physicalToLogical( const AMP::Geometry::Point &x ) const
{
    return x;
}
void PureLogicalMesh::coord( const MeshElementIndex &index, double *pos ) const
{
    AMP_ASSERT( index.type() == AMP::Mesh::GeomType::Vertex );
    for ( int d = 0; d < PhysicalDim; d++ )
        pos[d] = index.index( d );
}
std::unique_ptr<Mesh> PureLogicalMesh::clone() const
{
    return std::make_unique<PureLogicalMesh>( *this );
}


/****************************************************************
 * Check if two meshes are equal                                 *
 ****************************************************************/
bool PureLogicalMesh::operator==( const Mesh &rhs ) const
{
    // Check base class variables
    if ( !BoxMesh::operator==( rhs ) )
        return false;
    // Check if we can cast to a PureLogicalMesh
    auto mesh = dynamic_cast<const PureLogicalMesh *>( &rhs );
    if ( !mesh )
        return false;
    return true;
}


} // namespace Mesh
} // namespace AMP
