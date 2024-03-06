#include "AMP/mesh/structured/PureLogicalMesh.h"
#include "AMP/mesh/MeshParameters.h"


namespace AMP::Mesh {


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
    auto db = params->getDatabase();
    AMP_INSIST( db.get(), "Database must exist" );
    if ( db->keyExists( "commSize" ) ) {
        d_rank = db->getScalar<int>( "commRank" );
        d_size = db->getScalar<int>( "commSize" );
        d_comm = AMP_COMM_NULL;
    } else {
        AMP_INSIST( !d_comm.isNull(), "Communicator must be set" );
    }
    // Fill basic mesh information
    auto size = db->getVector<int>( "Size" );
    auto per  = db->getWithDefault<std::vector<bool>>( "Periodic",
                                                      std::vector<bool>( size.size(), false ) );
    AMP_INSIST( size.size() >= 1u && size.size() <= 3u, "bad value for Size" );
    AMP_ASSERT( per.size() == size.size() );
    PhysicalDim = size.size();
    GeomDim     = static_cast<AMP::Mesh::GeomType>( size.size() );
    d_max_gcw   = db->getWithDefault<int>( "GCW", 2 );
    AMP_ASSERT( PhysicalDim == db->getWithDefault<int>( "dim", PhysicalDim ) );
    for ( size_t d = 0; d < size.size(); d++ ) {
        d_globalSize[d] = size[d];
        d_isPeriodic[d] = per[d];
        if ( !d_isPeriodic[d] ) {
            d_surfaceId[2 * d + 0] = 2 * d + 0;
            d_surfaceId[2 * d + 1] = 2 * d + 1;
        } else {
            d_surfaceId[2 * d + 0] = -1;
            d_surfaceId[2 * d + 1] = -1;
        }
    }
    // Initialize the logical mesh
    BoxMesh::initialize( db->getWithDefault<std::vector<int>>( "LoadBalanceMinSize", {} ) );
    BoxMesh::finalize( db->getString( "MeshName" ), getDisplacement( db ) );
}
PureLogicalMesh::PureLogicalMesh( const PureLogicalMesh &mesh ) = default;


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


/****************************************************************
 * Write/Read restart data                                       *
 ****************************************************************/
void PureLogicalMesh::writeRestart( int64_t fid ) const { BoxMesh::writeRestart( fid ); }
PureLogicalMesh::PureLogicalMesh( int64_t fid, AMP::IO::RestartManager *manager )
    : BoxMesh( fid, manager )
{
}


} // namespace AMP::Mesh
