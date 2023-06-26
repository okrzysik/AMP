#include "AMP/mesh/structured/StructuredGeometryMesh.h"
#include "AMP/IO/HDF5.hpp"
#include "AMP/IO/RestartManager.h"
#include "AMP/mesh/MeshParameters.h"


namespace AMP::Mesh {


/****************************************************************
 * Constructors                                                 *
 ****************************************************************/
StructuredGeometryMesh::StructuredGeometryMesh( std::shared_ptr<const MeshParameters> params )
    : BoxMesh( params ), d_pos_hash( 0 )
{
    // Basic defaults
    d_globalSize.fill( 1 );
    d_isPeriodic.fill( false );
    d_numBlocks.fill( 1 );
    // Check for valid inputs
    AMP_INSIST( params.get(), "Params must not be null" );
    auto db = params->getDatabase();
    AMP_INSIST( db.get(), "Database must exist" );
    AMP_INSIST( d_comm != AMP_MPI( AMP_COMM_NULL ), "Communicator must be set" );
    // Construct the geometry
    auto db2 = db->cloneDatabase();
    db2->erase( "x_offset", false );
    db2->erase( "y_offset", false );
    db2->erase( "z_offset", false );
    d_geometry  = AMP::Geometry::Geometry::buildGeometry( std::move( db2 ) );
    d_geometry2 = std::dynamic_pointer_cast<AMP::Geometry::LogicalGeometry>( d_geometry );
    AMP_ASSERT( d_geometry2 );
    // Fill basic mesh information
    PhysicalDim = d_geometry2->getDim();
    GeomDim     = static_cast<AMP::Mesh::GeomType>( d_geometry2->getLogicalDim() );
    d_max_gcw   = db->getWithDefault<int>( "GCW", 2 );
    AMP_ASSERT( PhysicalDim == db->getWithDefault<int>( "dim", PhysicalDim ) );
    auto size = d_geometry2->getLogicalGridSize( db->getVector<int>( "Size" ) );
    AMP_ASSERT( size.size() == static_cast<size_t>( GeomDim ) );
    for ( size_t d = 0; d < size.size(); d++ )
        d_globalSize[d] = size[d];
    auto isPeriodic = d_geometry2->getPeriodicDim();
    for ( size_t d = 0; d < isPeriodic.size(); d++ )
        d_isPeriodic[d] = isPeriodic[d];
    auto surfaceIds = d_geometry2->getLogicalSurfaceIds();
    for ( size_t d = 0; d < surfaceIds.size(); d++ )
        d_surfaceId[d] = surfaceIds[d];
    // Initialize the logical mesh
    BoxMesh::initialize( db->getWithDefault<std::vector<int>>( "LoadBalanceMinSize", {} ) );
    BoxMesh::finalize( db->getString( "MeshName" ), getDisplacement( db ) );
}
StructuredGeometryMesh::StructuredGeometryMesh( const StructuredGeometryMesh &mesh )
    : BoxMesh( mesh )
{
    d_geometry2 = std::dynamic_pointer_cast<AMP::Geometry::LogicalGeometry>( d_geometry );
    d_pos_hash  = mesh.d_pos_hash;
}


/****************************************************************
 * Write/Read restart data                                       *
 ****************************************************************/
void StructuredGeometryMesh::writeRestart( int64_t fid ) const
{
    writeHDF5( fid, "MeshType", std::string( "StructuredGeometryMesh" ) );
    writeHDF5( fid, "MeshName", d_name );
    writeHDF5( fid, "MeshID", d_meshID );
    writeHDF5( fid, "comm", d_comm.hashRanks() );
    writeHDF5( fid, "geometry", d_geometry );
    writeHDF5( fid, "gcw", d_max_gcw );
    writeHDF5( fid, "size", d_globalSize );
    writeHDF5( fid, "periodic", d_isPeriodic );
    writeHDF5( fid, "surfaceId", d_surfaceId );
    writeHDF5( fid, "numBlocks", d_numBlocks );
    writeHDF5( fid, "startIndex[0]", d_startIndex[0] );
    writeHDF5( fid, "startIndex[1]", d_startIndex[1] );
    writeHDF5( fid, "startIndex[2]", d_startIndex[2] );
    writeHDF5( fid, "endIndex[0]", d_endIndex[0] );
    writeHDF5( fid, "endIndex[1]", d_endIndex[1] );
    writeHDF5( fid, "endIndex[2]", d_endIndex[2] );
    writeHDF5( fid, "d_pos_hash", d_pos_hash );
}
StructuredGeometryMesh::StructuredGeometryMesh( int64_t fid, AMP::IO::RestartManager *manager )
    : BoxMesh()
{
    // Set the data for Mesh
    readHDF5( fid, "MeshName", d_name );
    uint64_t commHash = 0;
    readHDF5( fid, "comm", commHash );
    d_comm = manager->getComm( commHash );
    // Basic defaults
    d_globalSize.fill( 1 );
    d_isPeriodic.fill( false );
    d_numBlocks.fill( 1 );
    // Load the geometry
    readHDF5( fid, "geometry", d_geometry );
    d_geometry2 = std::dynamic_pointer_cast<AMP::Geometry::LogicalGeometry>( d_geometry );
    AMP_ASSERT( d_geometry2 );
    // Fill basic mesh information
    // Note: we do not call initialize to avoid any changes to the parallel decomposition
    PhysicalDim = d_geometry2->getDim();
    GeomDim     = static_cast<AMP::Mesh::GeomType>( d_geometry2->getLogicalDim() );
    readHDF5( fid, "gcw", d_max_gcw );
    readHDF5( fid, "size", d_globalSize );
    readHDF5( fid, "periodic", d_isPeriodic );
    readHDF5( fid, "surfaceId", d_surfaceId );
    readHDF5( fid, "numBlocks", d_numBlocks );
    readHDF5( fid, "startIndex[0]", d_startIndex[0] );
    readHDF5( fid, "startIndex[1]", d_startIndex[1] );
    readHDF5( fid, "startIndex[2]", d_startIndex[2] );
    readHDF5( fid, "endIndex[0]", d_endIndex[0] );
    readHDF5( fid, "endIndex[1]", d_endIndex[1] );
    readHDF5( fid, "endIndex[2]", d_endIndex[2] );
    readHDF5( fid, "d_pos_hash", d_pos_hash );
    auto block   = getLocalBlock( d_rank );
    d_indexSize  = { block[1] - block[0] + 3, block[3] - block[2] + 3, block[5] - block[4] + 3 };
    d_localIndex = block;
    for ( int d = 0; d < 3; d++ ) {
        if ( d_localIndex[2 * d + 1] == d_globalSize[d] - 1 )
            d_localIndex[2 * d + 1] = d_globalSize[d];
        d_localIndex[2 * d + 1]++;
    }
    BoxMesh::finalize( d_name, {} );
}


/********************************************************
 * Return the class name                                 *
 ********************************************************/
std::string StructuredGeometryMesh::meshClass() const { return "StructuredGeometryMesh"; }


/****************************************************************
 * Basic functions                                               *
 ****************************************************************/
Mesh::Movable StructuredGeometryMesh::isMeshMovable() const { return Mesh::Movable::Displace; }
uint64_t StructuredGeometryMesh::positionHash() const { return d_pos_hash; }
void StructuredGeometryMesh::displaceMesh( const std::vector<double> &x )
{
    for ( int i = 0; i < PhysicalDim; i++ ) {
        d_box[2 * i + 0] += x[i];
        d_box[2 * i + 1] += x[i];
        d_box_local[2 * i + 0] += x[i];
        d_box_local[2 * i + 1] += x[i];
    }
    d_geometry2->displace( x.data() );
    d_pos_hash++;
}
void StructuredGeometryMesh::displaceMesh( std::shared_ptr<const AMP::LinearAlgebra::Vector> )
{
    AMP_ERROR( "displaceMesh (vector) violates StructuredGeometryMesh properties" );
}
AMP::Geometry::Point
StructuredGeometryMesh::physicalToLogical( const AMP::Geometry::Point &x ) const
{
    return d_geometry2->logical( x );
}
void StructuredGeometryMesh::coord( const MeshElementIndex &index, double *pos ) const
{
    AMP_ASSERT( index.type() == AMP::Mesh::GeomType::Vertex );
    double x = static_cast<double>( index.index( 0 ) ) / static_cast<double>( d_globalSize[0] );
    double y = static_cast<double>( index.index( 1 ) ) / static_cast<double>( d_globalSize[1] );
    double z = static_cast<double>( index.index( 2 ) ) / static_cast<double>( d_globalSize[2] );
    auto tmp = d_geometry2->physical( AMP::Geometry::Point( x, y, z ) );
    for ( int d = 0; d < PhysicalDim; d++ )
        pos[d] = tmp[d];
}
std::unique_ptr<Mesh> StructuredGeometryMesh::clone() const
{
    return std::make_unique<StructuredGeometryMesh>( *this );
}


/****************************************************************
 * Check if two meshes are equal                                 *
 ****************************************************************/
bool StructuredGeometryMesh::operator==( const Mesh &rhs ) const
{
    // Check base class variables
    if ( !BoxMesh::operator==( rhs ) )
        return false;
    // Check if we can cast to a MovableBoxMesh
    auto mesh = dynamic_cast<const StructuredGeometryMesh *>( &rhs );
    if ( !mesh )
        return false;
    // Perform basic comparison
    return d_geometry2 == mesh->d_geometry2;
}


} // namespace AMP::Mesh
