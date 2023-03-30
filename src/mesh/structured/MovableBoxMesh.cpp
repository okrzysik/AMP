#include "AMP/mesh/structured/MovableBoxMesh.h"
#include "AMP/IO/RestartManager.h"
#include "AMP/discretization/DOF_Manager.h"
#include "AMP/discretization/simpleDOF_Manager.h"
#include "AMP/mesh/MultiIterator.h"
#include "AMP/mesh/structured/BoxMesh.h"
#include "AMP/mesh/structured/structuredMeshElement.h"
#include "AMP/mesh/structured/structuredMeshIterator.h"
#include "AMP/utils/Utilities.hpp"
#include "AMP/vectors/Variable.h"
#include "AMP/vectors/Vector.h"
#include "AMP/vectors/VectorBuilder.h"


namespace AMP::Mesh {


/****************************************************************
 * Constructors                                                  *
 ****************************************************************/
MovableBoxMesh::MovableBoxMesh( const AMP::Mesh::BoxMesh &mesh ) : BoxMesh( mesh ), d_pos_hash( 0 )
{
    if ( dynamic_cast<const MovableBoxMesh *>( &mesh ) ) {
        // We are copying another MovableBoxMesh
        auto rhs = dynamic_cast<const MovableBoxMesh &>( mesh );
        d_index  = rhs.d_index;
        d_coord  = rhs.d_coord;
    } else {
        // Get a list of all nodes on the current processor
        MeshIterator nodeIterator = mesh.getIterator( GeomType::Vertex, d_max_gcw );
        d_index.reserve( nodeIterator.size() );
        for ( size_t i = 0; i < nodeIterator.size(); ++i, ++nodeIterator ) {
            auto element = dynamic_cast<structuredMeshElement *>( nodeIterator->getRawElement() );
            AMP_ASSERT( element != nullptr );
            d_index.emplace_back( element->getIndex() );
        }
        AMP::Utilities::quicksort( d_index );
        // Generate coordinates
        d_coord.resize( PhysicalDim, d_index.size() );
        d_coord.fill( 0 );
        for ( size_t i = 0; i < d_index.size(); i++ )
            mesh.coord( d_index[i], &d_coord( 0, i ) );
    }
}


/****************************************************************
 * Write/Read restart data                                       *
 ****************************************************************/
void MovableBoxMesh::writeRestart( int64_t fid ) const
{
    writeHDF5( fid, "MeshType", std::string( "MovableBoxMesh" ) );
    // Write BoxMesh data
    writeHDF5( fid, "MeshName", d_name );
    writeHDF5( fid, "MeshID", d_meshID );
    writeHDF5( fid, "comm", d_comm.hashRanks() );
    writeHDF5( fid, "PhysicalDim", PhysicalDim );
    writeHDF5( fid, "GeomDim", GeomDim );
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
    // Write MovableBoxMesh data
    writeHDF5( fid, "pos_hash", d_pos_hash );
    writeHDF5( fid, "index", d_index );
    writeHDF5( fid, "coord", d_coord );
    writeHDF5( fid, "ids", d_ids );
}
MovableBoxMesh::MovableBoxMesh( int64_t fid, AMP::IO::RestartManager *manager ) : BoxMesh()
{
    // Set the data for Mesh
    readHDF5( fid, "MeshName", d_name );
    uint64_t commHash;
    readHDF5( fid, "comm", commHash );
    d_comm = manager->getComm( commHash );
    // Basic defaults
    d_globalSize.fill( 1 );
    d_isPeriodic.fill( false );
    d_numBlocks.fill( 1 );
    // Fill basic mesh information
    // Note: we do not call initialize to avoid any changes to the parallel decomposition
    readHDF5( fid, "PhysicalDim", PhysicalDim );
    readHDF5( fid, "GeomDim", GeomDim );
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
    auto block   = getLocalBlock( d_rank );
    d_indexSize  = { block[1] - block[0] + 3, block[3] - block[2] + 3, block[5] - block[4] + 3 };
    d_localIndex = block;
    for ( int d = 0; d < 3; d++ ) {
        if ( d_localIndex[2 * d + 1] == d_globalSize[d] - 1 )
            d_localIndex[2 * d + 1] = d_globalSize[d];
        d_localIndex[2 * d + 1]++;
    }
    // Read MovableBoxMesh data
    readHDF5( fid, "pos_hash", d_pos_hash );
    readHDF5( fid, "index", d_index );
    readHDF5( fid, "coord", d_coord );
    readHDF5( fid, "ids", d_ids );
    BoxMesh::finalize( d_name, {} );
}


/****************************************************************
 * Functions to displace the mesh                                *
 ****************************************************************/
Mesh::Movable MovableBoxMesh::isMeshMovable() const { return Mesh::Movable::Deform; }
uint64_t MovableBoxMesh::positionHash() const { return d_pos_hash; }
void MovableBoxMesh::displaceMesh( const std::vector<double> &x )
{
    AMP_ASSERT( x.size() == PhysicalDim );
    for ( size_t i = 0; i < d_coord.size( 1 ); i++ ) {
        for ( int d = 0; d < PhysicalDim; d++ )
            d_coord( d, i ) += x[d];
    }
    for ( int d = 0; d < PhysicalDim; d++ ) {
        d_box[2 * d + 0] += x[d];
        d_box[2 * d + 1] += x[d];
        d_box_local[2 * d + 0] += x[d];
        d_box_local[2 * d + 1] += x[d];
    }
    if ( d_geometry != nullptr )
        d_geometry->displace( x.data() );
    d_pos_hash++;
}
void MovableBoxMesh::displaceMesh( const AMP::LinearAlgebra::Vector::const_shared_ptr x )
{
    // Clear the geometry if it exists to ensure consistency
    d_geometry.reset();
    // Create the position vector with the necessary ghost nodes
    auto DOFs = AMP::Discretization::simpleDOFManager::create(
        shared_from_this(),
        getIterator( AMP::Mesh::GeomType::Vertex, d_max_gcw ),
        getIterator( AMP::Mesh::GeomType::Vertex, 0 ),
        PhysicalDim );
    auto nodalVariable = std::make_shared<AMP::LinearAlgebra::Variable>( "tmp_pos" );
    auto displacement  = AMP::LinearAlgebra::createVector( DOFs, nodalVariable, false );
    std::vector<size_t> dofs1( PhysicalDim );
    std::vector<size_t> dofs2( PhysicalDim );
    auto cur  = getIterator( AMP::Mesh::GeomType::Vertex, 0 );
    auto end  = cur.end();
    auto DOFx = x->getDOFManager();
    std::vector<double> data( PhysicalDim );
    while ( cur != end ) {
        AMP::Mesh::MeshElementID id = cur->globalID();
        DOFx->getDOFs( id, dofs1 );
        DOFs->getDOFs( id, dofs2 );
        x->getValuesByGlobalID( PhysicalDim, &dofs1[0], &data[0] );
        displacement->setValuesByGlobalID( PhysicalDim, &dofs2[0], &data[0] );
        ++cur;
    }
    displacement->makeConsistent( AMP::LinearAlgebra::VectorData::ScatterType::CONSISTENT_SET );
    // Move all nodes (including the ghost nodes)
    std::vector<size_t> dofs( PhysicalDim );
    std::vector<double> disp( PhysicalDim );
    for ( size_t i = 0; i < d_coord.size( 1 ); i++ ) {
        MeshElementID id = structuredMeshElement( d_index[i], this ).globalID();
        DOFs->getDOFs( id, dofs );
        AMP_ASSERT( dofs.size() == PhysicalDim );
        displacement->getValuesByGlobalID( (int) PhysicalDim, &dofs[0], &disp[0] );
        for ( int d = 0; d < PhysicalDim; d++ )
            d_coord( d, i ) += disp[d];
    }
    // Compute the new bounding box of the mesh
    d_box_local = std::vector<double>( 2 * PhysicalDim );
    for ( int d = 0; d < PhysicalDim; d++ ) {
        d_box_local[2 * d + 0] = 1e100;
        d_box_local[2 * d + 1] = -1e100;
    }
    for ( size_t i = 0; i < d_coord.size( 1 ); i++ ) {
        for ( int d = 0; d < PhysicalDim; d++ ) {
            d_box_local[2 * d + 0] = std::min( d_box_local[2 * d + 0], d_coord( d, i ) );
            d_box_local[2 * d + 1] = std::max( d_box_local[2 * d + 1], d_coord( d, i ) );
        }
    }
    d_box = std::vector<double>( PhysicalDim * 2 );
    for ( int i = 0; i < PhysicalDim; i++ ) {
        d_box[2 * i + 0] = d_comm.minReduce( d_box_local[2 * i + 0] );
        d_box[2 * i + 1] = d_comm.maxReduce( d_box_local[2 * i + 1] );
    }
    d_pos_hash++;
}


/********************************************************
 * Return the class name                                 *
 ********************************************************/
std::string MovableBoxMesh::meshClass() const { return "MovableBoxMesh"; }


/****************************************************************
 * Copy the mesh                                                 *
 ****************************************************************/
std::unique_ptr<Mesh> MovableBoxMesh::clone() const
{
    return std::make_unique<MovableBoxMesh>( *this );
}


/****************************************************************
 * Return the coordinate                                         *
 ****************************************************************/
void MovableBoxMesh::coord( const MeshElementIndex &index, double *pos ) const
{
    AMP_ASSERT( index.type() == AMP::Mesh::GeomType::Vertex );
    size_t i = AMP::Utilities::findfirst( d_index, index );
    if ( d_index[i] != index ) {
        char msg[1024];
        snprintf( msg,
                  sizeof msg,
                  "Did not find element %s in mesh MovableBoxMesh(%i,%i,%i,%i%i%i)",
                  index.print().data(),
                  d_globalSize[0],
                  d_globalSize[1],
                  d_globalSize[2],
                  d_isPeriodic[0] ? 1 : 0,
                  d_isPeriodic[1] ? 1 : 0,
                  d_isPeriodic[2] ? 1 : 0 );
        AMP_ERROR( msg );
    }
    for ( int d = 0; d < PhysicalDim; d++ )
        pos[d] = d_coord( d, i );
}


/****************************************************************
 * Return the logical coordinates                                *
 ****************************************************************/
AMP::Geometry::Point MovableBoxMesh::physicalToLogical( const AMP::Geometry::Point & ) const
{
    AMP_ERROR( "physicalToLogical is not supported in MovableBoxMesh" );
    return AMP::Geometry::Point();
}


/****************************************************************
 * Check if two meshes are equal                                 *
 ****************************************************************/
bool MovableBoxMesh::operator==( const Mesh &rhs ) const
{
    // Check base class variables
    if ( !BoxMesh::operator==( rhs ) )
        return false;
    // Check if we can cast to a MovableBoxMesh
    auto mesh = dynamic_cast<const MovableBoxMesh *>( &rhs );
    if ( !mesh )
        return false;
    // Perform final comparisons
    bool test = d_index == mesh->d_index;
    test &= d_coord == mesh->d_coord;
    test &= d_ids == mesh->d_ids;
    return test;
}


} // namespace AMP::Mesh
