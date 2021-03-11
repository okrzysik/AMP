#include "AMP/ampmesh/structured/MovableBoxMesh.h"

#include "AMP/ampmesh/MultiIterator.h"
#include "AMP/ampmesh/structured/BoxMesh.h"
#include "AMP/ampmesh/structured/structuredMeshElement.h"
#include "AMP/ampmesh/structured/structuredMeshIterator.h"

#ifdef USE_AMP_VECTORS
#include "AMP/vectors/Variable.h"
#include "AMP/vectors/Vector.h"
#include "AMP/vectors/VectorBuilder.h"
#endif
#ifdef USE_AMP_DISCRETIZATION
#include "AMP/discretization/DOF_Manager.h"
#include "AMP/discretization/simpleDOF_Manager.h"
#endif

namespace AMP {
namespace Mesh {


/****************************************************************
 * Constructors                                                  *
 ****************************************************************/
MovableBoxMesh::MovableBoxMesh( const AMP::Mesh::BoxMesh &mesh ) : BoxMesh( mesh ), d_pos_hash( 0 )
{
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
    d_coord.resize( d_index.size(), { 0, 0, 0 } );
    for ( size_t i = 0; i < d_index.size(); i++ )
        mesh.coord( d_index[i], d_coord[i].data() );
}


/****************************************************************
 * Functions to displace the mesh                                *
 ****************************************************************/
Mesh::Movable MovableBoxMesh::isMeshMovable() const { return Mesh::Movable::Deform; }
uint64_t MovableBoxMesh::positionHash() const { return d_pos_hash; }
void MovableBoxMesh::displaceMesh( const std::vector<double> &x )
{
    AMP_ASSERT( x.size() == PhysicalDim );
    for ( size_t i = 0; i < d_coord.size(); i++ ) {
        for ( int d = 0; d < PhysicalDim; d++ )
            d_coord[i][d] += x[d];
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
#ifdef USE_AMP_VECTORS
void MovableBoxMesh::displaceMesh( const AMP::LinearAlgebra::Vector::const_shared_ptr x )
{
#ifdef USE_AMP_DISCRETIZATION
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
    for ( size_t i = 0; i < d_coord.size(); i++ ) {
        MeshElementID id = structuredMeshElement( d_index[i], this ).globalID();
        DOFs->getDOFs( id, dofs );
        AMP_ASSERT( dofs.size() == PhysicalDim );
        displacement->getValuesByGlobalID( (int) PhysicalDim, &dofs[0], &disp[0] );
        for ( int d = 0; d < PhysicalDim; d++ )
            d_coord[i][d] += disp[d];
    }
    // Compute the new bounding box of the mesh
    d_box_local = std::vector<double>( 2 * PhysicalDim );
    for ( int d = 0; d < PhysicalDim; d++ ) {
        d_box_local[2 * d + 0] = 1e100;
        d_box_local[2 * d + 1] = -1e100;
    }
    for ( size_t i = 0; i < d_coord.size(); i++ ) {
        for ( int d = 0; d < PhysicalDim; d++ ) {
            d_box_local[2 * d + 0] = std::min( d_box_local[2 * d + 0], d_coord[i][d] );
            d_box_local[2 * d + 1] = std::max( d_box_local[2 * d + 1], d_coord[i][d] );
        }
    }
    d_box = std::vector<double>( PhysicalDim * 2 );
    for ( int i = 0; i < PhysicalDim; i++ ) {
        d_box[2 * i + 0] = d_comm.minReduce( d_box_local[2 * i + 0] );
        d_box[2 * i + 1] = d_comm.maxReduce( d_box_local[2 * i + 1] );
    }
    d_pos_hash++;
#else
    AMP_ERROR( "displaceMesh requires DISCRETIZATION" );
#endif
}
#endif


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
    AMP_ASSERT( d_index[i] == index );
    for ( int d = 0; d < PhysicalDim; d++ )
        pos[d] = d_coord[i][d];
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
    return true;
}


} // namespace Mesh
} // namespace AMP
