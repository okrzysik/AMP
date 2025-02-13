#include "AMP/discretization/boxMeshDOFManager.h"
#include "AMP/mesh/structured/BoxMesh.h"
#include "AMP/utils/AMP_MPI.h"
#include "AMP/utils/Utilities.h"
#include <set>
#include <vector>


namespace AMP::Discretization {


/****************************************************************
 * Find a value, returning the index if found or -1              *
 ****************************************************************/
template<class T>
static int64_t find( int64_t N, const T *x, const T &value )
{
    if ( N == 0 )
        return -1;
    if ( x[0] == value )
        return 0;
    int64_t lower = 0;
    int64_t upper = N - 1;
    while ( ( upper - lower ) > 1 ) {
        auto index = ( upper + lower ) / 2;
        if ( x[index] >= value )
            upper = index;
        else
            lower = index;
    }
    if ( x[upper] == value )
        return upper;
    return -1;
}
template<class T>
static inline int64_t find( const std::vector<T> &x, const T &value )
{
    return find( x.size(), x.data(), value );
}


/****************************************************************
 * Constructors                                                  *
 ****************************************************************/
boxMeshDOFManager::boxMeshDOFManager( std::shared_ptr<const AMP::Mesh::Mesh> mesh,
                                      AMP::Mesh::GeomType type,
                                      int gcw,
                                      int DOFsPerObject )
    : simpleDOFManager(
          mesh, mesh->getIterator( type, 0 ), mesh->getIterator( type, gcw ), type, DOFsPerObject ),
      d_gcw( gcw ),
      d_boxMesh( std::dynamic_pointer_cast<const AMP::Mesh::BoxMesh>( mesh ) )
{
    AMP_INSIST( d_boxMesh, "Mesh must be a BoxMesh" );

    // Initialize local data
    d_ifirst.resize( d_comm.getSize() );
    d_boxSize.resize( d_ifirst.size(), { { 0 } } );
    for ( size_t rank = 0; rank < d_ifirst.size(); rank++ ) {
        auto box   = d_boxMesh->getLocalBlock( rank );
        auto range = d_boxMesh->getIteratorRange( box, type, 0 );
        for ( size_t i = 0; i < range.size(); i++ ) {
            d_ifirst[rank][i]     = range[i].first;
            d_boxSize[rank][i][0] = range[i].second.index( 0 ) - range[i].first.index( 0 ) + 1;
            d_boxSize[rank][i][1] = range[i].second.index( 1 ) - range[i].first.index( 1 ) + 1;
            d_boxSize[rank][i][2] = range[i].second.index( 2 ) - range[i].first.index( 2 ) + 1;
        }
    }
    size_t start = beginDOF() / DOFsPerObject;
    d_start      = d_comm.allGather( start );

    // Perform some verification tests
#if ( defined( DEBUG ) || defined( _DEBUG ) ) && !defined( NDEBUG )
    for ( size_t i = 0; i < d_local_id.size(); i++ ) {
        auto dof = convert( d_local_id[i] );
        AMP_ASSERT( dof == i + start );
        // auto id  = convert( dof );
        // AMP_ASSERT( dof == i + d_begin && id == d_local_id[i] );
    }
    for ( size_t i = 0; i < d_remote_id.size(); i++ ) {
        auto dof = convert( d_remote_id[i] );
        if ( dof != d_remote_dof[i] )
            dof = convert( d_remote_id[i] );
        AMP_ASSERT( dof == d_remote_dof[i] );
        // auto id  = convert( dof );
        // AMP_ASSERT( dof == d_remote_dof[i] && id == d_local_id[i] );
    }
#endif
}


/****************************************************************
 * Get the array size for the local variables                    *
 ****************************************************************/
ArraySize boxMeshDOFManager::getArraySize() const
{
    if ( d_DOFsPerElement == 0 || d_local_id.empty() )
        return {};
    auto box   = d_boxMesh->getLocalBlock( d_comm.getRank() );
    auto range = d_boxMesh->getIteratorRange( box, d_type, 0 );
    if ( range.size() != 1u )
        return {};
    auto first = range[0].first;
    auto last  = range[0].second;
    int ndim   = static_cast<int>( d_boxMesh->getGeomType() );
    size_t Nx  = last.index( 0 ) - first.index( 0 ) + 1;
    size_t Ny  = last.index( 1 ) - first.index( 1 ) + 1;
    size_t Nz  = last.index( 2 ) - first.index( 2 ) + 1;
    ArraySize size;
    if ( d_DOFsPerElement == 1 )
        size = ArraySize( { Nx, Ny, Nz }, ndim );
    else
        size = ArraySize( { (size_t) d_DOFsPerElement, Nx, Ny, Nz }, ndim + 1 );
    AMP_ASSERT( size.length() == d_DOFsPerElement * d_local_id.size() );
    return size;
}


/****************************************************************
 * Convert indicies                                              *
 ****************************************************************/
size_t boxMeshDOFManager::convert( const AMP::Mesh::MeshElementID &id ) const
{
    auto index = d_boxMesh->convert( id );
    int rank   = id.owner_rank();
    int side   = index.side();
    size_t dof = 0;
    for ( int s = 0; s < side; s++ ) {
        auto N = d_boxSize[rank][s];
        dof += N[0] * N[1] * N[2];
    }
    auto ijk   = index.index();
    auto N     = d_boxSize[rank][side];
    auto first = d_ifirst[rank][side].index();
    dof +=
        ( ijk[0] - first[0] ) + ( ijk[1] - first[1] ) * N[0] + ( ijk[2] - first[2] ) * N[0] * N[1];
    return dof + d_start[rank];
}
AMP::Mesh::MeshElementID boxMeshDOFManager::convert( size_t ) const
{
    AMP_ERROR( "Not finished" );
    /*using MeshElementIndex = AMP::Mesh::BoxMesh::MeshElementIndex;
    int rank = id.owner_rank();
    for ( size_t side = 0; side < d_globalLast.size(); side++ ) {
        size_t Nt = MeshElementIndex::numElements( d_globalLast[side] );
        if ( dof >= Nt ) {
            dof -= Nt;
        } else {
            auto last = d_globalLast[side].index();
            int N[3] = { last[0] + 1, last[1] + 1, last[2] + 1 };
            int i = dof % N[0];
            dof /= N[0];
            int j = dof % N[1];
            dof /= N[1];
            int k = dof % N[2];
            dof /= N[2];
            AMP::Mesh::BoxMesh::MeshElementIndex index( d_type, side, i, j, k );
            return d_boxMesh->convert( index );
        }
    }*/
    return AMP::Mesh::MeshElementID();
}


/****************************************************************
 * Get the DOFs for the element                                  *
 * Note complete but promising, failure is how to eliminate      *
 *    elements beyond the gcw (remote)                           *
 ****************************************************************/
size_t boxMeshDOFManager::appendDOFs( const AMP::Mesh::MeshElementID &id,
                                      size_t *dofs,
                                      size_t N0,
                                      size_t capacity ) const
{
    // Check if the element should exist in the dof manager
    bool check = id.type() == d_type && id.meshID() == d_meshID && ( id.is_local() || d_gcw >= 1 );
    if ( !check )
        return 0;
    auto dof = d_DOFsPerElement * convert( id );
    for ( size_t j = 0, k = N0; j < d_DOFsPerElement && k < capacity; j++, k++, dof++ )
        dofs[k] = dof;
    return d_DOFsPerElement;
}


} // namespace AMP::Discretization
