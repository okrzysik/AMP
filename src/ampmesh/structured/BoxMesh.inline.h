#ifndef included_AMP_BoxMeshInline
#define included_AMP_BoxMeshInline

#include "AMP/ampmesh/structured/BoxMesh.h"

#include <algorithm>


namespace AMP {
namespace Mesh {


/****************************************************************
 * Box                                                           *
 ****************************************************************/
constexpr BoxMesh::Box::Box() : first{ 0, 0, 0 }, last{ 0, 0, 0 } {}
constexpr BoxMesh::Box::Box( int ifirst, int ilast, int jfirst, int jlast, int kfirst, int klast )
    : first{ ifirst, jfirst, kfirst }, last{ ilast, jlast, klast }
{
}
constexpr ArraySize BoxMesh::Box::size() const
{
    return { (size_t) last[0] - first[0] + 1,
             (size_t) last[1] - first[1] + 1,
             (size_t) last[2] - first[2] + 1 };
}


/****************************************************************
 * MeshElementIndex                                              *
 ****************************************************************/
constexpr BoxMesh::MeshElementIndex::MeshElementIndex()
    : d_type( 0 ), d_side( 255 ), d_index{ 0, 0, 0 }
{
}
constexpr BoxMesh::MeshElementIndex::MeshElementIndex(
    GeomType type_in, uint8_t side_in, int x, int y, int z )
    : d_type( static_cast<uint8_t>( type_in ) ), d_side( side_in ), d_index{ x, y, z }
{
}
constexpr void
BoxMesh::MeshElementIndex::reset( GeomType type_in, uint8_t side_in, int x, int y, int z )
{
    d_type     = static_cast<uint8_t>( type_in );
    d_side     = side_in;
    d_index[0] = x;
    d_index[1] = y;
    d_index[2] = z;
}
constexpr bool BoxMesh::MeshElementIndex::operator==( const MeshElementIndex &rhs ) const
{
    return d_type == rhs.d_type && d_side == rhs.d_side && d_index[0] == rhs.d_index[0] &&
           d_index[1] == rhs.d_index[1] && d_index[2] == rhs.d_index[2];
}
constexpr bool BoxMesh::MeshElementIndex::operator!=( const MeshElementIndex &rhs ) const
{
    return d_type != rhs.d_type || d_side != rhs.d_side || d_index[0] != rhs.d_index[0] ||
           d_index[1] != rhs.d_index[1] || d_index[2] != rhs.d_index[2];
}
constexpr bool BoxMesh::MeshElementIndex::operator>( const MeshElementIndex &rhs ) const
{
    if ( d_type < rhs.d_type ) {
        return false;
    } else if ( d_type > rhs.d_type ) {
        return true;
    }
    if ( d_side < rhs.d_side ) {
        return false;
    } else if ( d_side > rhs.d_side ) {
        return true;
    }
    for ( int i = 2; i >= 0; i-- ) {
        if ( d_index[i] < rhs.d_index[i] ) {
            return false;
        } else if ( d_index[i] > rhs.d_index[i] ) {
            return true;
        }
    }
    return false;
}
constexpr bool BoxMesh::MeshElementIndex::operator>=( const MeshElementIndex &rhs ) const
{
    return this->operator>( rhs ) || this->operator==( rhs );
}
constexpr bool BoxMesh::MeshElementIndex::operator<( const MeshElementIndex &rhs ) const
{
    return !this->operator>( rhs ) && !this->operator==( rhs );
}
constexpr bool BoxMesh::MeshElementIndex::operator<=( const MeshElementIndex &rhs ) const
{
    return !this->operator>( rhs );
}
constexpr size_t BoxMesh::MeshElementIndex::numElements( const MeshElementIndex &first,
                                                         const MeshElementIndex &last )
{
    if ( last.d_index[0] < first.d_index[0] || last.d_index[1] < first.d_index[1] ||
         last.d_index[2] < first.d_index[2] )
        return 0;
    return ( last.d_index[0] - first.d_index[0] + 1 ) * ( last.d_index[1] - first.d_index[1] + 1 ) *
           ( last.d_index[2] - first.d_index[2] + 1 );
}


/****************************************************************
 * getGlobalBox, getLocalBox                                     *
 ****************************************************************/
inline std::vector<bool> BoxMesh::periodic() const
{
    std::vector<bool> per( static_cast<int>( GeomDim ) );
    for ( int d = 0; d < static_cast<int>( GeomDim ); d++ )
        per[d] = d_isPeriodic[d];
    return per;
}
inline std::vector<size_t> BoxMesh::size() const
{
    std::vector<size_t> size( static_cast<int>( GeomDim ) );
    for ( int d = 0; d < static_cast<int>( GeomDim ); d++ )
        size[d] = d_globalSize[d];
    return size;
}
inline BoxMesh::Box BoxMesh::getGlobalBox( int gcw ) const
{
    Box box;
    for ( int d = 0; d < static_cast<int>( GeomDim ); d++ ) {
        box.first[d] = -gcw;
        box.last[d]  = d_globalSize[d] + gcw - 1;
        if ( !d_isPeriodic[d] ) {
            box.first[d] = std::max( box.first[d], 0 );
            box.last[d]  = std::min( box.last[d], d_globalSize[d] - 1 );
        }
    }
    return box;
}

inline BoxMesh::Box BoxMesh::getLocalBox( int gcw ) const
{
    auto range = getLocalBlock( d_comm.getRank() );
    Box box;
    for ( int d = 0; d < static_cast<int>( GeomDim ); d++ ) {
        box.first[d] = range[2 * d + 0] - gcw;
        box.last[d]  = range[2 * d + 1] + gcw;
        if ( !d_isPeriodic[d] ) {
            box.first[d] = std::max( box.first[d], 0 );
            box.last[d]  = std::min( box.last[d], d_globalSize[d] - 1 );
        }
    }
    return box;
}


/****************************************************************
 * Helper function to return the indices of the local block      *
 * owned by the given processor                                  *
 ****************************************************************/
inline std::array<int, 6> BoxMesh::getLocalBlock( unsigned int rank ) const
{
    int p[3];
    p[0] = rank % d_numBlocks[0];
    p[1] = rank / d_numBlocks[0] % d_numBlocks[1];
    p[2] = rank / ( d_numBlocks[0] * d_numBlocks[1] );
    AMP_ASSERT( p[2] < d_numBlocks[2] );
    std::array<int, 6> range;
    range.fill( 0 );
    for ( int d = 0; d < static_cast<int>( GeomDim ); d++ ) {
        int size         = ( d_globalSize[d] + d_numBlocks[d] - 1 ) / d_numBlocks[d];
        range[2 * d + 0] = p[d] * size;
        range[2 * d + 1] = std::min( ( p[d] + 1 ) * size - 1, d_globalSize[d] - 1 );
    }
    return range;
}


/****************************************************************
 * Convert between the different id types                        *
 ****************************************************************/
inline MeshElementID BoxMesh::convert( const BoxMesh::MeshElementIndex &index ) const
{
    int i  = index.index( 0 );
    int j  = index.index( 1 );
    int k  = index.index( 2 );
    int px = std::min( i / d_blockSize[0], d_numBlocks[0] - 1 );
    int py = std::min( j / d_blockSize[1], d_numBlocks[1] - 1 );
    int pz = std::min( k / d_blockSize[2], d_numBlocks[2] - 1 );
    i -= d_blockSize[0] * px;
    j -= d_blockSize[1] * py;
    k -= d_blockSize[2] * pz;
    unsigned int local_id =
        i + ( d_blockSize[0] + 1 ) *
                ( j + ( d_blockSize[1] + 1 ) * ( k + ( d_blockSize[2] + 1 ) * index.side() ) );
    int owner_rank = px + d_numBlocks[0] * ( py + d_numBlocks[1] * pz );
    bool is_local  = owner_rank == d_comm.getRank();
    return MeshElementID( is_local, (GeomType) index.type(), local_id, owner_rank, d_meshID );
}
inline BoxMesh::MeshElementIndex BoxMesh::convert( const MeshElementID &id ) const
{
    int rank    = id.owner_rank();
    int proc[3] = { rank % d_numBlocks[0],
                    rank / d_numBlocks[0] % d_numBlocks[1],
                    rank / ( d_numBlocks[0] * d_numBlocks[1] ) };
    size_t ijk  = id.local_id();
    int i       = ijk % ( d_blockSize[0] + 1 );
    ijk /= ( d_blockSize[0] + 1 );
    int j = ijk % ( d_blockSize[1] + 1 );
    ijk /= ( d_blockSize[1] + 1 );
    int k = ijk % ( d_blockSize[2] + 1 );
    ijk /= ( d_blockSize[2] + 1 );
    int side = ijk;
    i += proc[0] * d_blockSize[0];
    j += proc[1] * d_blockSize[1];
    k += proc[2] * d_blockSize[2];
    return MeshElementIndex( id.type(), side, i, j, k );
}


} // namespace Mesh
} // namespace AMP

#endif
