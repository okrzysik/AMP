#ifndef included_AMP_BoxMeshInline
#define included_AMP_BoxMeshInline

#include "ampmesh/structured/BoxMesh.h"

#include <algorithm>


namespace AMP {
namespace Mesh {


// BoxMesh::Box
BoxMesh::Box::Box()
{
    first[0] = 0;
    first[1] = 0;
    first[2] = 0;
    last[0]  = 0;
    last[1]  = 0;
    last[2]  = 0;
}
BoxMesh::Box::Box( int ifirst, int ilast, int jfirst, int jlast, int kfirst, int klast )
{
    first[0] = ifirst;
    first[1] = jfirst;
    first[2] = kfirst;
    last[0]  = ilast;
    last[1]  = jlast;
    last[2]  = klast;
}


// BoxMesh::MeshElementIndex
BoxMesh::MeshElementIndex::MeshElementIndex() : d_type( 0 ), d_side( 0 )
{
    d_index[0] = 0;
    d_index[1] = 0;
    d_index[2] = 0;
}
BoxMesh::MeshElementIndex::MeshElementIndex(
    GeomType type_in, unsigned char side_in, int x, int y, int z )
    : d_type( static_cast<unsigned char>( type_in ) ), d_side( side_in )
{
    d_index[0] = x;
    d_index[1] = y;
    d_index[2] = z;
}
void BoxMesh::MeshElementIndex::reset(
    GeomType type_in, unsigned char side_in, int x, int y, int z )
{
    d_type = static_cast<unsigned char>( type_in );
    d_side = side_in;
    d_index[0] = x;
    d_index[1] = y;
    d_index[2] = z;
}
inline bool BoxMesh::MeshElementIndex::operator==( const MeshElementIndex &rhs ) const
{
    return d_type == rhs.d_type && d_side == rhs.d_side && d_index[0] == rhs.d_index[0] &&
           d_index[1] == rhs.d_index[1] && d_index[2] == rhs.d_index[2];
}
inline bool BoxMesh::MeshElementIndex::operator!=( const MeshElementIndex &rhs ) const
{
    return d_type != rhs.d_type || d_side != rhs.d_side || d_index[0] != rhs.d_index[0] ||
           d_index[1] != rhs.d_index[1] || d_index[2] != rhs.d_index[2];
}
inline bool BoxMesh::MeshElementIndex::operator>( const MeshElementIndex &rhs ) const
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
inline bool BoxMesh::MeshElementIndex::operator>=( const MeshElementIndex &rhs ) const
{
    return this->operator>( rhs ) || this->operator==( rhs );
}
inline bool BoxMesh::MeshElementIndex::operator<( const MeshElementIndex &rhs ) const
{
    return !this->operator>( rhs ) && !this->operator==( rhs );
}
inline bool BoxMesh::MeshElementIndex::operator<=( const MeshElementIndex &rhs ) const
{
    return !this->operator>( rhs );
}
inline size_t BoxMesh::MeshElementIndex::numElements( const MeshElementIndex& first, const MeshElementIndex& last )
{
    if ( last.d_index[0]<first.d_index[0] || last.d_index[1]<first.d_index[1] || last.d_index[2]<first.d_index[2] )
        return 0;
    return (last.d_index[0]-first.d_index[0]+1) *
           (last.d_index[1]-first.d_index[1]+1) *
           (last.d_index[2]-first.d_index[2]+1);
}


// getGlobalBox, getLocalBox
inline std::vector<bool> BoxMesh::periodic() const
{
    std::vector<bool> per( static_cast<int>( GeomDim ) );
    for ( int d = 0; d < static_cast<int>( GeomDim ); d++ )
        per[d]  = d_isPeriodic[d];
    return per;
}
inline std::vector<size_t> BoxMesh::size() const
{
    std::vector<size_t> size( static_cast<int>( GeomDim ) );
    for ( int d = 0; d < static_cast<int>( GeomDim ); d++ )
        size[d]  = d_globalSize[d];
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
inline std::array<int,6> BoxMesh::getLocalBlock( unsigned int rank ) const
{
    int p[3];
    p[0] = rank % d_numBlocks[0];
    p[1] = rank/d_numBlocks[0] % d_numBlocks[1];
    p[2] = rank/(d_numBlocks[0]*d_numBlocks[1]);
    AMP_ASSERT( p[2] < d_numBlocks[2] );
    std::array<int,6> range;
    range.fill(0);
    for ( int d = 0; d < GeomDim; d++ ) {
        int size = (d_globalSize[d]+d_numBlocks[d]-1)/d_numBlocks[d];
        range[2*d+0] = p[d]*size;
        range[2*d+1] = std::min((p[d]+1)*size-1,d_globalSize[d]-1);
    }
    return range;
}


/****************************************************************
* Convert between the different id types                        *
****************************************************************/
inline MeshElementID BoxMesh::convert( const BoxMesh::MeshElementIndex& index ) const
{
    int size[3] = { (d_globalSize[0]+d_numBlocks[0]-1)/d_numBlocks[0],
                    (d_globalSize[1]+d_numBlocks[1]-1)/d_numBlocks[1],
                    (d_globalSize[2]+d_numBlocks[2]-1)/d_numBlocks[2] };
    int i = index.index(0);
    int j = index.index(1);
    int k = index.index(2);
    int px = std::min(i/size[0],d_numBlocks[0]-1);
    int py = std::min(j/size[1],d_numBlocks[1]-1);
    int pz = std::min(k/size[2],d_numBlocks[2]-1);
    i -= size[0]*px;
    j -= size[1]*py;
    k -= size[2]*pz;
    unsigned int local_id = i + (size[0]+1)*( j + (size[1]+1)*( k + (size[2]+1)*index.side() ) );
    int owner_rank = px + py*d_numBlocks[0] + pz*d_numBlocks[0]*d_numBlocks[1];
    bool is_local = owner_rank == d_comm.getRank();
    return MeshElementID( is_local, (GeomType) index.type(), local_id, owner_rank, d_meshID );
}
inline BoxMesh::MeshElementIndex BoxMesh::convert( const MeshElementID& id ) const
{
    int rank = id.owner_rank();
    int proc[3] = { rank % d_numBlocks[0],
                    rank/d_numBlocks[0] % d_numBlocks[1],
                    rank/(d_numBlocks[0]*d_numBlocks[1]) };
    int size[3] = { (d_globalSize[0]+d_numBlocks[0]-1)/d_numBlocks[0],
                    (d_globalSize[1]+d_numBlocks[1]-1)/d_numBlocks[1],
                    (d_globalSize[2]+d_numBlocks[2]-1)/d_numBlocks[2] };
    size_t ijk = id.local_id();
    int i = ijk % (size[0]+1);  ijk /= (size[0]+1);
    int j = ijk % (size[1]+1);  ijk /= (size[1]+1);
    int k = ijk % (size[2]+1);  ijk /= (size[2]+1);
    int side = ijk;
    i += proc[0]*size[0];
    j += proc[1]*size[1];
    k += proc[2]*size[2];
    return MeshElementIndex( id.type(), side, i, j, k );
}


} // Mesh namespace
} // AMP namespace

#endif
