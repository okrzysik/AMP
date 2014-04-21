#ifndef included_AMP_BoxMeshInline
#define included_AMP_BoxMeshInline

#include "ampmesh/structured/BoxMesh.h"
#include "ampmesh/structured/structuredMeshElement.h"


namespace AMP {
namespace Mesh {


// BoxMesh::Box
BoxMesh::Box::Box( )
{
    first[0] = 0;
    first[1] = 0;
    first[2] = 0;
    last[0] = 0;
    last[1] = 0;
    last[2] = 0;
}
BoxMesh::Box::Box( int ifirst, int ilast, int jfirst, int jlast, int kfirst, int klast )
{
    first[0] = ifirst;
    first[1] = jfirst;
    first[2] = kfirst;
    last[0] = ilast;
    last[1] = jlast;
    last[2] = klast;
}


// BoxMesh::MeshElementIndex
BoxMesh::MeshElementIndex::MeshElementIndex():
    type(0),
    side(0)
{
    index[0] = 0;
    index[1] = 0;
    index[2] = 0;
}
BoxMesh::MeshElementIndex::MeshElementIndex(GeomType type_in, unsigned char side_in, int x, int y, int z):
    type(static_cast<unsigned char>(type_in)),
    side(side_in)
{
    index[0] = x;
    index[1] = y;
    index[2] = z;
}
inline bool BoxMesh::MeshElementIndex::operator== (const MeshElementIndex& rhs ) const 
{
    return type==rhs.type && side==rhs.side && index[0]==rhs.index[0] 
        && index[1]==rhs.index[1] && index[2]==rhs.index[2];
}
inline bool BoxMesh::MeshElementIndex::operator!= (const MeshElementIndex& rhs ) const 
{
    return type!=rhs.type || side!=rhs.side || index[0]!=rhs.index[0] 
        || index[1]!=rhs.index[1] || index[2]!=rhs.index[2];
}
inline bool BoxMesh::MeshElementIndex::operator> (const MeshElementIndex& rhs ) const 
{
    if ( type < rhs.type ) { return false; }
    else if ( type > rhs.type ) { return true; }
    if ( side < rhs.side ) { return false; }
    else if ( side > rhs.side ) { return true; }
    for (int i=2; i>=0; i--) {
        if ( index[i] < rhs.index[i] ) { return false; }
        else if ( index[i] > rhs.index[i] ) { return true; }
    }
    return false;
}
inline bool BoxMesh::MeshElementIndex::operator>= (const MeshElementIndex& rhs ) const 
{
    return this->operator>(rhs) || this->operator==(rhs);
}
inline bool BoxMesh::MeshElementIndex::operator< (const MeshElementIndex& rhs ) const 
{
    return !this->operator>(rhs) && !this->operator==(rhs);
}
inline bool BoxMesh::MeshElementIndex::operator<= (const MeshElementIndex& rhs ) const 
{
    return !this->operator>(rhs);
}


// getGlobalBox, getLocalBox
inline std::vector<bool> BoxMesh::periodic() const
{
    std::vector<bool> per(static_cast<int>(GeomDim));
    for (int d=0; d<static_cast<int>(GeomDim); d++)
        per[d] = d_isPeriodic[2*d+1];
    return per;
}
inline BoxMesh::Box BoxMesh::getGlobalBox( int gcw ) const
{
    Box box;
    for (int d=0; d<static_cast<int>(GeomDim); d++) {
        box.first[d] = -gcw;
        box.last[d] = d_size[d]+gcw-1;
    }
    return box;
}

inline BoxMesh::Box BoxMesh::getLocalBox( int gcw ) const
{
    std::vector<int> range = getLocalBlock(d_comm.getRank());
    Box box;
    for (int d=0; d<static_cast<int>(GeomDim); d++) {
        box.first[d] = range[2*d+0]-gcw;
        box.last[d]  = range[2*d+1]+gcw-1;
        if ( box.first[d]+d_size[d] <= box.last[d] ) {
            box.first[d] = 0;
            box.last[d] = d_size[d]-1;
        }
    }
    return box;
}


} // Mesh namespace
} // AMP namespace

#endif

