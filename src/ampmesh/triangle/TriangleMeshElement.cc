#include "AMP/ampmesh/triangle/TriangleMeshElement.h"
#include "AMP/ampmesh/triangle/TriangleMeshIterator.h"
#include "AMP/ampmesh/triangle/TriangleMesh.h"
#include "AMP/utils/Utilities.h"


namespace AMP {
namespace Mesh {


/********************************************************
 * Create a unique id for each class                     *
 ********************************************************/
template<>
constexpr uint32_t TriangleMeshElement<1,1>::getTypeID() { return AMP::Utilities::hash_char( "TriangleMeshElement<1,1>" ); }
template<>
constexpr uint32_t TriangleMeshElement<1,2>::getTypeID() { return AMP::Utilities::hash_char( "TriangleMeshElement<1,2>" ); }
template<>
constexpr uint32_t TriangleMeshElement<1,3>::getTypeID() { return AMP::Utilities::hash_char( "TriangleMeshElement<1,3>" ); }
template<>
constexpr uint32_t TriangleMeshElement<2,2>::getTypeID() { return AMP::Utilities::hash_char( "TriangleMeshElement<2,2>" ); }
template<>
constexpr uint32_t TriangleMeshElement<2,3>::getTypeID() { return AMP::Utilities::hash_char( "TriangleMeshElement<2,3>" ); }
template<>
constexpr uint32_t TriangleMeshElement<3,3>::getTypeID() { return AMP::Utilities::hash_char( "TriangleMeshElement<3,3>" ); }


/********************************************************
 * Constructors                                          *
 ********************************************************/
template<size_t NG, size_t NP>
TriangleMeshElement<NG,NP>::TriangleMeshElement()
{
    typeID     = getTypeID();
    element    = nullptr;
    d_mesh     = nullptr;
}
template<size_t NG, size_t NP>
TriangleMeshElement<NG,NP>::TriangleMeshElement( const MeshElementID& id, const TriangleMesh<NG,NP> *mesh )
{
    typeID     = getTypeID();
    element    = nullptr;
    d_globalID = id;
    d_mesh     = mesh;
}
template<size_t NG, size_t NP>
TriangleMeshElement<NG,NP>::TriangleMeshElement( const TriangleMeshElement& rhs ): MeshElement()
{
    typeID     = rhs.typeID;
    element    = nullptr;
    d_globalID = rhs.d_globalID;
    d_mesh     = rhs.d_mesh;
}
template<size_t NG, size_t NP>
TriangleMeshElement<NG,NP>::TriangleMeshElement( TriangleMeshElement&& rhs ): MeshElement()
{
    typeID     = rhs.typeID;
    element    = nullptr;
    d_globalID = rhs.d_globalID;
    d_mesh     = rhs.d_mesh;
}
template<size_t NG, size_t NP>
TriangleMeshElement<NG,NP>& TriangleMeshElement<NG,NP>::operator=( const TriangleMeshElement& rhs )
{
    if ( &rhs == this ) 
        return *this;
    typeID     = rhs.typeID;
    element    = nullptr;
    d_globalID = rhs.d_globalID;
    d_mesh     = rhs.d_mesh;
    return *this;
}
template<size_t NG, size_t NP>
TriangleMeshElement<NG,NP>& TriangleMeshElement<NG,NP>::operator=( TriangleMeshElement&& rhs )
{
    if ( &rhs == this ) 
        return *this;
    typeID     = rhs.typeID;
    element    = nullptr;
    d_globalID = rhs.d_globalID;
    d_mesh     = rhs.d_mesh;
    return *this;
}


/****************************************************************
 * Function to clone the element                                 *
 ****************************************************************/
template<size_t NG, size_t NP>
MeshElement *TriangleMeshElement<NG,NP>::clone() const { return new TriangleMeshElement<NG,NP>( *this ); }


/****************************************************************
 * Return the global rank of the owner rank                      *
 ****************************************************************/
template<size_t NG, size_t NP>
unsigned int TriangleMeshElement<NG,NP>::globalOwnerRank() const
{
    return d_mesh->getComm().globalRanks()[d_globalID.owner_rank()];
}


/****************************************************************
 * Function to get the elements composing the current element    *
 ****************************************************************/
template<size_t NG, size_t NP>
void TriangleMeshElement<NG,NP>::getElementsID( const GeomType type, std::vector<MeshElementID> &ID ) const
{
    // Number of elements composing a given type [current][desired]
    constexpr uint8_t N_elements[4][4] = { { 1, 0, 0, 0 }, { 2, 1, 0, 0 }, { 3, 3, 2, 0 }, { 3, 4, 6, 4 } };
    int N = N_elements[static_cast<int>(d_globalID.type())][static_cast<int>(type)];
    // Get the element ids
    ID.resize( N );
    d_mesh->getElementsIDs( d_globalID, type, ID.data() );
}
template<size_t NG, size_t NP>
void TriangleMeshElement<NG,NP>::getElements( const GeomType type, std::vector<MeshElement> &children ) const
{
    // Number of elements composing a given type [current][desired]
    constexpr uint8_t N_elements[4][4] = { { 1, 0, 0, 0 }, { 2, 1, 0, 0 }, { 3, 3, 2, 0 }, { 3, 4, 6, 4 } };
    int N = N_elements[static_cast<int>(d_globalID.type())][static_cast<int>(type)];
    // Get the element ids
    MeshElementID tmp[6];
    d_mesh->getElementsIDs( d_globalID, type, tmp );
    // Create the mesh elements
    children.resize( N );
    for ( int i=0; i<N; i++)
        children[i] = TriangleMeshElement<NG,NP>( tmp[i], d_mesh );
}


/****************************************************************
 * Function to get the neighboring elements                      *
 ****************************************************************/
template<size_t NG, size_t NP>
void TriangleMeshElement<NG,NP>::getNeighbors( std::vector<MeshElement::shared_ptr> &neighbors ) const
{
    NULL_USE( neighbors );
    AMP_ERROR( "Not finished" );
}


/****************************************************************
 * Functions to get basic element properties                     *
 ****************************************************************/
static inline double dot( const std::array<double,1>& x, const std::array<double,1>& y )
{
    return x[0] * y[0];
}
static inline double dot( const std::array<double,2>& x, const std::array<double,2>& y )
{
    return x[0] * y[0] + x[1] * y[1];
}
static inline double dot( const std::array<double,3>& x, const std::array<double,3>& y )
{
    return x[0] * y[0] + x[1] * y[1] + x[2] * y[2];
}
static inline std::array<double,1> operator-( const std::array<double,1>& x, const std::array<double,1>& y )
{
    return { x[0] - y[0] };
}
static inline std::array<double,2> operator-( const std::array<double,2>& x, const std::array<double,2>& y )
{
    return { x[0] - y[0], x[1] - y[1] };
}
static inline std::array<double,3> operator-( const std::array<double,3>& x, const std::array<double,3>& y )
{
    return { x[0] - y[0], x[1] - y[1], x[2] - y[2] };
}
template<size_t NG, size_t NP>
double TriangleMeshElement<NG,NP>::volume() const
{
    if ( NG == 0 )
        return 0;
    MeshElementID ids[NG+1];
    d_mesh->getElementsIDs( d_globalID, GeomType::Vertex, ids );
    std::array<double,NP> coord[NG+1];
    for (size_t i=0; i<=NG; i++)
        coord[i] = d_mesh->getPos( ids[i] );
    if ( NG == 1 ) {
        double dist = 0;
        for (size_t d=0; d<NP; d++)
            dist += ( coord[0][d] - coord[1][d] ) * ( coord[0][d] - coord[1][d] );
        return dist;
    }
    if ( NG == 2 ) {
        auto AB = coord[0] - coord[1];
        auto AC = coord[0] - coord[2];
        double t = dot( AB, AC );
        return 0.5 * sqrt( dot(AB,AB) * dot(AC,AC) - t*t );
    }
    AMP_ERROR("Not finished");
    return 0;
}
template<size_t NG, size_t NP>
void TriangleMeshElement<NG,NP>::coord( size_t &N, double *x ) const
{
    if ( d_globalID.type() == GeomType::Vertex ) {
        N = NP;
        auto& x2 = d_mesh->getPos( d_globalID );
        for ( size_t d = 0; d<NP; d++)
            x[d] = x2[d];
    } else {
        AMP_ERROR( "coord is only valid for verticies" );
    }
}
template<size_t NG, size_t NP>
void TriangleMeshElement<NG,NP>::centroid( size_t &N, double *x ) const
{
    N = NP;
    if ( d_globalID.type() == GeomType::Vertex ) {
        auto& x2 = d_mesh->getPos( d_globalID );
        for ( size_t d = 0; d<NP; d++)
            x[d] = x2[d];
        return;
    }
    MeshElementID ids[NG+1];
    d_mesh->getElementsIDs( d_globalID, GeomType::Vertex, ids );
    for ( size_t d = 0; d<NP; d++)
        x[d] = 0;
    for ( size_t i=0; i<=NG; i++) {
        auto& x2 = d_mesh->getPos( ids[i] );
        for ( size_t d=0; d<NP; d++)
            x[d] += x2[d];
    }
    for ( size_t d=0; d<NP; d++)
        x[d] /= ( NG + 1 );
}
template<size_t NG, size_t NP>
bool TriangleMeshElement<NG,NP>::containsPoint( const std::vector<double> &pos, double TOL ) const
{
    NULL_USE( pos );
    NULL_USE( TOL );
    AMP_ERROR( "Not finished" );
    return false;
}
template<size_t NG, size_t NP>
bool TriangleMeshElement<NG,NP>::isOnSurface() const
{
    AMP_ERROR( "Not finished" );
    return false;
}
template<size_t NG, size_t NP>
bool TriangleMeshElement<NG,NP>::isOnBoundary( int id ) const
{
    NULL_USE( id );
    AMP_ERROR( "Not finished" );
    return false;
}
template<size_t NG, size_t NP>
bool TriangleMeshElement<NG,NP>::isInBlock( int id ) const
{
    NULL_USE( id );
    AMP_ERROR( "Not finished" );
    return false;
}


/********************************************************
 *  Explicit instantiations of TriangleMeshElement       *
 ********************************************************/
template class TriangleMeshElement<1,1>;
template class TriangleMeshElement<1,2>;
template class TriangleMeshElement<1,3>;
template class TriangleMeshElement<2,2>;
template class TriangleMeshElement<2,3>;
template class TriangleMeshElement<3,3>;


} // namespace Mesh
} // namespace AMP
