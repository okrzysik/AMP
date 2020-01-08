#include "AMP/ampmesh/triangle/TriangleMeshElement.h"
#include "AMP/ampmesh/triangle/TriangleMesh.h"
#include "AMP/ampmesh/triangle/TriangleMeshIterator.h"
#include "AMP/utils/Utilities.h"


namespace AMP {
namespace Mesh {


/********************************************************
 * Create a unique id for each class                     *
 ********************************************************/
template<>
constexpr uint32_t TriangleMeshElement<1, 1>::getTypeID()
{
    return AMP::Utilities::hash_char( "TriangleMeshElement<1,1>" );
}
template<>
constexpr uint32_t TriangleMeshElement<1, 2>::getTypeID()
{
    return AMP::Utilities::hash_char( "TriangleMeshElement<1,2>" );
}
template<>
constexpr uint32_t TriangleMeshElement<1, 3>::getTypeID()
{
    return AMP::Utilities::hash_char( "TriangleMeshElement<1,3>" );
}
template<>
constexpr uint32_t TriangleMeshElement<2, 2>::getTypeID()
{
    return AMP::Utilities::hash_char( "TriangleMeshElement<2,2>" );
}
template<>
constexpr uint32_t TriangleMeshElement<2, 3>::getTypeID()
{
    return AMP::Utilities::hash_char( "TriangleMeshElement<2,3>" );
}
template<>
constexpr uint32_t TriangleMeshElement<3, 3>::getTypeID()
{
    return AMP::Utilities::hash_char( "TriangleMeshElement<3,3>" );
}


/********************************************************
 * Constructors                                          *
 ********************************************************/
template<size_t NG, size_t NP>
TriangleMeshElement<NG, NP>::TriangleMeshElement()
{
    typeID  = getTypeID();
    element = nullptr;
    d_mesh  = nullptr;
}
template<size_t NG, size_t NP>
TriangleMeshElement<NG, NP>::TriangleMeshElement( const MeshElementID &id,
                                                  const TriangleMesh<NG, NP> *mesh )
{
    typeID     = getTypeID();
    element    = nullptr;
    d_globalID = id;
    d_mesh     = mesh;
}
template<size_t NG, size_t NP>
TriangleMeshElement<NG, NP>::TriangleMeshElement( const TriangleMeshElement &rhs )
    : MeshElement(), d_mesh( rhs.d_mesh ), d_globalID( rhs.d_globalID )
{
    typeID  = getTypeID();
    element = rhs.element;
}
template<size_t NG, size_t NP>
TriangleMeshElement<NG, NP>::TriangleMeshElement( TriangleMeshElement &&rhs )
    : MeshElement(), d_mesh( rhs.d_mesh ), d_globalID{ rhs.d_globalID }
{
    typeID  = rhs.typeID;
    element = nullptr;
}
template<size_t NG, size_t NP>
TriangleMeshElement<NG, NP> &TriangleMeshElement<NG, NP>::
operator=( const TriangleMeshElement &rhs )
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
TriangleMeshElement<NG, NP> &TriangleMeshElement<NG, NP>::operator=( TriangleMeshElement &&rhs )
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
MeshElement *TriangleMeshElement<NG, NP>::clone() const
{
    return new TriangleMeshElement<NG, NP>( *this );
}


/****************************************************************
 * Return the global rank of the owner rank                      *
 ****************************************************************/
template<size_t NG, size_t NP>
unsigned int TriangleMeshElement<NG, NP>::globalOwnerRank() const
{
    return d_mesh->getComm().globalRanks()[d_globalID.owner_rank()];
}


/****************************************************************
 * Function to get the elements composing the current element    *
 ****************************************************************/
template<size_t NG, size_t NP>
void TriangleMeshElement<NG, NP>::getElementsID( const GeomType type,
                                                 std::vector<MeshElementID> &ID ) const
{
    // Number of elements composing a given type [current][desired]
    constexpr uint8_t N_elements[4][4] = {
        { 1, 0, 0, 0 }, { 2, 1, 0, 0 }, { 3, 3, 1, 0 }, { 4, 6, 4, 1 }
    };
    int N = N_elements[static_cast<int>( d_globalID.type() )][static_cast<int>( type )];
    // Get the element ids
    ElementID tmp[6];
    d_mesh->getElementsIDs( d_globalID.elemID(), type, tmp );
    ID.resize( N );
    for ( int i = 0; i < N; i++ )
        ID[i] = MeshElementID( d_globalID.meshID(), tmp[i] );
}
template<size_t NG, size_t NP>
void TriangleMeshElement<NG, NP>::getElements( const GeomType type,
                                               std::vector<MeshElement> &children ) const
{
    // Number of elements composing a given type [current][desired]
    constexpr uint8_t N_elements[4][4] = {
        { 1, 0, 0, 0 }, { 2, 1, 0, 0 }, { 3, 3, 1, 0 }, { 4, 6, 4, 1 }
    };
    int N = N_elements[static_cast<int>( d_globalID.type() )][static_cast<int>( type )];
    // Get the element ids
    ElementID tmp[6];
    d_mesh->getElementsIDs( d_globalID.elemID(), type, tmp );
    // Create the mesh elements
    auto meshID = d_globalID.meshID();
    children.resize( N );
    for ( int i = 0; i < N; i++ )
        children[i] = TriangleMeshElement<NG, NP>( MeshElementID( meshID, tmp[i] ), d_mesh );
}


/****************************************************************
 * Function to get the neighboring elements                      *
 ****************************************************************/
template<size_t NG, size_t NP>
void TriangleMeshElement<NG, NP>::getNeighbors(
    std::vector<MeshElement::shared_ptr> &neighbors ) const
{
    std::vector<ElementID> neighborIDs;
    d_mesh->getNeighborIDs( d_globalID.elemID(), neighborIDs );
    neighbors.resize( neighborIDs.size() );
    auto meshID = d_globalID.meshID();
    for ( size_t i = 0; i < neighborIDs.size(); i++ )
        neighbors[i].reset(
            new TriangleMeshElement<NG, NP>( MeshElementID( meshID, neighborIDs[i] ), d_mesh ) );
}


/****************************************************************
 * Functions to get basic element properties                     *
 ****************************************************************/
static inline double dot( const std::array<double, 1> &x, const std::array<double, 1> &y )
{
    return x[0] * y[0];
}
static inline double dot( const std::array<double, 2> &x, const std::array<double, 2> &y )
{
    return x[0] * y[0] + x[1] * y[1];
}
static inline double dot( const std::array<double, 3> &x, const std::array<double, 3> &y )
{
    return x[0] * y[0] + x[1] * y[1] + x[2] * y[2];
}
static inline std::array<double, 1> operator-( const std::array<double, 1> &x,
                                               const std::array<double, 1> &y )
{
    return { x[0] - y[0] };
}
static inline std::array<double, 2> operator-( const std::array<double, 2> &x,
                                               const std::array<double, 2> &y )
{
    return { x[0] - y[0], x[1] - y[1] };
}
static inline std::array<double, 3> operator-( const std::array<double, 3> &x,
                                               const std::array<double, 3> &y )
{
    return { x[0] - y[0], x[1] - y[1], x[2] - y[2] };
}
template<size_t NG, size_t NP>
double TriangleMeshElement<NG, NP>::volume() const
{
    if ( NG == 0 ) { // Replace if statements with if constexpr when supported
        return 0;
    } else if ( NG == 1 ) {
        ElementID ids[2];
        d_mesh->getElementsIDs( d_globalID.elemID(), GeomType::Vertex, ids );
        auto p1     = d_mesh->getPos( ids[0] );
        auto p2     = d_mesh->getPos( ids[1] );
        double dist = 0;
        for ( size_t d = 0; d < NP; d++ )
            dist += ( p1[d] - p2[d] ) * ( p1[d] - p2[d] );
        return dist;
    } else if ( NG == 2 ) {
        ElementID ids[3];
        d_mesh->getElementsIDs( d_globalID.elemID(), GeomType::Vertex, ids );
        auto p0  = d_mesh->getPos( ids[0] );
        auto AB  = p0 - d_mesh->getPos( ids[1] );
        auto AC  = p0 - d_mesh->getPos( ids[2] );
        double t = dot( AB, AC );
        return 0.5 * sqrt( dot( AB, AB ) * dot( AC, AC ) - t * t );
    } else {
        AMP_ERROR( "Not finished" );
        return 0;
    }
}
template<size_t NG, size_t NP>
Point TriangleMeshElement<NG, NP>::coord() const
{
    if ( d_globalID.type() == GeomType::Vertex ) {
        auto &x = d_mesh->getPos( d_globalID.elemID() );
        return Point( NP, x.data() );
    } else {
        AMP_ERROR( "coord is only valid for verticies" );
    }
    return Point();
}
template<size_t NG, size_t NP>
Point TriangleMeshElement<NG, NP>::centroid() const
{
    if ( d_globalID.type() == GeomType::Vertex ) {
        auto &x = d_mesh->getPos( d_globalID.elemID() );
        return Point( NP, x.data() );
    }
    ElementID ids[NG + 1];
    d_mesh->getElementsIDs( d_globalID.elemID(), GeomType::Vertex, ids );
    Point x( NP );
    for ( size_t i = 0; i <= NG; i++ ) {
        auto &x2 = d_mesh->getPos( ids[i] );
        for ( size_t d = 0; d < NP; d++ )
            x[d] += x2[d];
    }
    for ( size_t d = 0; d < NP; d++ )
        x[d] /= ( NG + 1 );
    return x;
}
template<size_t NG, size_t NP>
bool TriangleMeshElement<NG, NP>::containsPoint( const Point &pos, double TOL ) const
{
    NULL_USE( pos );
    NULL_USE( TOL );
    AMP_ERROR( "Not finished" );
    return false;
}
template<size_t NG, size_t NP>
bool TriangleMeshElement<NG, NP>::isOnSurface() const
{
    return d_mesh->isOnSurface( d_globalID.elemID() );
}
template<size_t NG, size_t NP>
bool TriangleMeshElement<NG, NP>::isOnBoundary( int id ) const
{
    return d_mesh->isOnBoundary( d_globalID.elemID(), id );
}
template<size_t NG, size_t NP>
bool TriangleMeshElement<NG, NP>::isInBlock( int id ) const
{
    return d_mesh->isInBlock( d_globalID.elemID(), id );
}


/********************************************************
 *  Explicit instantiations of TriangleMeshElement       *
 ********************************************************/
template class TriangleMeshElement<1, 1>;
template class TriangleMeshElement<1, 2>;
template class TriangleMeshElement<1, 3>;
template class TriangleMeshElement<2, 2>;
template class TriangleMeshElement<2, 3>;
template class TriangleMeshElement<3, 3>;


} // namespace Mesh
} // namespace AMP
