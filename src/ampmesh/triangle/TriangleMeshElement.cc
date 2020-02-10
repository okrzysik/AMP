#include "AMP/ampmesh/triangle/TriangleMeshElement.h"
#include "AMP/ampmesh/triangle/TriangleMesh.h"
#include "AMP/ampmesh/triangle/TriangleMeshIterator.h"
#include "AMP/utils/DelaunayHelpers.h"
#include "AMP/utils/Utilities.h"


namespace AMP {
namespace Mesh {


/****************************************************************
 * Get the number of n-Simplex elements of each type             *
 ****************************************************************/
// clang-format off
static constexpr uint8_t n_Simplex_elements[4][4] = {
    {  1, 0, 0, 0 },
    {  2, 1, 0, 0 },
    {  3, 3, 1, 0 },
    {  4, 6, 4, 1 },
};
// clang-format on


/********************************************************
 * Create a unique id for each class                     *
 ********************************************************/
template<size_t NG, size_t NP>
constexpr uint32_t TriangleMeshElement<NG, NP>::getTypeID()
{
    static_assert( NG > 0 && NP <= 3 && NP >= NG );
    if constexpr ( NG == 1 && NP == 1 )
        return AMP::Utilities::hash_char( "TriangleMeshElement<1,1>" );
    else if constexpr ( NG == 1 && NP == 2 )
        return AMP::Utilities::hash_char( "TriangleMeshElement<1,2>" );
    else if constexpr ( NG == 1 && NP == 3 )
        return AMP::Utilities::hash_char( "TriangleMeshElement<1,3>" );
    else if constexpr ( NG == 2 && NP == 2 )
        return AMP::Utilities::hash_char( "TriangleMeshElement<2,2>" );
    else if constexpr ( NG == 2 && NP == 3 )
        return AMP::Utilities::hash_char( "TriangleMeshElement<2,3>" );
    else if constexpr ( NG == 3 && NP == 3 )
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
static constexpr int get_N_elements( const GeomType &src, const GeomType &dst )
{
    return n_Simplex_elements[static_cast<int>( src )][static_cast<int>( dst )];
}
template<size_t NG, size_t NP>
void TriangleMeshElement<NG, NP>::getElementsID( const GeomType type,
                                                 std::vector<MeshElementID> &ID ) const
{
    // Number of elements composing a given type
    int N = get_N_elements( d_globalID.type(), type );
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
    // Number of elements composing a given type
    int N = get_N_elements( d_globalID.type(), type );
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
static constexpr double inv_factorial( int N )
{
    double x = 1;
    for ( int i = 2; i <= N; i++ )
        x *= i;
    return 1.0 / x;
}
template<size_t N>
static inline double dot( const std::array<double, N> &x, const std::array<double, N> &y )
{
    if constexpr ( N == 1 )
        return x[0] * y[0];
    else if constexpr ( N == 2 )
        return x[0] * y[0] + x[1] * y[1];
    else if constexpr ( N == 3 )
        return x[0] * y[0] + x[1] * y[1] + x[2] * y[2];
}
template<size_t N>
static inline std::array<double, N> operator-( const std::array<double, N> &x,
                                               const std::array<double, N> &y )
{
    if constexpr ( N == 1 )
        return { x[0] - y[0] };
    else if constexpr ( N == 2 )
        return { x[0] - y[0], x[1] - y[1] };
    else if constexpr ( N == 3 )
        return { x[0] - y[0], x[1] - y[1], x[2] - y[2] };
}
template<size_t NG, size_t NP>
double TriangleMeshElement<NG, NP>::volume() const
{
    if constexpr ( NG == 0 ) {
        return 0;
    } else if constexpr ( NG == 1 ) {
        ElementID ids[2];
        d_mesh->getElementsIDs( d_globalID.elemID(), GeomType::Vertex, ids );
        auto p1     = d_mesh->getPos( ids[0] );
        auto p2     = d_mesh->getPos( ids[1] );
        double dist = 0;
        for ( size_t d = 0; d < NP; d++ )
            dist += ( p1[d] - p2[d] ) * ( p1[d] - p2[d] );
        return dist;
    } else if constexpr ( NG == 2 ) {
        ElementID ids[3];
        d_mesh->getElementsIDs( d_globalID.elemID(), GeomType::Vertex, ids );
        auto p0  = d_mesh->getPos( ids[0] );
        auto AB  = p0 - d_mesh->getPos( ids[1] );
        auto AC  = p0 - d_mesh->getPos( ids[2] );
        double t = dot( AB, AC );
        return 0.5 * sqrt( dot( AB, AB ) * dot( AC, AC ) - t * t );
    } else if constexpr ( NG == NP ) {
        /* Calculate the volume of a N-dimensional simplex:
         *         1  |  x1-x4   x2-x4   x3-x4  |
         *    V = --  |  y1-y4   y2-y4   y3-y4  |   (3D)
         *        n!  |  z1-z4   z2-z4   z3-z4  |
         * Note: the sign of the volume depends on the order of the points.
         *   It will be positive for points stored in a clockwise manner
         * Note:  If the volume is zero, then the simplex is invalid
         *   Eg. a line in 2D or a plane in 3D.             */
        ElementID ids[NG + 1];
        d_mesh->getElementsIDs( d_globalID.elemID(), GeomType::Vertex, ids );
        double M[NG * NG];
        auto p0 = d_mesh->getPos( ids[NG] );
        for ( size_t i = 0; i < NG; i++ ) {
            auto p = d_mesh->getPos( ids[i] );
            for ( size_t d = 0; d < NG; d++ )
                M[d + i * NG] = p[d] - p0[d];
        }
        constexpr double C = inv_factorial( NG );
        return std::abs( C * DelaunayHelpers<NG>::det( M ) );
    } else {
        AMP_ERROR( "Not finished" );
        return 0;
    }
}
template<size_t NG, size_t NP>
Point TriangleMeshElement<NG, NP>::norm() const
{
    AMP_ERROR( "norm not implimented yet" );
    return Point();
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
