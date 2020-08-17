#include "AMP/ampmesh/triangle/TriangleMeshElement.h"
#include "AMP/ampmesh/MeshPoint.h"
#include "AMP/ampmesh/shapes/GeometryHelpers.h"
#include "AMP/ampmesh/triangle/TriangleMesh.h"
#include "AMP/ampmesh/triangle/TriangleMeshIterator.h"
#include "AMP/utils/DelaunayHelpers.h"
#include "AMP/utils/Utilities.h"

#include <limits>


/****************************************************************
 * Overload basic operations                                     *
 ****************************************************************/
template<size_t N>
static inline std::array<double, N> operator*( double x, const std::array<double, N> &y )
{
    if constexpr ( N == 1 )
        return { x * y[0] };
    else if constexpr ( N == 2 )
        return { x * y[0], x * y[1] };
    else if constexpr ( N == 3 )
        return { x * y[0], x * y[1], x * y[2] };
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
static constexpr double inv_factorial( int N )
{
    double x = 1;
    for ( int i = 2; i <= N; i++ )
        x *= i;
    return 1.0 / x;
}
template<size_t N>
static inline double abs( const std::array<double, N> &x )
{
    if constexpr ( N == 1 )
        return std::abs( x[0] );
    else if constexpr ( N == 2 )
        return sqrt( x[0] * x[0] + x[1] * x[1] );
    else if constexpr ( N == 3 )
        return sqrt( x[0] * x[0] + x[1] * x[1] + x[2] * x[2] );
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
static inline std::array<double, 3> cross( const std::array<double, 3> &x,
                                           const std::array<double, 3> &y )
{
    return { x[1] * y[2] - x[2] * y[1], x[2] * y[0] - x[0] * y[2], x[0] * y[1] - x[1] * y[0] };
}
template<size_t N>
static inline std::array<double, N> normalize( const std::array<double, N> &x )
{
    if constexpr ( N == 1 ) {
        return { 1.0 };
    } else if constexpr ( N == 2 ) {
        double tmp = 1.0 / sqrt( x[0] * x[0] + x[1] * x[1] );
        return { tmp * x[0], tmp * x[1] };
    } else if constexpr ( N == 3 ) {
        double tmp = 1.0 / sqrt( x[0] * x[0] + x[1] * x[1] + x[2] * x[2] );
        return { tmp * x[0], tmp * x[1], tmp * x[2] };
    }
}


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
template<uint8_t NG, uint8_t NP, uint8_t TYPE>
constexpr uint32_t TriangleMeshElement<NG, NP, TYPE>::getTypeID()
{
    char name[] = "TriangleMeshElement<0,0,0>";
    name[21]    = 48 + NG;
    name[23]    = 48 + NP;
    name[25]    = 48 + TYPE;
    return AMP::Utilities::hash_char( name );
}


/********************************************************
 * Constructors                                          *
 ********************************************************/
template<uint8_t NG, uint8_t NP, uint8_t TYPE>
TriangleMeshElement<NG, NP, TYPE>::TriangleMeshElement()
{
    typeID  = getTypeID();
    element = nullptr;
    d_mesh  = nullptr;
}
template<uint8_t NG, uint8_t NP, uint8_t TYPE>
TriangleMeshElement<NG, NP, TYPE>::TriangleMeshElement( const MeshElementID &id,
                                                        const TriangleMesh<NG, NP> *mesh )
{
    typeID     = getTypeID();
    element    = nullptr;
    d_globalID = id;
    d_mesh     = mesh;
#if ( defined( DEBUG ) || defined( _DEBUG ) ) && !defined( NDEBUG )
    auto type = static_cast<uint8_t>( id.type() );
    if ( type != TYPE && type != 255 )
        printf( "%i %i %i %i\n", NG, NP, TYPE, type );
    AMP_ASSERT( type == TYPE || type == 255 );
#endif
}
template<uint8_t NG, uint8_t NP, uint8_t TYPE>
TriangleMeshElement<NG, NP, TYPE>::TriangleMeshElement( const TriangleMeshElement &rhs )
    : MeshElement(), d_mesh( rhs.d_mesh ), d_globalID( rhs.d_globalID )
{
    typeID  = getTypeID();
    element = rhs.element;
}
template<uint8_t NG, uint8_t NP, uint8_t TYPE>
TriangleMeshElement<NG, NP, TYPE>::TriangleMeshElement( TriangleMeshElement &&rhs )
    : MeshElement(), d_mesh( rhs.d_mesh ), d_globalID{ rhs.d_globalID }
{
    typeID  = rhs.typeID;
    element = nullptr;
}
template<uint8_t NG, uint8_t NP, uint8_t TYPE>
TriangleMeshElement<NG, NP, TYPE> &TriangleMeshElement<NG, NP, TYPE>::
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
template<uint8_t NG, uint8_t NP, uint8_t TYPE>
TriangleMeshElement<NG, NP, TYPE> &TriangleMeshElement<NG, NP, TYPE>::
operator=( TriangleMeshElement &&rhs )
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
template<uint8_t NG, uint8_t NP, uint8_t TYPE>
MeshElement *TriangleMeshElement<NG, NP, TYPE>::clone() const
{
    return new TriangleMeshElement<NG, NP, TYPE>( *this );
}


/****************************************************************
 * Return the global rank of the owner rank                      *
 ****************************************************************/
template<uint8_t NG, uint8_t NP, uint8_t TYPE>
unsigned int TriangleMeshElement<NG, NP, TYPE>::globalOwnerRank() const
{
    return d_mesh->getComm().globalRanks()[d_globalID.owner_rank()];
}


/****************************************************************
 * Function to get the elements composing the current element    *
 ****************************************************************/
template<uint8_t NG, uint8_t NP, uint8_t TYPE>
void TriangleMeshElement<NG, NP, TYPE>::getElementsID( const GeomType type,
                                                       std::vector<MeshElementID> &ID ) const
{
    // Number of elements composing a given type
    int N = n_Simplex_elements[TYPE][static_cast<uint8_t>( type )];
    // Get the element ids
    ElementID tmp[6];
    d_mesh->getElementsIDs( d_globalID.elemID(), type, tmp );
    ID.resize( N );
    for ( int i = 0; i < N; i++ )
        ID[i] = MeshElementID( d_globalID.meshID(), tmp[i] );
}
template<uint8_t NG, uint8_t NP, uint8_t TYPE>
void TriangleMeshElement<NG, NP, TYPE>::getElements( const GeomType type,
                                                     std::vector<MeshElement> &children ) const
{
    // Number of elements composing a given type
    int N = n_Simplex_elements[TYPE][static_cast<uint8_t>( type )];
    // Get the element ids
    ElementID tmp[6];
    d_mesh->getElementsIDs( d_globalID.elemID(), type, tmp );
    // Create the mesh elements
    auto meshID = d_globalID.meshID();
    children.resize( N );
    for ( int i = 0; i < N; i++ )
        children[i] = d_mesh->getElement( MeshElementID( meshID, tmp[i] ) );
}


/****************************************************************
 * Function to get the neighboring elements                      *
 ****************************************************************/
template<uint8_t NG, uint8_t NP, uint8_t TYPE>
void TriangleMeshElement<NG, NP, TYPE>::getNeighbors(
    std::vector<MeshElement::shared_ptr> &neighbors ) const
{
    std::vector<ElementID> neighborIDs;
    d_mesh->getNeighborIDs( d_globalID.elemID(), neighborIDs );
    neighbors.resize( neighborIDs.size() );
    auto meshID = d_globalID.meshID();
    for ( size_t i = 0; i < neighborIDs.size(); i++ )
        neighbors[i].reset( d_mesh->getElement2( MeshElementID( meshID, neighborIDs[i] ) ) );
}


/****************************************************************
 * Get the coordinates of the verticies                          *
 ****************************************************************/
template<uint8_t NG, uint8_t NP, uint8_t TYPE>
inline void TriangleMeshElement<NG, NP, TYPE>::getVertexCoord( std::array<double, NP> *x ) const
{
    if constexpr ( TYPE == 0 ) {
        x[0] == d_mesh->getPos( d_globalID.elemID() );
    } else {
        ElementID ids[TYPE + 1];
        d_mesh->getElementsIDs( d_globalID.elemID(), GeomType::Vertex, ids );
        for ( size_t i = 0; i <= TYPE; i++ )
            x[i] = d_mesh->getPos( ids[i] );
    }
}


/****************************************************************
 * Functions to get basic element properties                     *
 ****************************************************************/
template<uint8_t NG, uint8_t NP, uint8_t TYPE>
double TriangleMeshElement<NG, NP, TYPE>::volume() const
{
    std::array<double, NP> x[TYPE + 1];
    getVertexCoord( x );
    if constexpr ( TYPE == 0 ) {
        return 0;
    } else if constexpr ( TYPE == 1 ) {
        return abs( x[1] - x[0] );
    } else if constexpr ( TYPE == 2 ) {
        auto AB  = x[1] - x[0];
        auto AC  = x[2] - x[0];
        double t = dot( AB, AC );
        return 0.5 * sqrt( dot( AB, AB ) * dot( AC, AC ) - t * t );
    } else if constexpr ( TYPE == NP ) {
        /* Calculate the volume of a N-dimensional simplex:
         *         1  |  x1-x4   x2-x4   x3-x4  |
         *    V = --  |  y1-y4   y2-y4   y3-y4  |   (3D)
         *        n!  |  z1-z4   z2-z4   z3-z4  |
         * Note: the sign of the volume depends on the order of the points.
         *   It will be positive for points stored in a clockwise manner
         * Note:  If the volume is zero, then the simplex is invalid
         *   Eg. a line in 2D or a plane in 3D.             */
        double M[TYPE * TYPE];
        for ( size_t i = 0; i < TYPE; i++ ) {
            for ( size_t d = 0; d < TYPE; d++ )
                M[d + i * TYPE] = x[i][d] - x[TYPE][d];
        }
        constexpr double C = inv_factorial( TYPE );
        double V           = std::abs( C * DelaunayHelpers<TYPE>::det( M ) );
        return V;
    } else {
        AMP_ERROR( "Not finished" );
        return 0;
    }
}
template<uint8_t NG, uint8_t NP, uint8_t TYPE>
MeshPoint<double> TriangleMeshElement<NG, NP, TYPE>::norm() const
{
    std::array<double, NP> x[TYPE + 1];
    getVertexCoord( x );
    if constexpr ( TYPE == 2 && NP == 3 ) {
        auto n = AMP::Geometry::GeometryHelpers::normal( x[0], x[1], x[2] );
        return { n[0], n[1], n[2] };
    } else {
        AMP_ERROR( "Not finished" );
    }
    return MeshPoint<double>();
}
template<uint8_t NG, uint8_t NP, uint8_t TYPE>
MeshPoint<double> TriangleMeshElement<NG, NP, TYPE>::coord() const
{
    if constexpr ( TYPE == 0 ) {
        auto x = d_mesh->getPos( d_globalID.elemID() );
        return MeshPoint<double>( NP, x.data() );
    } else {
        AMP_ERROR( "coord is only valid for verticies: " + std::to_string( (int) TYPE ) );
        return MeshPoint<double>();
    }
}
template<uint8_t NG, uint8_t NP, uint8_t TYPE>
MeshPoint<double> TriangleMeshElement<NG, NP, TYPE>::centroid() const
{
    if constexpr ( TYPE == 0 )
        return MeshPoint<double>( d_mesh->getPos( d_globalID.elemID() ) );
    std::array<double, NP> x2[TYPE + 1];
    getVertexCoord( x2 );
    ElementID ids[TYPE + 1];
    d_mesh->getElementsIDs( d_globalID.elemID(), GeomType::Vertex, ids );
    MeshPoint<double> x( (size_t) NP );
    for ( size_t i = 0; i <= TYPE; i++ ) {
        for ( size_t d = 0; d < NP; d++ )
            x[d] += x2[i][d];
    }
    for ( size_t d = 0; d < NP; d++ )
        x[d] /= ( TYPE + 1 );
    return x;
}
template<size_t N>
std::array<double, N> convert( const MeshPoint<double> &p )
{
    std::array<double, N> p2;
    for ( size_t i = 0; i < N; i++ )
        p2[i] = p[i];
    return p2;
}
template<uint8_t NG, uint8_t NP, uint8_t TYPE>
bool TriangleMeshElement<NG, NP, TYPE>::containsPoint( const MeshPoint<double> &pos,
                                                       double TOL ) const
{
    // Get the vertex coordinates
    std::array<double, NP> x[TYPE + 1];
    getVertexCoord( x );
    // Check if the point is in the triangle
    if constexpr ( TYPE == 2 && NP == 3 ) {
        // Compute barycentric coordinates
        auto L =
            AMP::Geometry::GeometryHelpers::barycentric<3, 3>( x, { pos.x(), pos.y(), pos.z() } );
        return ( L[0] >= -TOL ) && ( L[1] >= -TOL ) && ( L[2] >= -TOL );
    } else {
        AMP_ERROR( "Not finished" );
    }
    return false;
}
template<uint8_t NG, uint8_t NP, uint8_t TYPE>
bool TriangleMeshElement<NG, NP, TYPE>::isOnSurface() const
{
    return d_mesh->isOnSurface( d_globalID.elemID() );
}
template<uint8_t NG, uint8_t NP, uint8_t TYPE>
bool TriangleMeshElement<NG, NP, TYPE>::isOnBoundary( int id ) const
{
    return d_mesh->isOnBoundary( d_globalID.elemID(), id );
}
template<uint8_t NG, uint8_t NP, uint8_t TYPE>
bool TriangleMeshElement<NG, NP, TYPE>::isInBlock( int id ) const
{
    return d_mesh->isInBlock( d_globalID.elemID(), id );
}


/****************************************************************
 * Calculate the nearest point on the element                    *
 ****************************************************************/
template<uint8_t NG, uint8_t NP, uint8_t TYPE>
MeshPoint<double> TriangleMeshElement<NG, NP, TYPE>::nearest( const MeshPoint<double> &pos ) const
{
    // Get the vertex coordinates
    std::array<double, NP> v[TYPE + 1];
    getVertexCoord( v );
    if constexpr ( TYPE == 2 && NP == 3 ) {
        auto p = AMP::Geometry::GeometryHelpers::nearest( v, { pos.x(), pos.y(), pos.z() } );
        return { p[0], p[1], p[2] };
    } else {
        AMP_ERROR( "Not finished" );
    }
    return MeshPoint<double>();
}


/****************************************************************
 * Calculate the distance to the element                         *
 ****************************************************************/
template<uint8_t NG, uint8_t NP, uint8_t TYPE>
double TriangleMeshElement<NG, NP, TYPE>::distance( const MeshPoint<double> &pos,
                                                    const MeshPoint<double> &dir ) const
{
    // Get the vertex coordinates
    std::array<double, NP> x[TYPE + 1];
    getVertexCoord( x );
    if constexpr ( TYPE == 2 && NP == 3 ) {
        // Get the normal and a point on the plane containing the triangle
        auto n = AMP::Geometry::GeometryHelpers::normal( x[0], x[1], x[2] );
        // Find the point of intersection with the line and the plane
        std::array<double, 3> l0 = { pos[0], pos[1], pos[2] };
        std::array<double, 3> l  = { dir[0], dir[1], dir[2] };
        double d                 = dot( x[0] - l0, n ) / dot( l, n );
        if ( fabs( d ) < 1e-8 ) {
            // Point is either in the element or the line and plane do not intersect
            if ( containsPoint( pos, 1e-8 ) )
                return 0;
            else
                return std::numeric_limits<double>::infinity();
        }
        if ( d < 0 ) {
            // Ray is pointed away from element
            return std::numeric_limits<double>::infinity();
        }
        // Calculate point of intersection and check if it is in the element
        auto p2 = pos + d * dir;
        if ( containsPoint( { p2[0], p2[1], p2[2] }, 1e-8 ) )
            return d;
        else
            return std::numeric_limits<double>::infinity();
    } else {
        AMP_ERROR( "Not finished" );
    }
    return 0;
}


/********************************************************
 *  Explicit instantiations of TriangleMeshElement       *
 ********************************************************/
template class TriangleMeshElement<1, 1, 0>;
template class TriangleMeshElement<1, 1, 1>;
template class TriangleMeshElement<1, 2, 0>;
template class TriangleMeshElement<1, 2, 1>;
template class TriangleMeshElement<1, 3, 0>;
template class TriangleMeshElement<1, 3, 1>;
template class TriangleMeshElement<2, 2, 0>;
template class TriangleMeshElement<2, 2, 1>;
template class TriangleMeshElement<2, 2, 2>;
template class TriangleMeshElement<2, 3, 0>;
template class TriangleMeshElement<2, 3, 1>;
template class TriangleMeshElement<2, 3, 2>;
template class TriangleMeshElement<3, 3, 0>;
template class TriangleMeshElement<3, 3, 1>;
template class TriangleMeshElement<3, 3, 2>;
template class TriangleMeshElement<3, 3, 3>;


} // namespace Mesh
} // namespace AMP
