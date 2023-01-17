#include "AMP/mesh/libmesh/libmeshNodeIterator.h"
#include "AMP/mesh/libmesh/libmeshMeshElement.h"

// libMesh includes
#include "libmesh/elem.h"

namespace AMP::Mesh {


// unused global variable to prevent compiler warning
static MeshElement nullElement;


/********************************************************
 * Constructors                                          *
 ********************************************************/
libmeshNodeIterator::libmeshNodeIterator( const AMP::Mesh::libmeshMesh *mesh,
                                          int gcw,
                                          const libMesh::Mesh::node_iterator &begin,
                                          const libMesh::Mesh::node_iterator &end,
                                          const libMesh::Mesh::node_iterator &pos,
                                          int size,
                                          int pos2 )
    : d_gcw( gcw ),
      d_dim( mesh->getlibMesh()->mesh_dimension() ),
      d_rank( mesh->getComm().getRank() ),
      d_begin2( begin ),
      d_end2( end ),
      d_pos2( pos ),
      d_meshID( mesh->meshID() ),
      d_mesh( mesh )
{
    d_typeID       = getTypeID<decltype( *this )>();
    d_iteratorType = MeshIterator::Type::Forward;
    d_pos          = pos2;
    d_size         = size;
    d_element      = &d_cur_element;
    // Count the number of elements in the iterator
    if ( size == -1 ) {
        d_size = 0;
        libMesh::Mesh::node_iterator cur( begin );
        while ( cur != d_end2 ) {
            d_size++;
            ++cur;
        }
    }
    // Count the position
    if ( pos2 == -1 ) {
        d_pos = 0;
        libMesh::Mesh::node_iterator cur( begin );
        while ( cur != pos ) {
            d_pos++;
            ++cur;
        }
    }
    setCurrentElement();
}
libmeshNodeIterator::libmeshNodeIterator( const libmeshNodeIterator &rhs )
    : MeshIterator(), // Note: we never want to call the base copy constructor
      d_gcw( rhs.d_gcw ),
      d_dim( rhs.d_dim ),
      d_rank( rhs.d_rank ),
      d_begin2( rhs.d_begin2 ),
      d_end2( rhs.d_end2 ),
      d_pos2( rhs.d_pos2 ),
      d_meshID( rhs.d_meshID ),
      d_mesh( rhs.d_mesh )
{
    d_typeID       = rhs.d_typeID;
    d_iteratorType = rhs.d_iteratorType;
    d_pos          = rhs.d_pos;
    d_size         = rhs.d_size;
    d_element      = &d_cur_element;
    setCurrentElement();
}
libmeshNodeIterator &libmeshNodeIterator::operator=( const libmeshNodeIterator &rhs )
{
    if ( this == &rhs ) // protect against invalid self-assignment
        return *this;
    this->d_iterator     = nullptr;
    this->d_typeID       = getTypeID<decltype( *this )>();
    this->d_iteratorType = rhs.d_iteratorType;
    this->d_mesh         = rhs.d_mesh;
    this->d_gcw          = rhs.d_gcw;
    this->d_pos          = rhs.d_pos;
    this->d_size         = rhs.d_size;
    this->d_rank         = rhs.d_rank;
    this->d_meshID       = rhs.d_meshID;
    this->d_dim          = rhs.d_dim;
    this->d_element      = &d_cur_element;
    this->d_begin2       = rhs.d_begin2;
    this->d_end2         = rhs.d_end2;
    this->d_pos2         = rhs.d_pos2;
    setCurrentElement();
    return *this;
}


/********************************************************
 * Function to clone the iterator                        *
 ********************************************************/
MeshIterator *libmeshNodeIterator::clone() const { return new libmeshNodeIterator( *this ); }


/********************************************************
 * Return an iterator to the beginning or end            *
 ********************************************************/
MeshIterator libmeshNodeIterator::begin() const
{
    return libmeshNodeIterator( d_mesh, d_gcw, d_begin2, d_end2, d_begin2, d_size, 0 );
}
MeshIterator libmeshNodeIterator::end() const
{
    return libmeshNodeIterator( d_mesh, d_gcw, d_begin2, d_end2, d_end2, d_size, d_size );
}


/********************************************************
 * Increment/Decrement the iterator                      *
 ********************************************************/
MeshIterator &libmeshNodeIterator::operator++()
{
    // Prefix increment (increment and return this)
    ++d_pos;
    ++d_pos2;
    setCurrentElement();
    return *this;
}
MeshIterator libmeshNodeIterator::operator++( int )
{
    // Postfix increment (increment and return temporary object)
    libmeshNodeIterator tmp( *this ); // Create a temporary variable
    this->operator++();               // apply operator
    return tmp;                       // return temporary result
}
MeshIterator &libmeshNodeIterator::operator--()
{
    // Prefix decrement (decrement and return this)
    AMP_ERROR( "Decrementing libmeshMesh iterators is not supported" );
    return *this;
}
MeshIterator libmeshNodeIterator::operator--( int )
{
    // Postfix decrement (increment and return temporary object)
    libmeshNodeIterator tmp( *this ); // Create a temporary variable
    --( *this );                      // apply operator
    return tmp;                       // return temporary result
}


/********************************************************
 * Random access incrementors                            *
 ********************************************************/
MeshIterator libmeshNodeIterator::operator+( int n ) const
{
    libmeshNodeIterator tmp( *this ); // Create a temporary iterator
    tmp.operator+=( n );              // Increment temporary iterator
    return tmp;                       // return temporary result
}
MeshIterator &libmeshNodeIterator::operator+=( int n )
{
    // Check the input
    if ( n >= 0 ) {
        if ( d_pos + n > d_size )
            AMP_ERROR( "Iterated past end of iterator" );
    } else { // decrement *this
        if ( -n > (int64_t) d_pos )
            AMP_ERROR( "Iterated past beginning of iterator" );
    }
    // Perform the increment and return
    if ( n >= 0 ) {
        for ( int i = 0; i < n; i++ )
            ++d_pos2;
    } else {
        AMP_ERROR( "Decrementing libmeshMesh iterators is not supported" );
    }
    d_pos += n;
    setCurrentElement();
    return *this;
}


/********************************************************
 * Compare two iterators                                 *
 ********************************************************/
bool libmeshNodeIterator::operator==( const MeshIterator &rhs ) const
{
    const libmeshNodeIterator *rhs2 = nullptr;
    // Convert rhs to a libmeshNodeIterator* so we can access the base class members
    auto *tmp = reinterpret_cast<const libmeshNodeIterator *>( &rhs );
    if ( tmp->d_typeID == getTypeID<decltype( *this )>() ) {
        rhs2 = tmp; // We can safely cast rhs to a libmeshNodeIterator
    } else if ( tmp->d_iterator != nullptr ) {
        tmp = reinterpret_cast<const libmeshNodeIterator *>( tmp->d_iterator );
        if ( tmp->d_typeID == getTypeID<decltype( *this )>() )
            rhs2 = tmp; // We can safely cast rhs.iterator to a libmeshNodeIterator
    }
    // Perform direct comparisions if we are dealing with two libmeshNodeIterators;
    if ( rhs2 != nullptr )
        return d_pos2 == rhs2->d_pos2;
    /* We are comparing a libmeshNodeIterator to an arbitrary iterator
     * The iterators are the same if they point to the same position and iterate
     * over the same elements in the same order
     */
    // Check the size
    if ( this->size() != rhs.size() )
        return false;
    // Check the current position
    if ( this->position() != rhs.position() )
        return false;
    // Check that the elements match
    MeshIterator it1    = this->begin();
    MeshIterator it2    = rhs.begin();
    bool elements_match = true;
    for ( size_t i = 0; i < it1.size(); ++i, ++it1, ++it2 ) {
        if ( it1->globalID() != it2->globalID() )
            elements_match = false;
    }
    return elements_match;
}
bool libmeshNodeIterator::operator!=( const MeshIterator &rhs ) const
{
    return !( ( *this ) == rhs );
}


/********************************************************
 * Dereference the iterator to get the element           *
 ********************************************************/
void libmeshNodeIterator::setCurrentElement()
{
    if ( d_pos >= d_size )
        d_cur_element = libmeshMeshElement();
    else
        d_cur_element =
            libmeshMeshElement( d_dim, GeomType::Vertex, *d_pos2, d_rank, d_meshID, d_mesh );
}


} // namespace AMP::Mesh
