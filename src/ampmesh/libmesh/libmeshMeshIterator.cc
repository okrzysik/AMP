#include "AMP/ampmesh/libmesh/libmeshMeshIterator.h"
#include "AMP/ampmesh/libmesh/libmeshMeshElement.h"
#include "AMP/utils/Utilities.h"

// libMesh includes
#include "libmesh/elem.h"

namespace AMP {
namespace Mesh {


// unused global variable to prevent compiler warning
static MeshElement nullElement;


/********************************************************
 * Constructors                                          *
 ********************************************************/
libmeshMeshIterator::libmeshMeshIterator()
{
    d_typeID   = getTypeID();
    d_iterator = nullptr;
    d_begin2   = nullptr;
    d_end2     = nullptr;
    d_pos2     = nullptr;
    d_pos      = -1;
    d_type     = -1;
    d_size     = 0;
    d_rank     = 0;
    d_gcw      = 0;
    d_dim      = 0;
    d_mesh     = nullptr;
}
libmeshMeshIterator::libmeshMeshIterator( int type,
                                  const AMP::Mesh::libmeshMesh *mesh,
                                  int gcw,
                                  void *begin,
                                  void *end,
                                  void *pos,
                                  int size,
                                  int pos2 )
{
    d_typeID   = getTypeID();
    d_iterator = nullptr;
    d_type     = type;
    d_mesh     = mesh;
    d_gcw      = gcw;
    d_pos      = pos2;
    d_size     = size;
    d_rank     = mesh->getComm().getRank();
    d_meshID   = d_mesh->meshID();
    d_dim      = d_mesh->getlibMesh()->mesh_dimension();
    d_element  = &d_cur_element;
    if ( d_type == 0 ) {
        // Node iterator
        d_begin2 = (void *) new libMesh::Mesh::node_iterator( *( (libMesh::Mesh::node_iterator *) begin ) );
        d_end2   = (void *) new libMesh::Mesh::node_iterator( *( (libMesh::Mesh::node_iterator *) end ) );
        d_pos2   = (void *) new libMesh::Mesh::node_iterator( *( (libMesh::Mesh::node_iterator *) pos ) );
        // Count the number of nodes in the iterator
        if ( size == -1 ) {
            libMesh::Mesh::node_iterator cur_node =
                libMesh::Mesh::node_iterator( *( (libMesh::Mesh::node_iterator *) begin ) );
            libMesh::Mesh::node_iterator end_node =
                libMesh::Mesh::node_iterator( *( (libMesh::Mesh::node_iterator *) end ) );
            d_size = 0;
            while ( cur_node != end_node ) {
                d_size++;
                ++cur_node;
            }
        }
    } else if ( d_type == 1 ) {
        // Element iterator
        d_begin2 = (void *) new libMesh::Mesh::element_iterator( *( (libMesh::Mesh::element_iterator *) begin ) );
        d_end2   = (void *) new libMesh::Mesh::element_iterator( *( (libMesh::Mesh::element_iterator *) end ) );
        d_pos2   = (void *) new libMesh::Mesh::element_iterator( *( (libMesh::Mesh::element_iterator *) pos ) );
        // Count the number of elements in the iterator
        if ( size == -1 ) {
            libMesh::Mesh::element_iterator cur_element =
                libMesh::Mesh::element_iterator( *( (libMesh::Mesh::element_iterator *) begin ) );
            libMesh::Mesh::element_iterator end_element =
                libMesh::Mesh::element_iterator( *( (libMesh::Mesh::element_iterator *) end ) );
            d_size = 0;
            while ( cur_element != end_element ) {
                d_size++;
                ++cur_element;
            }
        }
    } else {
        AMP_ERROR( "libmeshMesh does not support iterators over this (unknown) type" );
    }
    // Count the position
    if ( pos2 == -1 ) {
        d_pos            = 0;
        MeshIterator tmp = this->begin();
        while ( this->operator!=( tmp ) ) {
            d_pos++;
            ++tmp;
        }
    }
    setCurrentElement();
}
libmeshMeshIterator::libmeshMeshIterator( const libmeshMeshIterator &rhs )
    : MeshIterator(), // Note: we never want to call the base copy constructor
      d_meshID( rhs.d_meshID )
{
    d_typeID   = getTypeID();
    d_iterator = nullptr;
    d_type     = rhs.d_type;
    d_mesh     = rhs.d_mesh;
    d_gcw      = rhs.d_gcw;
    d_pos      = rhs.d_pos;
    d_size     = rhs.d_size;
    d_rank     = rhs.d_rank;
    d_dim      = rhs.d_dim;
    d_element  = &d_cur_element;
    if ( d_type == 0 ) {
        // Node iterator
        d_begin2 =
            (void *) new libMesh::Mesh::node_iterator( *( (libMesh::Mesh::node_iterator *) rhs.d_begin2 ) );
        d_end2 = (void *) new libMesh::Mesh::node_iterator( *( (libMesh::Mesh::node_iterator *) rhs.d_end2 ) );
        d_pos2 = (void *) new libMesh::Mesh::node_iterator( *( (libMesh::Mesh::node_iterator *) rhs.d_pos2 ) );
    } else if ( d_type == 1 ) {
        // Element iterator
        d_begin2 =
            (void *) new libMesh::Mesh::element_iterator( *( (libMesh::Mesh::element_iterator *) rhs.d_begin2 ) );
        d_end2 =
            (void *) new libMesh::Mesh::element_iterator( *( (libMesh::Mesh::element_iterator *) rhs.d_end2 ) );
        d_pos2 =
            (void *) new libMesh::Mesh::element_iterator( *( (libMesh::Mesh::element_iterator *) rhs.d_pos2 ) );
    } else {
        AMP_ERROR( "libmeshMesh does not support iterators over this (unknown) type" );
    }
    setCurrentElement();
}
libmeshMeshIterator &libmeshMeshIterator::operator=( const libmeshMeshIterator &rhs )
{
    if ( this == &rhs ) // protect against invalid self-assignment
        return *this;
    this->d_typeID   = getTypeID();
    this->d_iterator = nullptr;
    this->d_type     = rhs.d_type;
    this->d_mesh     = rhs.d_mesh;
    this->d_gcw      = rhs.d_gcw;
    this->d_pos      = rhs.d_pos;
    this->d_size     = rhs.d_size;
    this->d_rank     = rhs.d_rank;
    this->d_meshID   = rhs.d_meshID;
    this->d_dim      = rhs.d_dim;
    this->d_element  = &d_cur_element;
    if ( this->d_type == 0 ) {
        // Node iterator
        this->d_begin2 =
            (void *) new libMesh::Mesh::node_iterator( *( (libMesh::Mesh::node_iterator *) rhs.d_begin2 ) );
        this->d_end2 =
            (void *) new libMesh::Mesh::node_iterator( *( (libMesh::Mesh::node_iterator *) rhs.d_end2 ) );
        this->d_pos2 =
            (void *) new libMesh::Mesh::node_iterator( *( (libMesh::Mesh::node_iterator *) rhs.d_pos2 ) );
    } else if ( this->d_type == 1 ) {
        // Element iterator
        this->d_begin2 =
            (void *) new libMesh::Mesh::element_iterator( *( (libMesh::Mesh::element_iterator *) rhs.d_begin2 ) );
        this->d_end2 =
            (void *) new libMesh::Mesh::element_iterator( *( (libMesh::Mesh::element_iterator *) rhs.d_end2 ) );
        this->d_pos2 =
            (void *) new libMesh::Mesh::element_iterator( *( (libMesh::Mesh::element_iterator *) rhs.d_pos2 ) );
    } else {
        AMP_ERROR( "libmeshMesh does not support iterators over this (unknown) type" );
    }
    setCurrentElement();
    return *this;
}


/********************************************************
 * Function to clone the iterator                        *
 ********************************************************/
MeshIterator *libmeshMeshIterator::clone() const { return new libmeshMeshIterator( *this ); }


/********************************************************
 * De-constructor                                        *
 ********************************************************/
libmeshMeshIterator::~libmeshMeshIterator()
{
    if ( d_pos2 != nullptr ) {
        if ( d_type == 0 ) {
            // Node iterator
            delete (libMesh::Mesh::node_iterator *) d_pos2;
            delete (libMesh::Mesh::node_iterator *) d_begin2;
            delete (libMesh::Mesh::node_iterator *) d_end2;
        } else if ( d_type == 1 ) {
            // Element iterator
            delete (libMesh::Mesh::element_iterator *) d_pos2;
            delete (libMesh::Mesh::element_iterator *) d_begin2;
            delete (libMesh::Mesh::element_iterator *) d_end2;
        }
    }
    d_iterator = nullptr;
}


/********************************************************
 * Return an iterator to the beginning or end            *
 ********************************************************/
MeshIterator libmeshMeshIterator::begin() const
{
    return libmeshMeshIterator( d_type, d_mesh, d_gcw, d_begin2, d_end2, d_begin2, d_size, 0 );
}
MeshIterator libmeshMeshIterator::end() const
{
    return libmeshMeshIterator( d_type, d_mesh, d_gcw, d_begin2, d_end2, d_end2, d_size, d_size );
}


/********************************************************
 * Increment/Decrement the iterator                      *
 ********************************************************/
MeshIterator &libmeshMeshIterator::operator++()
{
    // Prefix increment (increment and return this)
    d_pos++;
    if ( d_type == 0 ) {
        // Node iterator
        auto *it = (libMesh::Mesh::node_iterator *) d_pos2;
        it->operator++();
    } else if ( d_type == 1 ) {
        // Element iterator
        auto *it = (libMesh::Mesh::element_iterator *) d_pos2;
        it->operator++();
    } else {
        AMP_ERROR( "libmeshMesh does not support iterators over this (unknown) type" );
    }
    setCurrentElement();
    return *this;
}
MeshIterator libmeshMeshIterator::operator++( int )
{
    // Postfix increment (increment and return temporary object)
    libmeshMeshIterator tmp( *this ); // Create a temporary variable
    this->operator++();           // apply operator
    return tmp;                   // return temporary result
}
MeshIterator &libmeshMeshIterator::operator--()
{
    // Prefix decrement (increment and return this)
    AMP_ERROR( "Decrementing libmeshMesh iterators is not supported" );
    return *this;
}
MeshIterator libmeshMeshIterator::operator--( int )
{
    // Postfix decrement (increment and return temporary object)
    libmeshMeshIterator tmp( *this ); // Create a temporary variable
    --( *this );                  // apply operator
    return tmp;                   // return temporary result
}


/********************************************************
 * Random access incrementors                            *
 ********************************************************/
MeshIterator libmeshMeshIterator::operator+( int n ) const
{
    libmeshMeshIterator tmp( *this ); // Create a temporary iterator
    tmp.operator+=( n );          // Increment temporary iterator
    return tmp;
}
MeshIterator &libmeshMeshIterator::operator+=( int n )
{
    // Check the input
    if ( n >= 0 ) {
        if ( d_pos + n > d_size )
            AMP_ERROR( "Iterated past end of iterator" );
    } else { // decrement *this
        if ( -n > (int64_t) d_pos )
            AMP_ERROR( "Iterated past beginning of iterator" );
    }
    // Prform the increment and return
    if ( d_type == 0 ) {
        // Node iterator
        auto *it = (libMesh::Mesh::node_iterator *) d_pos2;
        if ( n >= 0 ) {
            for ( int i = 0; i < n; i++ )
                it->operator++();
        } else {
            AMP_ERROR( "Decrementing libmeshMesh iterators is not supported" );
        }
    } else if ( d_type == 1 ) {
        // Element iterator
        auto *it = (libMesh::Mesh::element_iterator *) d_pos2;
        if ( n >= 0 ) {
            for ( int i = 0; i < n; i++ )
                it->operator++();
        } else {
            AMP_ERROR( "Decrementing libmeshMesh iterators is not supported" );
        }
    } else {
        AMP_ERROR( "libmeshMesh does not support iterators over this (unknown) type" );
    }
    d_pos += n;
    setCurrentElement();
    return *this;
}


/********************************************************
 * Compare two iterators                                 *
 ********************************************************/
bool libmeshMeshIterator::operator==( const MeshIterator &rhs ) const
{
    const libmeshMeshIterator *rhs2 = nullptr;
    // Convert rhs to a libmeshMeshIterator* so we can access the base class members
    auto *tmp = reinterpret_cast<const libmeshMeshIterator *>( &rhs );
    if ( tmp->d_typeID == getTypeID() ) {
        rhs2 = tmp; // We can safely cast rhs to a libmeshMeshIterator
    } else if ( tmp->d_iterator != nullptr ) {
        tmp = reinterpret_cast<const libmeshMeshIterator *>( tmp->d_iterator );
        if ( tmp->d_typeID == getTypeID() )
            rhs2 = tmp; // We can safely cast rhs.iterator to a libmeshMeshIterator
    }
    // Perform direct comparisions if we are dealing with two libmeshMeshIterators
    if ( rhs2 != nullptr ) {
        if ( d_type != rhs2->d_type ) {
            // We are dealing with different types of elements
            return false;
        }
        if ( d_type == 0 ) {
            // Node iterator
            return ( *( (libMesh::Mesh::node_iterator *) d_pos2 ) ) ==
                   ( *( (libMesh::Mesh::node_iterator *) rhs2->d_pos2 ) );
        } else if ( d_type == 1 ) {
            // Element iterator
            return ( *( (libMesh::Mesh::element_iterator *) d_pos2 ) ) ==
                   ( *( (libMesh::Mesh::element_iterator *) rhs2->d_pos2 ) );
        } else {
            AMP_ERROR( "libmeshMesh does not support iterators over this (unknown) type" );
        }
    }
    /* We are comparing a libmeshMeshIterator to an arbitrary iterator
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
bool libmeshMeshIterator::operator!=( const MeshIterator &rhs ) const { return !( ( *this ) == rhs ); }


/********************************************************
 * Dereference the iterator to get the element           *
 ********************************************************/
void libmeshMeshIterator::setCurrentElement()
{
    if ( d_pos >= d_size ) {
        d_cur_element = libmeshMeshElement();
    } else if ( d_type == 0 ) {
        // Node iterator
        auto *it     = (libMesh::Mesh::node_iterator *) d_pos2;
        libMesh::Node *node = it->operator*();
        d_cur_element =
            libmeshMeshElement( d_dim, GeomType::Vertex, (void *) node, d_rank, d_meshID, d_mesh );
    } else if ( d_type == 1 ) {
        // Element iterator
        auto *it     = (libMesh::Mesh::element_iterator *) d_pos2;
        libMesh::Elem *elem = it->operator*();
        d_cur_element =
            libmeshMeshElement( d_dim, (GeomType) d_dim, (void *) elem, d_rank, d_meshID, d_mesh );
    } else {
        AMP_ERROR( "libmeshMesh does not support iterators over this (unknown) type" );
    }
}


} // namespace Mesh
} // namespace AMP
