#include "ampmesh/libmesh/libMeshIterator.h"
#include "ampmesh/libmesh/libMeshElement.h"
#include "utils/Utilities.h"

// libMesh includes
#include "libmesh/elem.h"

namespace AMP {
namespace Mesh {


// Create a unique id for this class
static unsigned int libMeshIteratorTypeID = TYPE_HASH( libMeshIterator );

// unused global variable to prevent compiler warning
static MeshElement nullElement;


/********************************************************
* Constructors                                          *
********************************************************/
libMeshIterator::libMeshIterator()
{
    d_typeID   = libMeshIteratorTypeID;
    d_iterator = nullptr;
    d_begin2   = nullptr;
    d_end2     = nullptr;
    d_pos2     = nullptr;
    d_pos      = -1;
    d_type     = -1;
    d_size     = 0;
    d_rank     = 0;
    d_element  = &d_cur_element;
}
libMeshIterator::libMeshIterator( int type,
                                  const AMP::Mesh::libMesh *mesh,
                                  int gcw,
                                  void *begin,
                                  void *end,
                                  void *pos,
                                  int size,
                                  int pos2 )
{
    d_typeID   = libMeshIteratorTypeID;
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
        d_begin2 = (void *) new ::Mesh::node_iterator( *( (::Mesh::node_iterator *) begin ) );
        d_end2   = (void *) new ::Mesh::node_iterator( *( (::Mesh::node_iterator *) end ) );
        d_pos2   = (void *) new ::Mesh::node_iterator( *( (::Mesh::node_iterator *) pos ) );
        // Count the number of nodes in the iterator
        if ( size == -1 ) {
            ::Mesh::node_iterator cur_node =
                ::Mesh::node_iterator( *( (::Mesh::node_iterator *) begin ) );
            ::Mesh::node_iterator end_node =
                ::Mesh::node_iterator( *( (::Mesh::node_iterator *) end ) );
            d_size = 0;
            while ( cur_node != end_node ) {
                d_size++;
                ++cur_node;
            }
        }
    } else if ( d_type == 1 ) {
        // Element iterator
        d_begin2 = (void *) new ::Mesh::element_iterator( *( (::Mesh::element_iterator *) begin ) );
        d_end2   = (void *) new ::Mesh::element_iterator( *( (::Mesh::element_iterator *) end ) );
        d_pos2   = (void *) new ::Mesh::element_iterator( *( (::Mesh::element_iterator *) pos ) );
        // Count the number of elements in the iterator
        if ( size == -1 ) {
            ::Mesh::element_iterator cur_element =
                ::Mesh::element_iterator( *( (::Mesh::element_iterator *) begin ) );
            ::Mesh::element_iterator end_element =
                ::Mesh::element_iterator( *( (::Mesh::element_iterator *) end ) );
            d_size = 0;
            while ( cur_element != end_element ) {
                d_size++;
                ++cur_element;
            }
        }
    } else {
        AMP_ERROR( "libMesh does not support iterators over this (unknown) type" );
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
libMeshIterator::libMeshIterator( const libMeshIterator &rhs )
    : MeshIterator() // Note: we never want to call the base copy constructor
{
    d_typeID   = libMeshIteratorTypeID;
    d_iterator = nullptr;
    d_type     = rhs.d_type;
    d_mesh     = rhs.d_mesh;
    d_gcw      = rhs.d_gcw;
    d_pos      = rhs.d_pos;
    d_size     = rhs.d_size;
    d_rank     = rhs.d_rank;
    d_meshID   = rhs.d_meshID;
    d_dim      = rhs.d_dim;
    d_element  = &d_cur_element;
    if ( d_type == 0 ) {
        // Node iterator
        d_begin2 = (void *) new ::Mesh::node_iterator( *( (::Mesh::node_iterator *) rhs.d_begin2 ) );
        d_end2   = (void *) new ::Mesh::node_iterator( *( (::Mesh::node_iterator *) rhs.d_end2 ) );
        d_pos2   = (void *) new ::Mesh::node_iterator( *( (::Mesh::node_iterator *) rhs.d_pos2 ) );
    } else if ( d_type == 1 ) {
        // Element iterator
        d_begin2 =
            (void *) new ::Mesh::element_iterator( *( (::Mesh::element_iterator *) rhs.d_begin2 ) );
        d_end2 =
            (void *) new ::Mesh::element_iterator( *( (::Mesh::element_iterator *) rhs.d_end2 ) );
        d_pos2 =
            (void *) new ::Mesh::element_iterator( *( (::Mesh::element_iterator *) rhs.d_pos2 ) );
    } else {
        AMP_ERROR( "libMesh does not support iterators over this (unknown) type" );
    }
    setCurrentElement();
}
libMeshIterator &libMeshIterator::operator=( const libMeshIterator &rhs )
{
    if ( this == &rhs ) // protect against invalid self-assignment
        return *this;
    this->d_typeID   = libMeshIteratorTypeID;
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
            (void *) new ::Mesh::node_iterator( *( (::Mesh::node_iterator *) rhs.d_begin2 ) );
        this->d_end2 =
            (void *) new ::Mesh::node_iterator( *( (::Mesh::node_iterator *) rhs.d_end2 ) );
        this->d_pos2 =
            (void *) new ::Mesh::node_iterator( *( (::Mesh::node_iterator *) rhs.d_pos2 ) );
    } else if ( this->d_type == 1 ) {
        // Element iterator
        this->d_begin2 =
            (void *) new ::Mesh::element_iterator( *( (::Mesh::element_iterator *) rhs.d_begin2 ) );
        this->d_end2 =
            (void *) new ::Mesh::element_iterator( *( (::Mesh::element_iterator *) rhs.d_end2 ) );
        this->d_pos2 =
            (void *) new ::Mesh::element_iterator( *( (::Mesh::element_iterator *) rhs.d_pos2 ) );
    } else {
        AMP_ERROR( "libMesh does not support iterators over this (unknown) type" );
    }
    setCurrentElement();
    return *this;
}


/********************************************************
* Function to clone the iterator                        *
********************************************************/
MeshIterator *libMeshIterator::clone() const { return new libMeshIterator( *this ); }


/********************************************************
* De-constructor                                        *
********************************************************/
libMeshIterator::~libMeshIterator()
{
    if ( d_pos2 != nullptr ) {
        if ( d_type == 0 ) {
            // Node iterator
            delete (::Mesh::node_iterator *) d_pos2;
            delete (::Mesh::node_iterator *) d_begin2;
            delete (::Mesh::node_iterator *) d_end2;
        } else if ( d_type == 1 ) {
            // Element iterator
            delete (::Mesh::element_iterator *) d_pos2;
            delete (::Mesh::element_iterator *) d_begin2;
            delete (::Mesh::element_iterator *) d_end2;
        }
    }
    d_iterator = nullptr;
}


/********************************************************
* Return an iterator to the beginning or end            *
********************************************************/
MeshIterator libMeshIterator::begin() const
{
    return libMeshIterator( d_type, d_mesh, d_gcw, d_begin2, d_end2, d_begin2, d_size, 0 );
}
MeshIterator libMeshIterator::end() const
{
    return libMeshIterator( d_type, d_mesh, d_gcw, d_begin2, d_end2, d_end2, d_size, d_size );
}


/********************************************************
* Increment/Decrement the iterator                      *
********************************************************/
MeshIterator &libMeshIterator::operator++()
{
    // Prefix increment (increment and return this)
    d_pos++;
    if ( d_type == 0 ) {
        // Node iterator
        ::Mesh::node_iterator *it = (::Mesh::node_iterator *) d_pos2;
        it->operator++();
    } else if ( d_type == 1 ) {
        // Element iterator
        ::Mesh::element_iterator *it = (::Mesh::element_iterator *) d_pos2;
        it->operator++();
    } else {
        AMP_ERROR( "libMesh does not support iterators over this (unknown) type" );
    }
    setCurrentElement();
    return *this;
}
MeshIterator libMeshIterator::operator++( int )
{
    // Postfix increment (increment and return temporary object)
    libMeshIterator tmp( *this ); // Create a temporary variable
    this->operator++();           // apply operator
    return tmp;                   // return temporary result
}
MeshIterator &libMeshIterator::operator--()
{
    // Prefix decrement (increment and return this)
    AMP_ERROR( "Decrementing libMesh iterators is not supported" );
    return *this;
}
MeshIterator libMeshIterator::operator--( int )
{
    // Postfix decrement (increment and return temporary object)
    libMeshIterator tmp( *this ); // Create a temporary variable
    --( *this );                  // apply operator
    return tmp;                   // return temporary result
}


/********************************************************
* Random access incrementors                            *
********************************************************/
MeshIterator libMeshIterator::operator+( int n ) const
{
    libMeshIterator tmp( *this ); // Create a temporary iterator
    tmp.operator+=( n );          // Increment temporary iterator
    return tmp;
}
MeshIterator &libMeshIterator::operator+=( int n )
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
        ::Mesh::node_iterator *it = (::Mesh::node_iterator *) d_pos2;
        if ( n >= 0 ) {
            for ( int i = 0; i < n; i++ )
                it->operator++();
        } else {
            AMP_ERROR( "Decrementing libMesh iterators is not supported" );
        }
    } else if ( d_type == 1 ) {
        // Element iterator
        ::Mesh::element_iterator *it = (::Mesh::element_iterator *) d_pos2;
        if ( n >= 0 ) {
            for ( int i = 0; i < n; i++ )
                it->operator++();
        } else {
            AMP_ERROR( "Decrementing libMesh iterators is not supported" );
        }
    } else {
        AMP_ERROR( "libMesh does not support iterators over this (unknown) type" );
    }
    d_pos += n;
    setCurrentElement();
    return *this;
}


/********************************************************
* Compare two iterators                                 *
********************************************************/
bool libMeshIterator::operator==( const MeshIterator &rhs ) const
{
    libMeshIterator *rhs2 = nullptr;
    libMeshIterator *tmp =
        (libMeshIterator
             *) &rhs; // Convert rhs to a libMeshIterator* so we can access the base class members
    if ( typeid( rhs ) == typeid( libMeshIterator ) ) {
        rhs2 = tmp; // We can safely cast rhs to a libMeshIterator
    } else if ( tmp->d_typeID == libMeshIteratorTypeID ) {
        rhs2 = tmp; // We can safely cast rhs.iterator to a libMeshIterator
    } else if ( ( (libMeshIterator *) tmp->d_iterator )->d_typeID == libMeshIteratorTypeID ) {
        rhs2 = (libMeshIterator *) tmp->d_iterator;
    }
    // Perform direct comparisions if we are dealing with two libMeshIterators
    if ( rhs2 != nullptr ) {
        if ( d_type != rhs2->d_type ) {
            // We are dealing with different types of elements
            return false;
        }
        if ( d_type == 0 ) {
            // Node iterator
            return ( *( (::Mesh::node_iterator *) d_pos2 ) ) ==
                   ( *( (::Mesh::node_iterator *) rhs2->d_pos2 ) );
        } else if ( d_type == 1 ) {
            // Element iterator
            return ( *( (::Mesh::element_iterator *) d_pos2 ) ) ==
                   ( *( (::Mesh::element_iterator *) rhs2->d_pos2 ) );
        } else {
            AMP_ERROR( "libMesh does not support iterators over this (unknown) type" );
        }
    }
    /* We are comparing a libMeshIterator to an arbitrary iterator
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
bool libMeshIterator::operator!=( const MeshIterator &rhs ) const { return !( ( *this ) == rhs ); }


/********************************************************
* Dereference the iterator to get the element           *
********************************************************/
void libMeshIterator::setCurrentElement()
{
    if ( d_pos >= d_size ) {
        d_cur_element = libMeshElement();
    } else if ( d_type == 0 ) {
        // Node iterator
        ::Mesh::node_iterator *it = (::Mesh::node_iterator *) d_pos2;
        ::Node *node              = it->operator*();
        d_cur_element = libMeshElement( d_dim, Vertex, (void *) node, d_rank, d_meshID, d_mesh );
    } else if ( d_type == 1 ) {
        // Element iterator
        ::Mesh::element_iterator *it = (::Mesh::element_iterator *) d_pos2;
        ::Elem *elem                 = it->operator*();
        d_cur_element =
            libMeshElement( d_dim, (GeomType) d_dim, (void *) elem, d_rank, d_meshID, d_mesh );
    } else {
        AMP_ERROR( "libMesh does not support iterators over this (unknown) type" );
    }
}


} // Mesh namespace
} // AMP namespace
