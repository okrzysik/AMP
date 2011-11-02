#include "ampmesh/libmesh/libMeshIterator.h"
#include "ampmesh/libmesh/libMeshElement.h"
#include "utils/Utilities.h"

// libMesh includes
#include "elem.h"

namespace AMP {
namespace Mesh {


// Create a unique id for this class
static int libMeshIteratorTypeID = 2;

// unused global variable to prevent compiler warning
static MeshElement nullElement;


/********************************************************
* Constructors                                          *
********************************************************/
libMeshIterator::libMeshIterator()
{
    typeID = libMeshIteratorTypeID;
    iterator = this;
    d_begin=NULL;
    d_end=NULL;
    d_pos=NULL;
    d_type = -1;
}
libMeshIterator::libMeshIterator(int type, boost::shared_ptr< ::Mesh> mesh, int gcw, void *begin, void *end, void *pos)
{
    typeID = libMeshIteratorTypeID;
    iterator = this;
    d_type = type;
    d_libMesh = mesh;
    d_gcw = gcw;
    if ( d_type==0 ) {
        // Node iterator
        d_begin = (void*) new ::Mesh::node_iterator( *((::Mesh::node_iterator*)begin) );
        d_end   = (void*) new ::Mesh::node_iterator( *((::Mesh::node_iterator*)end) );
        d_pos   = (void*) new ::Mesh::node_iterator( *((::Mesh::node_iterator*)pos) );
    } else if ( d_type==1 ) {
        // Element iterator
        d_begin = (void*) new ::Mesh::element_iterator( *((::Mesh::element_iterator*)begin) );
        d_end   = (void*) new ::Mesh::element_iterator( *((::Mesh::element_iterator*)end) );
        d_pos   = (void*) new ::Mesh::element_iterator( *((::Mesh::element_iterator*)pos) );
    } else {
        AMP_ERROR("libMesh does not support iterators over this (unknown) type");
    }
}
libMeshIterator::libMeshIterator(const libMeshIterator& rhs)
{
    typeID = libMeshIteratorTypeID;
    iterator = this;
    d_type = rhs.d_type;
    d_libMesh = rhs.d_libMesh;
    d_gcw = rhs.d_gcw;
    if ( d_type==0 ) {
        // Node iterator
        d_begin = (void*) new ::Mesh::node_iterator( *((::Mesh::node_iterator*)rhs.d_begin) );
        d_end   = (void*) new ::Mesh::node_iterator( *((::Mesh::node_iterator*)rhs.d_end) );
        d_pos   = (void*) new ::Mesh::node_iterator( *((::Mesh::node_iterator*)rhs.d_pos) );
    } else if ( d_type==1 ) {
        // Element iterator
        d_begin = (void*) new ::Mesh::element_iterator( *((::Mesh::element_iterator*)rhs.d_begin) );
        d_end   = (void*) new ::Mesh::element_iterator( *((::Mesh::element_iterator*)rhs.d_end) );
        d_pos   = (void*) new ::Mesh::element_iterator( *((::Mesh::element_iterator*)rhs.d_pos) );
    } else {
        AMP_ERROR("libMesh does not support iterators over this (unknown) type");
    }
}
libMeshIterator& libMeshIterator::operator=(const libMeshIterator& rhs)
{
    if (this == &rhs) // protect against invalid self-assignment
        return *this;
    this->typeID = libMeshIteratorTypeID;
    this->iterator = this;
    this->d_type = rhs.d_type;
    this->d_libMesh = rhs.d_libMesh;
    this->d_gcw = rhs.d_gcw;
    if ( this->d_type==0 ) {
        // Node iterator
        this->d_begin = (void*) new ::Mesh::node_iterator( *((::Mesh::node_iterator*)rhs.d_begin) );
        this->d_end   = (void*) new ::Mesh::node_iterator( *((::Mesh::node_iterator*)rhs.d_end) );
        this->d_pos   = (void*) new ::Mesh::node_iterator( *((::Mesh::node_iterator*)rhs.d_pos) );
    } else if ( this->d_type==1 ) {
        // Element iterator
        this->d_begin = (void*) new ::Mesh::element_iterator( *((::Mesh::element_iterator*)rhs.d_begin) );
        this->d_end   = (void*) new ::Mesh::element_iterator( *((::Mesh::element_iterator*)rhs.d_end) );
        this->d_pos   = (void*) new ::Mesh::element_iterator( *((::Mesh::element_iterator*)rhs.d_pos) );
    } else {
        AMP_ERROR("libMesh does not support iterators over this (unknown) type");
    }
    return *this;
}


/********************************************************
* Function to clone the iterator                        *
********************************************************/
MeshIterator* libMeshIterator::clone() const
{
    return new libMeshIterator(*this);
}


/********************************************************
* De-constructor                                        *
********************************************************/
libMeshIterator::~libMeshIterator()
{
    if ( d_pos!=NULL ) {
        if ( d_type==0 ) {
            // Node iterator
            delete (::Mesh::node_iterator*) d_pos;
            delete (::Mesh::node_iterator*) d_begin;
            delete (::Mesh::node_iterator*) d_end;
        } else if ( d_type==1 ) {
            // Element iterator
            delete (::Mesh::element_iterator*) d_pos;
            delete (::Mesh::element_iterator*) d_begin;
            delete (::Mesh::element_iterator*) d_end;
        }
    }
    iterator = NULL;
}


/********************************************************
* Return an iterator to the beginning or end            *
********************************************************/
MeshIterator libMeshIterator::begin()
{
    return libMeshIterator( d_type, d_libMesh, d_gcw, d_begin, d_end, d_begin );
}
MeshIterator libMeshIterator::end()
{
    return libMeshIterator( d_type, d_libMesh, d_gcw, d_begin, d_end, d_end );
}


/********************************************************
* Increment/Decrement the iterator                      *
********************************************************/
MeshIterator& libMeshIterator::operator++()
{
    // Prefix increment (increment and return this)
    if ( d_type==0 ) {
        // Node iterator
        ::Mesh::node_iterator* it = (::Mesh::node_iterator*) d_pos;
        it->operator++();
    } else if ( d_type==1 ) {
        // Element iterator
        ::Mesh::element_iterator* it = (::Mesh::element_iterator*) d_pos;
        it->operator++();
    } else {
        AMP_ERROR("libMesh does not support iterators over this (unknown) type");
    }
    return *this;
}
MeshIterator libMeshIterator::operator++(int)
{
    // Postfix increment (increment and return temporary object)
    libMeshIterator tmp(*this);     // Create a temporary variable
    this->operator++();             // apply operator
    return tmp;                     // return temporary result
}
MeshIterator& libMeshIterator::operator--()
{
    // Prefix decrement (increment and return this)
    AMP_ERROR("Decrementing libMesh iterators is not supported");
    return *this;
}
MeshIterator libMeshIterator::operator--(int)
{
    // Postfix decrement (increment and return temporary object)
    libMeshIterator tmp(*this);      // Create a temporary variable
    --(*this);                      // apply operator
    return tmp;                     // return temporary result
}


/********************************************************
* Compare two iterators                                 *
********************************************************/
bool libMeshIterator::operator==(const MeshIterator& rhs)
{
    libMeshIterator* rhs2 = NULL;
    libMeshIterator* tmp = (libMeshIterator*) &rhs;     // Convert rhs to a libMeshIterator* so we can access the base class members
    if ( typeid(rhs)==typeid(libMeshIterator) ) {
        rhs2 = tmp;     // We can safely cast rhs to a libMeshIterator
    } else if ( tmp->typeID==libMeshIteratorTypeID ) {
        rhs2 = tmp;     // We can safely cast rhs.iterator to a libMeshIterator
    } else if ( ((libMeshIterator*)tmp->iterator)->typeID==libMeshIteratorTypeID ) {
        rhs2 = (libMeshIterator*) tmp->iterator;
    } else {
        AMP_ERROR("Error, comparing a libMesh iterator to an unknown iterator");
    }
    if ( rhs2 != NULL ) {
        if ( d_type==0 ) {
            // Node iterator
            return (*((::Mesh::node_iterator*)d_pos)) == (*((::Mesh::node_iterator*)rhs2->d_pos));
        } else if ( d_type==1 ) {
            // Element iterator
            return (*((::Mesh::element_iterator*)d_pos)) == (*((::Mesh::element_iterator*)rhs2->d_pos));
        } else {
            AMP_ERROR("libMesh does not support iterators over this (unknown) type");
        }
    }
    return false;
}
bool libMeshIterator::operator!=(const MeshIterator& rhs)
{
    return !((*this)==rhs);
}


/********************************************************
* Dereference the iterator to get the element           *
********************************************************/
MeshElement& libMeshIterator::operator*()
{
    int dim = d_libMesh->mesh_dimension();
    if ( d_type==0 ) {
        // Node iterator
        ::Mesh::node_iterator* it = (::Mesh::node_iterator*) d_pos;
        ::Node *node = it->operator*();
        d_cur_element = libMeshElement( dim, d_type, (void*) node );
    } else if ( d_type==1 ) {
        // Element iterator
        ::Mesh::element_iterator* it = (::Mesh::element_iterator*) d_pos;
        ::Elem *elem = it->operator*();
        d_cur_element = libMeshElement( dim, d_type, (void*) elem );
    } else {
        AMP_ERROR("libMesh does not support iterators over this (unknown) type");
    }
    return d_cur_element;
}
MeshElement* libMeshIterator::operator->()
{
    d_cur_element = this->operator*();
    return &d_cur_element;
}


} // Mesh namespace
} // AMP namespace

