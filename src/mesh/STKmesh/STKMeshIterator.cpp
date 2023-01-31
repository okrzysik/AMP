#include "AMP/mesh/STKmesh/STKMeshIterator.h"
#include "AMP/mesh/STKmesh/STKMeshElement.h"

#include "stk_mesh/base/BulkData.hpp"

namespace AMP::Mesh {


// unused global variable to prevent compiler warning
static MeshElement nullElement;

/********************************************************
 * Constructors                                          *
 ********************************************************/
STKMeshIterator::STKMeshIterator()
    : MeshIterator(),
      d_gcw( 0 ),
      d_dim( 0 ),
      d_rank( 0 ),
      d_meshID(),
      d_mesh( 0 ),
      d_entries(),
      d_pos(),
      d_cur_element()
{
    iterator = NULL;
    typeID   = getTypeID<decltype( *this )>();
}

STKMeshIterator::STKMeshIterator( const AMP::Mesh::STKMesh *mesh,
                                  int gcw,
                                  std::vector<stk::mesh::Entity *> &entries )
    : MeshIterator(),
      d_gcw( gcw ),
      d_dim( mesh->getSTKMeshMeta()->spatial_dimension() ),
      d_rank( mesh->getComm().getRank() ),
      d_meshID( mesh->meshID() ),
      d_mesh( mesh ),
      d_entries( new std::vector<stk::mesh::Entity *>( entries ) ),
      d_pos( d_entries->begin() ),
      d_cur_element()
{
    iterator = NULL;
    typeID   = getTypeID<decltype( *this )>();
}
STKMeshIterator::STKMeshIterator( const AMP::Mesh::STKMesh *mesh,
                                  int gcw,
                                  std::shared_ptr<std::vector<stk::mesh::Entity *>> entries )
    : MeshIterator(),
      d_gcw( gcw ),
      d_dim( mesh->getSTKMeshMeta()->spatial_dimension() ),
      d_rank( mesh->getComm().getRank() ),
      d_meshID( mesh->meshID() ),
      d_mesh( mesh ),
      d_entries( entries ),
      d_pos( d_entries->begin() ),
      d_cur_element()
{
    iterator = NULL;
    typeID   = getTypeID<decltype( *this )>();
}
STKMeshIterator::STKMeshIterator( const STKMeshIterator &rhs )
    : MeshIterator(),
      d_gcw( rhs.d_gcw ),
      d_dim( rhs.d_dim ),
      d_rank( rhs.d_rank ),
      d_meshID( rhs.d_meshID ),
      d_mesh( rhs.d_mesh ),
      d_entries( rhs.d_entries ),
      d_pos( d_entries->begin() + rhs.position() ),
      d_cur_element( rhs.d_cur_element )
{
    iterator = NULL;
    typeID   = getTypeID<decltype( *this )>();
}
STKMeshIterator &STKMeshIterator::operator=( const STKMeshIterator &rhs )
{
    this->iterator = NULL;
    if ( this == &rhs ) // protect against invalid self-assignment
        return *this;
    this->typeID        = getTypeID<decltype( *this )>();
    this->d_gcw         = rhs.d_gcw;
    this->d_dim         = rhs.d_dim;
    this->d_rank        = rhs.d_rank;
    this->d_meshID      = rhs.d_meshID;
    this->d_mesh        = rhs.d_mesh;
    this->d_entries     = rhs.d_entries;
    this->d_pos         = this->d_entries->begin() + rhs.position();
    this->d_cur_element = rhs.d_cur_element;
    return *this;
}


/********************************************************
 * Function to clone the iterator                        *
 ********************************************************/
MeshIterator *STKMeshIterator::clone() const { return new STKMeshIterator( *this ); }


/********************************************************
 * De-constructor                                        *
 ********************************************************/
STKMeshIterator::~STKMeshIterator() {}


/********************************************************
 * Return an iterator to the beginning or end            *
 ********************************************************/
MeshIterator STKMeshIterator::begin() const
{
    STKMeshIterator e( d_mesh, d_gcw, d_entries );
    return e;
}
MeshIterator STKMeshIterator::end() const
{
    STKMeshIterator e( d_mesh, d_gcw, d_entries );
    e.d_pos = d_entries->end();
    return e;
}


/********************************************************
 * Return the number of elements in the iterator         *
 ********************************************************/
size_t STKMeshIterator::size() const { return d_entries->size(); }
size_t STKMeshIterator::position() const
{
    const int size = d_pos - d_entries->begin();
    return size;
}


/********************************************************
 * Increment/Decrement the iterator                      *
 ********************************************************/
MeshIterator &STKMeshIterator::operator++()
{
    // Prefix increment (increment and return this)
    ++d_pos;
    return *this;
}
MeshIterator STKMeshIterator::operator++( int )
{
    // Postfix increment (increment and return temporary object)
    STKMeshIterator tmp( *this ); // Create a temporary variable
    this->operator++();           // apply operator
    return tmp;                   // return temporary result
}
MeshIterator &STKMeshIterator::operator--()
{
    // Prefix decrement (increment and return this)
    --d_pos;
    return *this;
}
MeshIterator STKMeshIterator::operator--( int )
{
    // Postfix decrement (increment and return temporary object)
    STKMeshIterator tmp( *this ); // Create a temporary variable
    --( *this );                  // apply operator
    return tmp;                   // return temporary result
}


/********************************************************
 * Compare two iterators                                 *
 ********************************************************/
bool STKMeshIterator::operator==( const MeshIterator &rhs ) const
{
    const STKMeshIterator *rhs2 = nullptr;
    // Convert rhs to a STKMeshIterator* so we can access the base class members
    const auto *tmp = reinterpret_cast<const STKMeshIterator *>( &rhs );
    if ( tmp->typeID == getTypeID<decltype( *this )>() ) {
        rhs2 = tmp // We can safely cast rhs to a STKMeshIterator
    } else if ( tmp->d_iterator != nullptr ) {
        tmp = reinterpret_cast<const STKMeshIterator *>( tmp->d_iterator );
        if ( tmp->d_typeID == getTypeID<decltype( *this )>() )
            rhs2 = tmp; // We can safely cast rhs.iterator to a STKMeshIterator
    }
    // Perform direct comparisions if we are dealing with two STKMeshIterators
    if ( rhs2 != NULL )
        return ( d_pos == rhs2->d_pos );

    /* We are comparing a STKMeshIterator to an arbitrary iterator
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
    MeshIterator iterator1 = this->begin();
    MeshIterator iterator2 = rhs.begin();
    bool elements_match    = true;
    for ( size_t i = 0; i < this->size(); i++ ) {
        if ( iterator1->globalID() != iterator2->globalID() )
            elements_match = false;
        ++iterator1;
        ++iterator2;
    }
    return elements_match;
}
bool STKMeshIterator::operator!=( const MeshIterator &rhs ) const { return !( ( *this ) == rhs ); }


/********************************************************
 * Dereference the iterator to get the element           *
 ********************************************************/
MeshElement &STKMeshIterator::operator*()
{
    this->operator->(); // Initialize d_cur_element
    return d_cur_element;
}
MeshElement *STKMeshIterator::operator->()
{
    d_cur_element = STKMeshElement( d_dim, *d_pos, d_rank, d_meshID, d_mesh );
    return &d_cur_element;
}


} // namespace AMP::Mesh
