
#include "DTKAMPMeshEntityIterator.h"

namespace AMP {
namespace Operator {


//---------------------------------------------------------------------------//
/*!
 * \brief Default constructor.
 */
AMPMeshEntityIterator::AMPMeshEntityIterator() { this->b_iterator_impl = NULL; }

//---------------------------------------------------------------------------//
/**
 * Constructor.
 */
AMPMeshEntityIterator::AMPMeshEntityIterator(
    const AMP::shared_ptr<std::unordered_map<int, int>> &rank_map,
    const AMP::shared_ptr<std::map<AMP::Mesh::MeshElementID, DataTransferKit::EntityId>> &id_map,
    const AMP::Mesh::MeshIterator &iterator,
    const std::function<bool( DataTransferKit::Entity )> &predicate )
    : d_amp_iterator( iterator.begin() ), d_rank_map( rank_map ), d_id_map( id_map )
{
    this->b_iterator_impl = NULL;
    this->b_predicate     = predicate;
}

//---------------------------------------------------------------------------//
/*!
 * \brief Copy constructor.
 */
AMPMeshEntityIterator::AMPMeshEntityIterator( const AMPMeshEntityIterator &rhs )
    : d_amp_iterator( rhs.d_amp_iterator )
{
    this->d_rank_map      = rhs.d_rank_map;
    this->d_id_map        = rhs.d_id_map;
    this->b_iterator_impl = NULL;
    this->b_predicate     = rhs.b_predicate;
}

//---------------------------------------------------------------------------//
/*!
 * \brief Assignment operator.
 */
AMPMeshEntityIterator &AMPMeshEntityIterator::operator=( const AMPMeshEntityIterator &rhs )
{
    this->d_rank_map      = rhs.d_rank_map;
    this->d_id_map        = rhs.d_id_map;
    this->b_iterator_impl = NULL;
    this->b_predicate     = rhs.b_predicate;
    if ( &rhs == this ) {
        return *this;
    }
    d_amp_iterator = rhs.d_amp_iterator;
    return *this;
}

//---------------------------------------------------------------------------//
// Destructor.
AMPMeshEntityIterator::~AMPMeshEntityIterator() { this->b_iterator_impl = NULL; }

//---------------------------------------------------------------------------//
// Pre-increment operator.
DataTransferKit::EntityIterator &AMPMeshEntityIterator::operator++()
{
    ++d_amp_iterator;
    return *this;
}

//---------------------------------------------------------------------------//
// Dereference operator.
DataTransferKit::Entity &AMPMeshEntityIterator::operator*( void )
{
    this->operator->();
    return d_current_entity;
}

//---------------------------------------------------------------------------//
// Dereference operator.
DataTransferKit::Entity *AMPMeshEntityIterator::operator->( void )
{
    d_current_entity = AMPMeshEntity( *d_amp_iterator, *d_rank_map, *d_id_map );
    return &d_current_entity;
}

//---------------------------------------------------------------------------//
// Equal comparison operator.
bool AMPMeshEntityIterator::operator==( const DataTransferKit::EntityIterator &rhs ) const
{
    const AMPMeshEntityIterator *rhs_it = static_cast<const AMPMeshEntityIterator *>( &rhs );
    const AMPMeshEntityIterator *rhs_it_impl =
        static_cast<const AMPMeshEntityIterator *>( rhs_it->b_iterator_impl.get() );
    return ( rhs_it_impl->d_amp_iterator == d_amp_iterator );
}

//---------------------------------------------------------------------------//
// Not equal comparison operator.
bool AMPMeshEntityIterator::operator!=( const DataTransferKit::EntityIterator &rhs ) const
{
    const AMPMeshEntityIterator *rhs_it = static_cast<const AMPMeshEntityIterator *>( &rhs );
    const AMPMeshEntityIterator *rhs_it_impl =
        static_cast<const AMPMeshEntityIterator *>( rhs_it->b_iterator_impl.get() );
    return ( rhs_it_impl->d_amp_iterator != d_amp_iterator );
}

//---------------------------------------------------------------------------//
// An iterator assigned to the first valid element in the iterator.
DataTransferKit::EntityIterator AMPMeshEntityIterator::begin() const
{
    return AMPMeshEntityIterator( d_rank_map, d_id_map, d_amp_iterator, this->b_predicate );
}

//---------------------------------------------------------------------------//
// An iterator assigned to the end of all elements under the iterator.
DataTransferKit::EntityIterator AMPMeshEntityIterator::end() const
{
    AMPMeshEntityIterator end_it( d_rank_map, d_id_map, d_amp_iterator, this->b_predicate );
    end_it.d_amp_iterator = d_amp_iterator.end();
    return end_it;
}

//---------------------------------------------------------------------------//
// Create a clone of the iterator. We need this for the copy constructor
// and assignment operator to pass along the underlying implementation.
std::unique_ptr<DataTransferKit::EntityIterator> AMPMeshEntityIterator::clone() const
{
    return std::unique_ptr<DataTransferKit::EntityIterator>( new AMPMeshEntityIterator( *this ) );
}

//---------------------------------------------------------------------------//
} // namespace Operator
} // namespace AMP
