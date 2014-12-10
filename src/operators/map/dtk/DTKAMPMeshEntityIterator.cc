
#include "DTKAMPMeshEntityIterator.h"

namespace AMP {
namespace Operator {


//---------------------------------------------------------------------------//
/*!
 * \brief Default constructor.
 */
AMPMeshEntityIterator::AMPMeshEntityIterator()
{ 
    this->b_iterator_impl = NULL;
}

//---------------------------------------------------------------------------//
/**
 * Constructor.
 */
AMPMeshEntityIterator::AMPMeshEntityIterator( 
    const AMP::Mesh::MeshIterator& iterator,
    const std::function<bool(DataTransferKit::Entity)>& predicate )
    : d_amp_iterator( iterator.begin() )
{
    this->b_iterator_impl = NULL;
    this->b_predicate = predicate;
}

//---------------------------------------------------------------------------//
/*!
 * \brief Copy constructor.
 */
AMPMeshEntityIterator::AMPMeshEntityIterator( const AMPMeshEntityIterator& rhs )
    : d_amp_iterator( rhs.d_amp_iterator )
{
    this->b_iterator_impl = NULL;
    this->b_predicate = rhs.b_predicate;
}

//---------------------------------------------------------------------------//
/*!
 * \brief Assignment operator.
 */
AMPMeshEntityIterator& AMPMeshEntityIterator::operator=( 
    const AMPMeshEntityIterator& rhs )
{
    this->b_iterator_impl = NULL;
    this->b_predicate = rhs.b_predicate;
    if ( &rhs == this )
    {
	return *this;
    }
    d_amp_iterator = rhs.d_amp_iterator;
    return *this;
}

//---------------------------------------------------------------------------//
// Destructor.
AMPMeshEntityIterator::~AMPMeshEntityIterator()
{
    this->b_iterator_impl = NULL;
}

//---------------------------------------------------------------------------//
// Pre-increment operator.
DataTransferKit::EntityIterator& AMPMeshEntityIterator::operator++()
{
    ++d_amp_iterator;
    return *this;
}

//---------------------------------------------------------------------------//
// Dereference operator.
DataTransferKit::Entity& AMPMeshEntityIterator::operator*(void)
{
    this->operator->();
    return d_current_entity;
}

//---------------------------------------------------------------------------//
// Dereference operator.
DataTransferKit::Entity* AMPMeshEntityIterator::operator->(void)
{
    d_current_entity = AMPMeshEntity( *d_amp_iterator );
    return &d_current_entity;
}

//---------------------------------------------------------------------------//
// Equal comparison operator.
bool AMPMeshEntityIterator::operator==( 
    const DataTransferKit::EntityIterator& rhs ) const
{
    const AMPMeshEntityIterator* rhs_it = 
	static_cast<const AMPMeshEntityIterator*>(&rhs);
    const AMPMeshEntityIterator* rhs_it_impl = 
	static_cast<const AMPMeshEntityIterator*>(rhs_it->b_iterator_impl);
    return ( rhs_it_impl->d_amp_iterator == d_amp_iterator );
}

//---------------------------------------------------------------------------//
// Not equal comparison operator.
bool AMPMeshEntityIterator::operator!=( 
    const DataTransferKit::EntityIterator& rhs ) const
{
    const AMPMeshEntityIterator* rhs_it = 
    	static_cast<const AMPMeshEntityIterator*>(&rhs);
    const AMPMeshEntityIterator* rhs_it_impl = 
    	static_cast<const AMPMeshEntityIterator*>(rhs_it->b_iterator_impl);
    return ( rhs_it_impl->d_amp_iterator != d_amp_iterator );
}

//---------------------------------------------------------------------------//
// An iterator assigned to the first valid element in the iterator.
DataTransferKit::EntityIterator AMPMeshEntityIterator::begin() const
{
    return AMPMeshEntityIterator( d_amp_iterator, this->b_predicate );
}

//---------------------------------------------------------------------------//
// An iterator assigned to the end of all elements under the iterator.
DataTransferKit::EntityIterator AMPMeshEntityIterator::end() const
{
    AMPMeshEntityIterator end_it( d_amp_iterator, this->b_predicate );
    end_it.d_amp_iterator = d_amp_iterator.end();
    return end_it;
}

//---------------------------------------------------------------------------//
// Create a clone of the iterator. We need this for the copy constructor
// and assignment operator to pass along the underlying implementation.
DataTransferKit::EntityIterator* AMPMeshEntityIterator::clone() const
{
    return new AMPMeshEntityIterator( *this );
}

//---------------------------------------------------------------------------//

}
}


