
#ifndef included_AMP_DTK_AMPMeshEntityIterator
#define included_AMP_DTK_AMPMeshEntityIterator

#include <functional>

#include "utils/AMP_MPI.h"

#include "ampmesh/MeshIterator.h"

#include "DTKAMPMeshEntity.h"
#include <DTK_EntityIterator.hpp>

namespace AMP {
namespace Operator {


/**
  * AMP Mesh element implementation for DTK EntityIterator interface.
*/
class AMPMeshEntityIterator : public DataTransferKit::EntityIterator
{
public:
    /*!
     * \brief Default constructor.
     */
    AMPMeshEntityIterator();

    /**
     * Constructor.
     */
    explicit AMPMeshEntityIterator(
        const AMP::Mesh::MeshIterator &iterator,
        const std::function<bool( DataTransferKit::Entity )> &predicate );

    /*!
     * \brief Copy constructor.
     */
    AMPMeshEntityIterator( const AMPMeshEntityIterator &rhs );

    /*!
     * \brief Assignment operator.
     */
    AMPMeshEntityIterator &operator=( const AMPMeshEntityIterator &rhs );

    /*!
     * \brief Destructor.
     */
    ~AMPMeshEntityIterator();

    // Pre-increment operator.
    DataTransferKit::EntityIterator &operator++();

    // Dereference operator.
    DataTransferKit::Entity &operator*( void );

    // Dereference operator.
    DataTransferKit::Entity *operator->( void );

    // Equal comparison operator.
    bool operator==( const DataTransferKit::EntityIterator &rhs ) const;

    // Not equal comparison operator.
    bool operator!=( const DataTransferKit::EntityIterator &rhs ) const;

    // An iterator assigned to the first valid element in the iterator.
    DataTransferKit::EntityIterator begin() const;

    // An iterator assigned to the end of all elements under the iterator.
    DataTransferKit::EntityIterator end() const;

    // Create a clone of the iterator. We need this for the copy constructor
    // and assignment operator to pass along the underlying implementation.
    DataTransferKit::EntityIterator *clone() const;

private:
    // AMP iterator.
    AMP::Mesh::MeshIterator d_amp_iterator;

    // Current AMP entity.
    DataTransferKit::Entity d_current_entity;
};


} // end namespace Operator
} // end namespace AMP

#endif
