
#ifndef included_AMP_DTK_AMPMeshEntityIterator
#define included_AMP_DTK_AMPMeshEntityIterator

#include "utils/shared_ptr.h"
#include "matrices/Matrix.h"
#include "operators/Operator.h"
#include "operators/OperatorParameters.h"
#include "vectors/Vector.h"

#include <DTK_EntityIterator.hpp>
#include "DTKAMPMeshEntity.h"

namespace AMP {
namespace Operator {


/**
  * AMP Mesh element implementation for DTK EntityIteratorExtraData interface.
*/
class AMPMeshEntityIterator : public DataTransferKit::EntityIterator
{
public :

    /*!
     * \brief Default constructor.
     */
    AMPMeshEntityIterator();

    /**
     * Constructor.
     */
    AMPMeshEntityIterator( const AMP::Mesh::MeshIterator& iterator,
			   const std::function<bool(Entity)>& predicate );

    /*!
     * \brief Copy constructor.
     */
    AMPMeshEntityIterator( const AMPMeshEntityIterator& rhs );

    /*!
     * \brief Assignment operator.
     */
    AMPMeshEntityIterator& operator=( const AMPMeshEntityIterator& rhs );

    /*!
     * \brief Destructor.
     */
    ~AMPMeshEntityIterator();

    // Pre-increment operator.
    DataTransferKit::EntityIterator& operator++();

    // Dereference operator.
    DataTransferKit::Entity& operator*(void);

    // Dereference operator.
    DataTransferKit::Entity* operator->(void);

    // Equal comparison operator.
    bool operator==( const DataTransferKit::EntityIterator& rhs ) const;

    // Not equal comparison operator.
    bool operator!=( const DataTransferKit::EntityIterator& rhs ) const;

    // An iterator assigned to the first valid element in the iterator.
    DataTransferKit::EntityIterator begin() const;

    // An iterator assigned to the end of all elements under the iterator.
    DataTransferKit::EntityIterator end() const;

    // Create a clone of the iterator. We need this for the copy constructor
    // and assignment operator to pass along the underlying implementation.
    DataTransferKit::EntityIterator* clone() const;

  private:

    // Set the current entity.
    void setCurrentEntity();

  private:

    // AMP iterator.
    AMP::Mesh::MeshIterator d_iterator;

    // Current AMP entity.
    AMPMeshEntity d_current_entity;
};


}
}

#endif


