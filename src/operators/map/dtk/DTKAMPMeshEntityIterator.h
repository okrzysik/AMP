
#ifndef included_AMP_DTK_AMPMeshEntityIterator
#define included_AMP_DTK_AMPMeshEntityIterator

#include <functional>
#include <unordered_map>

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
        const AMP::shared_ptr<std::unordered_map<int, int>> &rank_map,
        const AMP::shared_ptr<std::map<AMP::Mesh::MeshElementID, DataTransferKit::EntityId>>
            &id_map,
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
    DataTransferKit::EntityIterator &operator++() override;

    // Dereference operator.
    DataTransferKit::Entity &operator*(void) override;

    // Dereference operator.
    DataTransferKit::Entity *operator->(void) override;

    // Equal comparison operator.
    bool operator==( const DataTransferKit::EntityIterator &rhs ) const override;

    // Not equal comparison operator.
    bool operator!=( const DataTransferKit::EntityIterator &rhs ) const override;

    // An iterator assigned to the first valid element in the iterator.
    DataTransferKit::EntityIterator begin() const override;

    // An iterator assigned to the end of all elements under the iterator.
    DataTransferKit::EntityIterator end() const override;

    // Create a clone of the iterator. We need this for the copy constructor
    // and assignment operator to pass along the underlying implementation.
    std::unique_ptr<DataTransferKit::EntityIterator> clone() const override;

private:
    // AMP iterator.
    AMP::Mesh::MeshIterator d_amp_iterator;

    // Current AMP entity.
    DataTransferKit::Entity d_current_entity;

    // Global rank map.
    AMP::shared_ptr<std::unordered_map<int, int>> d_rank_map;

    // Mesh id map.
    AMP::shared_ptr<std::map<AMP::Mesh::MeshElementID, DataTransferKit::EntityId>> d_id_map;
};


} // end namespace Operator
} // end namespace AMP

#endif
