
#ifndef included_AMP_DTK_AMPMeshEntitySet
#define included_AMP_DTK_AMPMeshEntitySet

#include "ampmesh/Mesh.h"

#include "utils/AMP_MPI.h"

#include <DTK_EntitySet.hpp>

namespace AMP {
namespace Operator {


/**
  * AMP Mesh element implementation for DTK EntitySet interface.
*/
class AMPMeshEntitySet : public DataTransferKit::EntitySet
{
public:
    /**
     * Constructor.
     */
    explicit AMPMeshEntitySet( const AMP::shared_ptr<AMP::Mesh::Mesh> &mesh );

    //@{
    //! Parallel functions.
    /*!
     * \brief Get the parallel communicator for the entity set.
     * \return A reference-counted pointer to the parallel communicator.
     */
    Teuchos::RCP<const Teuchos::Comm<int>> communicator() const override;
    //@}

    //@{
    //! Geometric data functions.
    /*!
     * \brief Return the largest physical dimension of the entities in the
     * set.
     * \return The physical dimension of the set.
     */
    int physicalDimension() const override;
    //@}

    //@{
    //! Entity access functions.
    /*!
     * \brief Given an EntityId, get the entity.
     * \param entity_id Get the entity with this id.
     * \param entity The entity with the given id.
     */
    void getEntity( const DataTransferKit::EntityId entity_id,
		    const int topological_dimension,
                    DataTransferKit::Entity &entity ) const override;

    /*!
     * \brief Get a iterator of the given entity type that satisfy the given
     * predicate.
     * \param entity_type The type of entity to get a iterator for.
     * \param predicate The selection predicate.
     * \return A iterator of entities of the given type.
     */
    DataTransferKit::EntityIterator entityIterator( 
	const int topological_dimension,
	const DataTransferKit::PredicateFunction& predicate ) const override;

    /*!
     * \brief Given an entity, get the entities of the given type that are
     * adjacent to it.
     */
    void getAdjacentEntities( 
	const DataTransferKit::Entity &entity,
	const int topological_dimension,
	Teuchos::Array<DataTransferKit::Entity> &adjacent_entities ) const override;
    //@}

private:
    // Given a DTK entity type, get an AMP GeomType.
    AMP::Mesh::GeomType
    getGeomTypeFromEntityType( const int topological_dimension ) const;

private:
    // AMP mesh.
    AMP::shared_ptr<AMP::Mesh::Mesh> d_amp_mesh;
};
}
}

#endif
