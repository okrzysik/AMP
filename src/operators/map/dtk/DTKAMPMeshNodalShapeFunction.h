
#ifndef included_AMP_DTK_AMPMeshNodalShapeFunction
#define included_AMP_DTK_AMPMeshNodalShapeFunction

#include "utils/AMP_MPI.h"

#include "discretization/DOF_Manager.h"

#include <DTK_EntityShapeFunction.hpp>

namespace AMP {
namespace Operator {


/**
  * AMP Mesh nodal shape function implementation for DTK EntityShapeFunction
  * interface. For now this API is only for Hex-8 elements and nodes.
*/
class AMPMeshNodalShapeFunction : public DataTransferKit::EntityShapeFunction
{
public:
    /**
     * Constructor.
     */
    explicit AMPMeshNodalShapeFunction(
        const AMP::shared_ptr<AMP::Discretization::DOFManager> &dof_manager );

    //! Destructor
    ~AMPMeshNodalShapeFunction() {}

    /*!
     * \brief Given an entity, get the ids of the degrees of freedom in the
     * vector space supporting its shape function.
     * \param entity Get the degrees of freedom for this entity.
     * \param dof_ids Return the ids of the degrees of freedom in the parallel
     * vector space supporting the entities.
     */
    void entityDOFIds( const DataTransferKit::Entity &entity,
                       Teuchos::Array<std::size_t> &dof_ids ) const;

    /*!
     * \brief Given an entity and a reference point, evaluate the shape
     * function of the entity at that point.
     * \param entity Evaluate the shape function of this entity.
     * \param reference_point Evaluate the shape function at this point
     * given in reference coordinates.
     * \param values Entity shape function evaluated at the reference
     * point.
     */
    void evaluateValue( const DataTransferKit::Entity &entity,
                        const Teuchos::ArrayView<const double> &reference_point,
                        Teuchos::Array<double> &values ) const;

    /*!
     * \brief Given an entity and a reference point, evaluate the gradient of
     * the shape function of the entity at that point.
     * \param entity Evaluate the shape function of this entity.
     * \param reference_point Evaluate the shape function at this point
     * given in reference coordinates.
     * \param gradients Entity shape function gradients evaluated at the reference
     * point. Return these ordered with respect to those return by
     * getDOFIds() such that gradients[N][D] gives the gradient value of the
     * Nth DOF in the Dth spatial dimension.
     */
    void evaluateGradient( const DataTransferKit::Entity &entity,
                           const Teuchos::ArrayView<const double> &reference_point,
                           Teuchos::Array<Teuchos::Array<double>> &gradients ) const;

private:
    // DOF manager.
    AMP::shared_ptr<AMP::Discretization::DOFManager> d_dof_manager;
};
}
}

#endif
