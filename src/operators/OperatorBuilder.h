#ifndef included_AMP_OperatorBuilder
#define included_AMP_OperatorBuilder

#include "AMP/mesh/Mesh.h"
#include "AMP/operators/Operator.h"
#include "AMP/operators/boundary/BoundaryOperator.h"
#include "AMP/utils/Database.h"
#include "ElementPhysicsModel.h"
#include "ElementPhysicsModelFactory.h"

namespace AMP {
namespace Operator {


/**
 * \class OperatorBuilder
 * \brief A class used to create operators
 * \details  This class provides routines for creating operators.
 *  This is a helper class that simplifies operator creation in tests.
 */
class OperatorBuilder
{
public:
    /**
     * \brief Create operator from parameters
     * \details  This function will create a new operator given the parameters
     * \param in_params  Parameters for constructing the operator
     */
    static std::shared_ptr<Operator>
    createOperator( std::shared_ptr<OperatorParameters> in_params );

    /**
     * \brief Create operator from database
     * \details  This function will create a new operator given mesh, input database,
     * elementPhysicsModel, and
     * localModelFactory
     * \param mesh                  Mesh for the operator
     * \param operatorName          Name of the operator to create
     * \param input_db              Input database
     * \param elementPhysicsModel   Element physics model to use
     * \param localModelFactory     Local model factor to use
     */
    static std::shared_ptr<Operator>
    createOperator( AMP::Mesh::Mesh::shared_ptr mesh,
                    std::string operatorName,
                    std::shared_ptr<AMP::Database> input_db,
                    std::shared_ptr<AMP::Operator::ElementPhysicsModel> elementPhysicsModel =
                        std::shared_ptr<AMP::Operator::ElementPhysicsModel>(),
                    std::shared_ptr<AMP::Operator::ElementPhysicsModelFactory> localModelFactory =
                        std::shared_ptr<AMP::Operator::ElementPhysicsModelFactory>() );

    /**
     * \brief Create operator from database
     * \details  This function will create a new operator given mesh, input database
     * \param mesh1                 Mesh1 for the operator
     * \param mesh2                 Mesh2 for the operator
     * \param comm                  Comm to use for the operator
     * \param input_db              Input database
     */
    static std::shared_ptr<Operator> createOperator( AMP::Mesh::Mesh::shared_ptr mesh1,
                                                     AMP::Mesh::Mesh::shared_ptr mesh2,
                                                     const AMP::AMP_MPI &comm,
                                                     std::shared_ptr<AMP::Database> input_db );

    /**
     * \brief Create operator from database
     * \details  This function will create a new operator given mesh, input database,
     * elementPhysicsModel, and
     * localModelFactory
     * \param meshAdapter           Mesh for the operator
     * \param boundaryOperatorName  Name of the operator to create
     * \param input_db              Input database
     * \param volumeOperator        GeomType::Volume operator to use
     * \param elementPhysicsModel   Element physics model to use
     * \param localModelFactory     Local model factor to use
     */
    static std::shared_ptr<BoundaryOperator> createColumnBoundaryOperator(
        AMP::Mesh::Mesh::shared_ptr meshAdapter,
        std::string boundaryOperatorName,
        std::shared_ptr<AMP::Database> input_db,
        AMP::Operator::Operator::shared_ptr volumeOperator,
        std::shared_ptr<AMP::Operator::ElementPhysicsModel> elementPhysicsModel =
            std::shared_ptr<AMP::Operator::ElementPhysicsModel>(),
        std::shared_ptr<AMP::Operator::ElementPhysicsModelFactory> localModelFactory =
            std::shared_ptr<AMP::Operator::ElementPhysicsModelFactory>() );

    /**
     * \brief Create operator from database
     * \details  This function will create a new operator given mesh, input database,
     * elementPhysicsModel, and
     * localModelFactory
     * \param meshAdapter           Mesh for the operator
     * \param boundaryOperatorName  Name of the operator to create
     * \param input_db              Input database
     * \param volumeOperator        GeomType::Volume operator to use
     * \param elementPhysicsModel   Element physics model to use
     * \param localModelFactory     Local model factor to use
     */
    static std::shared_ptr<BoundaryOperator> createBoundaryOperator(
        AMP::Mesh::Mesh::shared_ptr meshAdapter,
        std::string boundaryOperatorName,
        std::shared_ptr<AMP::Database> input_db,
        AMP::Operator::Operator::shared_ptr volumeOperator,
        std::shared_ptr<AMP::Operator::ElementPhysicsModel> elementPhysicsModel =
            std::shared_ptr<AMP::Operator::ElementPhysicsModel>(),
        std::shared_ptr<AMP::Operator::ElementPhysicsModelFactory> localModelFactory =
            std::shared_ptr<AMP::Operator::ElementPhysicsModelFactory>() );

protected:
    OperatorBuilder() {}
    ~OperatorBuilder() {}


    static std::shared_ptr<Operator>
    createIdentityOperator( AMP::Mesh::Mesh::shared_ptr meshAdapter,
                            std::shared_ptr<AMP::Database> input_db );

    static std::shared_ptr<Operator>
    createFlowFrapconOperator( AMP::Mesh::Mesh::shared_ptr meshAdapter,
                               std::shared_ptr<AMP::Database> input_db );

    static std::shared_ptr<Operator>
    createFlowFrapconJacobian( AMP::Mesh::Mesh::shared_ptr meshAdapter,
                               std::shared_ptr<AMP::Database> input_db );

    static std::shared_ptr<Operator> createSubchannelTwoEqLinearOperator(
        AMP::Mesh::Mesh::shared_ptr meshAdapter,
        std::shared_ptr<AMP::Database> input_db,
        std::shared_ptr<AMP::Operator::ElementPhysicsModel> elementPhysicsModel );

    static std::shared_ptr<Operator> createSubchannelTwoEqNonlinearOperator(
        AMP::Mesh::Mesh::shared_ptr meshAdapter,
        std::shared_ptr<AMP::Database> input_db,
        std::shared_ptr<AMP::Operator::ElementPhysicsModel> elementPhysicsModel );

    static std::shared_ptr<Operator> createSubchannelFourEqNonlinearOperator(
        AMP::Mesh::Mesh::shared_ptr meshAdapter,
        std::shared_ptr<AMP::Database> input_db,
        std::shared_ptr<AMP::Operator::ElementPhysicsModel> elementPhysicsModel );

    static std::shared_ptr<Operator>
    createNeutronicsRhsOperator( AMP::Mesh::Mesh::shared_ptr meshAdapter,
                                 std::shared_ptr<AMP::Database> input_db );

    static std::shared_ptr<Operator> createVolumeIntegralOperator(
        AMP::Mesh::Mesh::shared_ptr meshAdapter,
        std::shared_ptr<AMP::Database> input_db,
        std::shared_ptr<AMP::Operator::ElementPhysicsModel> elementPhysicsModel );

    static std::shared_ptr<Operator> createMassLinearFEOperator(
        AMP::Mesh::Mesh::shared_ptr meshAdapter,
        std::shared_ptr<AMP::Database> input_db,
        std::shared_ptr<AMP::Operator::ElementPhysicsModel> elementPhysicsModel );

    static std::shared_ptr<Operator> createLinearDiffusionOperator(
        AMP::Mesh::Mesh::shared_ptr meshAdapter,
        std::shared_ptr<AMP::Database> input_db,
        std::shared_ptr<AMP::Operator::ElementPhysicsModel> elementPhysicsModel );

    static std::shared_ptr<Operator> createNonlinearDiffusionOperator(
        AMP::Mesh::Mesh::shared_ptr meshAdapter,
        std::shared_ptr<AMP::Database> input_db,
        std::shared_ptr<AMP::Operator::ElementPhysicsModel> elementPhysicsModel );

    static std::shared_ptr<Operator> createNonlinearFickSoretOperator(
        AMP::Mesh::Mesh::shared_ptr meshAdapter,
        std::string operatorName,
        std::shared_ptr<AMP::Database> input_db,
        std::shared_ptr<AMP::Operator::ElementPhysicsModel> elementPhysicsModel,
        std::shared_ptr<AMP::Operator::ElementPhysicsModelFactory> localModelFactory );

    static std::shared_ptr<Operator> createGapConductanceOperator(
        AMP::Mesh::Mesh::shared_ptr meshAdapter,
        std::shared_ptr<AMP::Database> input_db,
        std::shared_ptr<AMP::Operator::ElementPhysicsModel> elementPhysicsModel );

    static std::shared_ptr<Operator> createLinearMechanicsOperator(
        AMP::Mesh::Mesh::shared_ptr meshAdapter,
        std::shared_ptr<AMP::Database> input_db,
        std::shared_ptr<AMP::Operator::ElementPhysicsModel> elementPhysicsModel );

    static std::shared_ptr<Operator> createNonlinearMechanicsOperator(
        AMP::Mesh::Mesh::shared_ptr meshAdapter,
        std::shared_ptr<AMP::Database> input_db,
        std::shared_ptr<AMP::Operator::ElementPhysicsModel> elementPhysicsModel );

    static std::shared_ptr<Operator> createLinearNavierStokesLSWFOperator(
        AMP::Mesh::Mesh::shared_ptr meshAdapter,
        std::shared_ptr<AMP::Database> input_db,
        std::shared_ptr<AMP::Operator::ElementPhysicsModel> elementPhysicsModel );

    static std::shared_ptr<Operator> createNonlinearNavierStokesLSWFOperator(
        AMP::Mesh::Mesh::shared_ptr meshAdapter,
        std::shared_ptr<AMP::Database> input_db,
        std::shared_ptr<AMP::Operator::ElementPhysicsModel> elementPhysicsModel );

    static std::shared_ptr<Operator> createLinearBVPOperator(
        AMP::Mesh::Mesh::shared_ptr meshAdapter,
        std::string operatorName,
        std::shared_ptr<AMP::Database> input_db,
        std::shared_ptr<AMP::Operator::ElementPhysicsModel> elementPhysicsModel,
        std::shared_ptr<AMP::Operator::ElementPhysicsModelFactory> localModelFactory );

    static std::shared_ptr<Operator> createNonlinearBVPOperator(
        AMP::Mesh::Mesh::shared_ptr meshAdapter,
        std::string operatorName,
        std::shared_ptr<AMP::Database> input_db,
        std::shared_ptr<AMP::Operator::ElementPhysicsModel> elementPhysicsModel,
        std::shared_ptr<AMP::Operator::ElementPhysicsModelFactory> localModelFactory );


    static std::shared_ptr<BoundaryOperator> createDirichletMatrixCorrection(
        AMP::Mesh::Mesh::shared_ptr meshAdapter,
        std::shared_ptr<AMP::Database> input_db,
        AMP::Operator::Operator::shared_ptr volumeOperator,
        std::shared_ptr<AMP::Operator::ElementPhysicsModel> elementPhysicsModel );

    static std::shared_ptr<BoundaryOperator> createMassMatrixCorrection(
        AMP::Mesh::Mesh::shared_ptr meshAdapter,
        std::shared_ptr<AMP::Database> input_db,
        AMP::Operator::Operator::shared_ptr volumeOperator,
        std::shared_ptr<AMP::Operator::ElementPhysicsModel> elementPhysicsModel );

    static std::shared_ptr<BoundaryOperator> createRobinMatrixCorrection(
        AMP::Mesh::Mesh::shared_ptr meshAdapter,
        std::shared_ptr<AMP::Database> input_db,
        AMP::Operator::Operator::shared_ptr volumeOperator,
        std::shared_ptr<AMP::Operator::ElementPhysicsModel> elementPhysicsModel );

    static std::shared_ptr<BoundaryOperator> createRobinVectorCorrection(
        AMP::Mesh::Mesh::shared_ptr meshAdapter,
        std::shared_ptr<AMP::Database> input_db,
        AMP::Operator::Operator::shared_ptr volumeOperator,
        std::shared_ptr<AMP::Operator::ElementPhysicsModel> elementPhysicsModel );

    static std::shared_ptr<BoundaryOperator> createNeumannVectorCorrection(
        AMP::Mesh::Mesh::shared_ptr meshAdapter,
        std::shared_ptr<AMP::Database> input_db,
        AMP::Operator::Operator::shared_ptr volumeOperator,
        std::shared_ptr<AMP::Operator::ElementPhysicsModel> elementPhysicsModel );

    static std::shared_ptr<BoundaryOperator> createDirichletVectorCorrection(
        AMP::Mesh::Mesh::shared_ptr meshAdapter,
        std::shared_ptr<AMP::Database> input_db,
        AMP::Operator::Operator::shared_ptr volumeOperator,
        std::shared_ptr<AMP::Operator::ElementPhysicsModel> elementPhysicsModel );

    static std::shared_ptr<BoundaryOperator>
    createDirichletVectorCorrection( AMP::Mesh::Mesh::shared_ptr meshAdapter,
                                     std::shared_ptr<AMP::Database> input_db,
                                     std::shared_ptr<AMP::Operator::ElementPhysicsModel> );

    static std::shared_ptr<BoundaryOperator>
    createPressureBoundaryOperator( AMP::Mesh::Mesh::shared_ptr meshAdapter,
                                    std::shared_ptr<AMP::Database> input_db,
                                    std::shared_ptr<AMP::Operator::ElementPhysicsModel> );
};
} // namespace Operator
} // namespace AMP

#endif
