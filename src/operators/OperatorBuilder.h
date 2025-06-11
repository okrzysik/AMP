#ifndef included_AMP_OperatorBuilder
#define included_AMP_OperatorBuilder

#include "AMP/mesh/Mesh.h"
#include "AMP/operators/ElementPhysicsModel.h"
#include "AMP/operators/ElementPhysicsModelFactory.h"
#include "AMP/operators/Operator.h"
#include "AMP/operators/boundary/BoundaryOperator.h"
#include "AMP/utils/Database.h"


/**
 * \namespace OperatorBuilder
 * \brief A class used to create operators
 * \details  This class provides routines for creating operators.
 *  This is a helper class that simplifies operator creation in tests.
 */
namespace AMP::Operator::OperatorBuilder {


/**
 * \brief Create operator from database
 * \details  This function will create a new operator given mesh, input database,
 *      elementPhysicsModel, and localModelFactory
 * \param mesh                  Mesh for the operator
 * \param operatorName          Name of the operator to create
 * \param input_db              Input database
 */
std::shared_ptr<Operator> createOperator( std::shared_ptr<AMP::Mesh::Mesh> mesh,
                                          const std::string &operatorName,
                                          std::shared_ptr<AMP::Database> input_db );

/**
 * \brief Create operator from database
 * \details  This function will create a new operator given mesh, input database,
 *      elementPhysicsModel, and localModelFactory
 * \param mesh                  Mesh for the operator
 * \param operatorName          Name of the operator to create
 * \param input_db              Input database
 * \param elementPhysicsModel   Element physics model to use
 * \param localModelFactory     Local model factor to use
 */
std::shared_ptr<Operator> createOperator(
    std::shared_ptr<AMP::Mesh::Mesh> mesh,
    const std::string &operatorName,
    std::shared_ptr<AMP::Database> input_db,
    std::shared_ptr<AMP::Operator::ElementPhysicsModel> &elementPhysicsModel,
    std::shared_ptr<AMP::Operator::ElementPhysicsModelFactory> localModelFactory = nullptr );

/**
 * \brief Create operator from database
 * \details  This function will create a new operator given mesh, input database
 * \param mesh1                 Mesh1 for the operator
 * \param mesh2                 Mesh2 for the operator
 * \param comm                  Comm to use for the operator
 * \param input_db              Input database
 */
std::shared_ptr<Operator> createOperator( std::shared_ptr<AMP::Mesh::Mesh> mesh1,
                                          std::shared_ptr<AMP::Mesh::Mesh> mesh2,
                                          const AMP::AMP_MPI &comm,
                                          std::shared_ptr<AMP::Database> input_db );

/**
 * \brief Create operator from database
 * \details  This function will create a new operator given mesh, input database,
 * elementPhysicsModel, and
 * localModelFactory
 * \param mesh                  Mesh for the operator
 * \param boundaryOperatorName  Name of the operator to create
 * \param input_db              Input database
 * \param volumeOperator        GeomType::Cell operator to use
 * \param elementPhysicsModel   Element physics model to use
 * \param localModelFactory     Local model factor to use
 */
std::shared_ptr<BoundaryOperator> createColumnBoundaryOperator(
    std::shared_ptr<AMP::Mesh::Mesh> mesh,
    std::string boundaryOperatorName,
    std::shared_ptr<AMP::Database> input_db,
    AMP::Operator::Operator::shared_ptr volumeOperator,
    std::shared_ptr<AMP::Operator::ElementPhysicsModel> elementPhysicsModel      = nullptr,
    std::shared_ptr<AMP::Operator::ElementPhysicsModelFactory> localModelFactory = nullptr );

/**
 * \brief Create operator from database
 * \details  This function will create a new operator given mesh, input database,
 * elementPhysicsModel, and
 * localModelFactory
 * \param meshAdapter           Mesh for the operator
 * \param boundaryOperatorName  Name of the operator to create
 * \param input_db              Input database
 * \param volumeOperator        GeomType::Cell operator to use
 * \param elementPhysicsModel   Element physics model to use
 * \param localModelFactory     Local model factor to use
 */
std::shared_ptr<BoundaryOperator> createBoundaryOperator(
    std::shared_ptr<AMP::Mesh::Mesh> meshAdapter,
    std::string boundaryOperatorName,
    std::shared_ptr<AMP::Database> input_db,
    AMP::Operator::Operator::shared_ptr volumeOperator,
    std::shared_ptr<AMP::Operator::ElementPhysicsModel> elementPhysicsModel      = nullptr,
    std::shared_ptr<AMP::Operator::ElementPhysicsModelFactory> localModelFactory = nullptr );


//! Advanced functions

std::shared_ptr<Operator> createSubchannelTwoEqLinearOperator(
    std::shared_ptr<AMP::Mesh::Mesh> meshAdapter,
    std::shared_ptr<AMP::Database> input_db,
    std::shared_ptr<AMP::Operator::ElementPhysicsModel> &elementPhysicsModel );

std::shared_ptr<Operator> createSubchannelTwoEqNonlinearOperator(
    std::shared_ptr<AMP::Mesh::Mesh> meshAdapter,
    std::shared_ptr<AMP::Database> input_db,
    std::shared_ptr<AMP::Operator::ElementPhysicsModel> &elementPhysicsModel );

std::shared_ptr<Operator> createSubchannelFourEqNonlinearOperator(
    std::shared_ptr<AMP::Mesh::Mesh> meshAdapter,
    std::shared_ptr<AMP::Database> input_db,
    std::shared_ptr<AMP::Operator::ElementPhysicsModel> &elementPhysicsModel );

std::shared_ptr<Operator> createNeutronicsRhsOperator( std::shared_ptr<AMP::Mesh::Mesh> meshAdapter,
                                                       std::shared_ptr<AMP::Database> input_db );

std::shared_ptr<Operator> createVolumeIntegralOperator(
    std::shared_ptr<AMP::Mesh::Mesh> meshAdapter,
    std::shared_ptr<AMP::Database> input_db,
    std::shared_ptr<AMP::Operator::ElementPhysicsModel> &elementPhysicsModel );

std::shared_ptr<Operator> createMassLinearFEOperator(
    std::shared_ptr<AMP::Mesh::Mesh> meshAdapter,
    std::shared_ptr<AMP::Database> input_db,
    std::shared_ptr<AMP::Operator::ElementPhysicsModel> &elementPhysicsModel );

std::shared_ptr<Operator> createLinearDiffusionOperator(
    std::shared_ptr<AMP::Mesh::Mesh> meshAdapter,
    std::shared_ptr<AMP::Database> input_db,
    std::shared_ptr<AMP::Operator::ElementPhysicsModel> &elementPhysicsModel );

std::shared_ptr<Operator> createNonlinearDiffusionOperator(
    std::shared_ptr<AMP::Mesh::Mesh> meshAdapter,
    std::shared_ptr<AMP::Database> input_db,
    std::shared_ptr<AMP::Operator::ElementPhysicsModel> &elementPhysicsModel );

std::shared_ptr<Operator> createNonlinearFickSoretOperator(
    std::shared_ptr<AMP::Mesh::Mesh> meshAdapter,
    std::string operatorName,
    std::shared_ptr<AMP::Database> input_db,
    std::shared_ptr<AMP::Operator::ElementPhysicsModelFactory> localModelFactory );

std::shared_ptr<Operator> createGapConductanceOperator(
    std::shared_ptr<AMP::Mesh::Mesh> meshAdapter,
    std::shared_ptr<AMP::Database> input_db,
    std::shared_ptr<AMP::Operator::ElementPhysicsModel> &elementPhysicsModel );

std::shared_ptr<Operator> createLinearMechanicsOperator(
    std::shared_ptr<AMP::Mesh::Mesh> meshAdapter,
    std::shared_ptr<AMP::Database> input_db,
    std::shared_ptr<AMP::Operator::ElementPhysicsModel> &elementPhysicsModel );

std::shared_ptr<Operator> createNonlinearMechanicsOperator(
    std::shared_ptr<AMP::Mesh::Mesh> meshAdapter,
    std::shared_ptr<AMP::Database> input_db,
    std::shared_ptr<AMP::Operator::ElementPhysicsModel> &elementPhysicsModel );

std::shared_ptr<Operator> createLinearNavierStokesLSWFOperator(
    std::shared_ptr<AMP::Mesh::Mesh> meshAdapter,
    std::shared_ptr<AMP::Database> input_db,
    std::shared_ptr<AMP::Operator::ElementPhysicsModel> &elementPhysicsModel );

std::shared_ptr<Operator> createNonlinearNavierStokesLSWFOperator(
    std::shared_ptr<AMP::Mesh::Mesh> meshAdapter,
    std::shared_ptr<AMP::Database> input_db,
    std::shared_ptr<AMP::Operator::ElementPhysicsModel> &elementPhysicsModel );

std::shared_ptr<Operator> createLinearBVPOperator(
    std::shared_ptr<AMP::Mesh::Mesh> meshAdapter,
    std::string operatorName,
    std::shared_ptr<AMP::Database> input_db,
    std::shared_ptr<AMP::Operator::ElementPhysicsModel> &elementPhysicsModel,
    std::shared_ptr<AMP::Operator::ElementPhysicsModelFactory> localModelFactory );

std::shared_ptr<Operator> createNonlinearBVPOperator(
    std::shared_ptr<AMP::Mesh::Mesh> meshAdapter,
    std::string operatorName,
    std::shared_ptr<AMP::Database> input_db,
    std::shared_ptr<AMP::Operator::ElementPhysicsModel> &elementPhysicsModel,
    std::shared_ptr<AMP::Operator::ElementPhysicsModelFactory> localModelFactory );


std::shared_ptr<BoundaryOperator>
createDirichletMatrixCorrection( std::shared_ptr<AMP::Mesh::Mesh> meshAdapter,
                                 std::shared_ptr<AMP::Database> input_db,
                                 AMP::Operator::Operator::shared_ptr volumeOperator );

std::shared_ptr<BoundaryOperator>
createMassMatrixCorrection( std::shared_ptr<AMP::Mesh::Mesh> meshAdapter,
                            std::shared_ptr<AMP::Database> input_db,
                            AMP::Operator::Operator::shared_ptr volumeOperator );

std::shared_ptr<BoundaryOperator> createRobinMatrixCorrection(
    std::shared_ptr<AMP::Mesh::Mesh> meshAdapter,
    std::shared_ptr<AMP::Database> input_db,
    AMP::Operator::Operator::shared_ptr volumeOperator,
    std::shared_ptr<AMP::Operator::ElementPhysicsModel> &elementPhysicsModel );

std::shared_ptr<BoundaryOperator> createRobinVectorCorrection(
    std::shared_ptr<AMP::Mesh::Mesh> meshAdapter,
    std::shared_ptr<AMP::Database> input_db,
    AMP::Operator::Operator::shared_ptr volumeOperator,
    std::shared_ptr<AMP::Operator::ElementPhysicsModel> &elementPhysicsModel );

std::shared_ptr<BoundaryOperator> createNeumannVectorCorrection(
    std::shared_ptr<AMP::Mesh::Mesh> meshAdapter,
    std::shared_ptr<AMP::Database> input_db,
    AMP::Operator::Operator::shared_ptr volumeOperator,
    std::shared_ptr<AMP::Operator::ElementPhysicsModel> &elementPhysicsModel );

std::shared_ptr<BoundaryOperator>
createDirichletVectorCorrection( std::shared_ptr<AMP::Mesh::Mesh> meshAdapter,
                                 std::shared_ptr<AMP::Database> input_db,
                                 AMP::Operator::Operator::shared_ptr volumeOperator );

std::shared_ptr<BoundaryOperator>
createDirichletVectorCorrection( std::shared_ptr<AMP::Mesh::Mesh> meshAdapter,
                                 std::shared_ptr<AMP::Database> input_db );


void setNestedOperatorMemoryLocations( std::shared_ptr<AMP::Database> input_db,
                                       std::string outerOperatorName,
                                       std::vector<std::string> nestedOperatorNames );

std::vector<std::string> getActiveVariables( std::shared_ptr<const AMP::Database> db,
                                             const std::string &key );

} // namespace AMP::Operator::OperatorBuilder


#endif
