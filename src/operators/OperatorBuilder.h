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
 * \details  This function will create a new operator given mesh and input database,
 * \param mesh                  Mesh for the operator
 * \param operatorName          Name of the operator to create
 * \param input_db              Input database
 */
std::shared_ptr<Operator> createOperator( std::shared_ptr<AMP::Mesh::Mesh> mesh,
                                          const std::string &operatorName,
                                          std::shared_ptr<AMP::Database> input_db );

/**
 * \brief Create operator from database
 * \details  This function will create a new operator given mesh, input database, and
 *      elementPhysicsModel
 * \param mesh                  Mesh for the operator
 * \param operatorName          Name of the operator to create
 * \param input_db              Input database
 * \param elementPhysicsModel   Element physics model to use
 */
[[deprecated]] std::shared_ptr<Operator>
createOperator( std::shared_ptr<AMP::Mesh::Mesh> mesh,
                const std::string &operatorName,
                std::shared_ptr<AMP::Database> input_db,
                std::shared_ptr<AMP::Operator::ElementPhysicsModel> elementPhysicsModel );

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
 * \details  This function will create a new operator given mesh, and input database
 * \param mesh                  Mesh for the operator
 * \param boundaryOperatorName  Name of the operator to create
 * \param input_db              Input database
 * \param volumeOperator        GeomType::Cell operator to use
 */
std::shared_ptr<BoundaryOperator>
createBoundaryOperator( std::shared_ptr<AMP::Mesh::Mesh> mesh,
                        std::string boundaryOperatorName,
                        std::shared_ptr<AMP::Database> input_db,
                        AMP::Operator::Operator::shared_ptr volumeOperator );


} // namespace AMP::Operator::OperatorBuilder


#endif
