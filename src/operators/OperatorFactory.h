#ifndef included_AMP_OperatorFactory
#define included_AMP_OperatorFactory

/* AMP files */
#include "AMP/ampmesh/MeshManager.h"
#include "AMP/utils/Database.h"

namespace AMP {
namespace Operator {

template<typename OPERATOR>
class OperatorFactory
{
public:
    typedef OPERATOR Operator_t;
    typedef typename Operator_t::Parameters OperatorParameters;
    typedef typename Operator_t::Jacobian Jacobian_t;
    typedef typename Jacobian_t::Parameters JacobianParameters;

    static Operator::shared_ptr getOperator( AMP::shared_ptr<AMP::Database> input_db,
                                             AMP::Mesh::MeshManager::Adapter::shared_ptr mesh =
                                                 AMP::Mesh::MeshManager::Adapter::shared_ptr( 0 ) );

    static Operator::shared_ptr getJacobian( Operator::shared_ptr oper,
                                             const AMP::LinearAlgebra::Vector::shared_ptr &vec,
                                             AMP::Mesh::MeshManager::Adapter::shared_ptr mesh =
                                                 AMP::Mesh::MeshManager::Adapter::shared_ptr( 0 ) );
};


template<typename OPERATOR>
Operator::shared_ptr
OperatorFactory<OPERATOR>::getOperator( AMP::shared_ptr<AMP::Database> input_db,
                                        AMP::Mesh::MeshManager::Adapter::shared_ptr mesh )
{
    AMP::shared_ptr<OperatorParameters> params(
        new OperatorParameters( input_db->getDatabase( Operator_t::DBName() ) ) );
    params->d_MeshAdapter = mesh;
    Operator::shared_ptr retVal( new Operator_t( params ) );
    retVal->setInputVariable( AMP::LinearAlgebra::Variable::shared_ptr(
        new typename Operator_t::InputVariable( "factory input" ) ) );
    retVal->setOutputVariable( AMP::LinearAlgebra::Variable::shared_ptr(
        new typename Operator_t::OutputVariable( "factory output" ) ) );
    return retVal;
}

template<typename OPERATOR>
Operator::shared_ptr
OperatorFactory<OPERATOR>::getJacobian( Operator::shared_ptr oper,
                                        const AMP::LinearAlgebra::Vector::shared_ptr &vec,
                                        AMP::Mesh::MeshManager::Adapter::shared_ptr mesh )
{
    AMP::shared_ptr<JacobianParameters> params =
        AMP::dynamic_pointer_cast<JacobianParameters>( oper->getJacobianParameters( vec ) );
    params->d_MeshAdapter = mesh;
    Operator::shared_ptr retVal( new Jacobian_t( params ) );
    retVal->setInputVariable( AMP::LinearAlgebra::Variable::shared_ptr(
        new typename Jacobian_t::InputVariable( "factory jacobian input" ) ) );
    retVal->setOutputVariable( AMP::LinearAlgebra::Variable::shared_ptr(
        new typename Jacobian_t::OutputVariable( "factory jacobian output" ) ) );
    return retVal;
}
} // namespace Operator
} // namespace AMP

#endif
