#ifdef USE_AMP_DISCRETIZATION
#ifndef included_AMP_VectorBuider_hpp
#define included_AMP_VectorBuider_hpp

#include "AMP/discretization/DOF_Manager.h"

#include "math.h"


namespace AMP {
namespace LinearAlgebra {


/****************************************************************
 * Constructors                                                  *
 ****************************************************************/
template<typename TYPE, typename OPS, typename DATA>
Vector::shared_ptr createSimpleVector( size_t localSize, const std::string &name )
{
    AMP_MPI comm( AMP_COMM_SELF );
    auto var = std::make_shared<Variable>( name );
    return createSimpleVector<TYPE, OPS, DATA>( localSize, var, comm );
}

template<typename TYPE, typename OPS, typename DATA>
Vector::shared_ptr createSimpleVector( size_t localSize, Variable::shared_ptr var )
{
    AMP_MPI comm( AMP_COMM_SELF );
    return createSimpleVector<TYPE, OPS, DATA>( localSize, var, comm );
}
template<typename TYPE, typename OPS, typename DATA>
Vector::shared_ptr createSimpleVector( size_t localSize, Variable::shared_ptr var, AMP_MPI comm )
{
    auto DOFs = std::make_shared<AMP::Discretization::DOFManager>( localSize, comm );
    auto ops  = std::make_shared<OPS>();
    auto data =
        std::make_shared<DATA>( DOFs->beginDOF(), DOFs->numLocalDOF(), DOFs->numGlobalDOF() );
    data->setCommunicationList(
        AMP::LinearAlgebra::CommunicationList::createEmpty( localSize, comm ) );
    return std::make_shared<Vector>( data, ops, var, DOFs );
}
template<typename TYPE, typename OPS, typename DATA>
Vector::shared_ptr createSimpleVector( Variable::shared_ptr var,
                                       AMP::Discretization::DOFManager::shared_ptr DOFs,
                                       AMP::LinearAlgebra::CommunicationList::shared_ptr commlist )
{
    auto ops = std::make_shared<OPS>();
    auto data =
        std::make_shared<DATA>( DOFs->beginDOF(), DOFs->numLocalDOF(), DOFs->numGlobalDOF() );
    data->setCommunicationList( commlist );
    return std::make_shared<Vector>( data, ops, var, DOFs );
}


} // namespace LinearAlgebra
} // namespace AMP

#endif
#endif
