#include "math.h"

#include "AMP/discretization/DOF_Manager.h"
#include "SimpleVector.h"

namespace AMP {
namespace LinearAlgebra {


/****************************************************************
 * Constructors                                                  *
 ****************************************************************/
template<typename TYPE, typename OPS, typename DATA>
SimpleVector<TYPE, OPS, DATA>::SimpleVector() : Vector()
{
    d_VectorOps = std::make_shared<OPS>();
}

template<typename TYPE, typename OPS, typename DATA>
Vector::shared_ptr SimpleVector<TYPE, OPS, DATA>::create( size_t localSize, const std::string &var )
{
    return create( localSize, std::make_shared<Variable>( var ) );
}

template<typename TYPE, typename OPS, typename DATA>
Vector::shared_ptr SimpleVector<TYPE, OPS, DATA>::create( size_t localSize,
                                                          Variable::shared_ptr var )
{
    AMP_MPI comm( AMP_COMM_SELF );
    auto DOFs = std::make_shared<AMP::Discretization::DOFManager>( localSize, comm );
    auto ptr  = new SimpleVector<TYPE, OPS, DATA>();
    ptr->setVariable( var );
    ptr->allocateVectorData( DOFs->beginDOF(), DOFs->numLocalDOF(), DOFs->numGlobalDOF() );
    ptr->setComm( comm );
    ptr->setDOFManager( DOFs );
    ptr->setCommunicationList(
        AMP::LinearAlgebra::CommunicationList::createEmpty( localSize, comm ) );
    return Vector::shared_ptr( ptr );
}
template<typename TYPE, typename OPS, typename DATA>
Vector::shared_ptr
SimpleVector<TYPE, OPS, DATA>::create( size_t localSize, Variable::shared_ptr var, AMP_MPI comm )
{
    auto DOFs = std::make_shared<AMP::Discretization::DOFManager>( localSize, comm );
    auto ptr  = new SimpleVector<TYPE, OPS, DATA>();
    ptr->setVariable( var );
    ptr->allocateVectorData( DOFs->beginDOF(), DOFs->numLocalDOF(), DOFs->numGlobalDOF() );
    ptr->setComm( comm );
    ptr->setDOFManager( DOFs );
    ptr->setCommunicationList(
        AMP::LinearAlgebra::CommunicationList::createEmpty( localSize, comm ) );
    return Vector::shared_ptr( ptr );
}
template<typename TYPE, typename OPS, typename DATA>
Vector::shared_ptr
SimpleVector<TYPE, OPS, DATA>::create( Variable::shared_ptr var,
                                       AMP::Discretization::DOFManager::shared_ptr DOFs,
                                       AMP::LinearAlgebra::CommunicationList::shared_ptr commlist )
{
    auto ptr = new SimpleVector<TYPE, OPS, DATA>();
    ptr->setVariable( var );
#if 1
    ptr->allocateVectorData( DOFs->beginDOF(), DOFs->numLocalDOF(), DOFs->numGlobalDOF() );
    ptr->setDOFManager( DOFs );
    auto comm = DOFs->getComm();
    ptr->setComm( comm );
#else
    ptr->allocate( DOFs->beginDOF(), DOFs->numLocalDOF(), DOFs->numGlobalDOF() );
    ptr->d_DOFManager = DOFs;
    ptr->d_comm       = DOFs->getComm();
#endif
    ptr->setCommunicationList( commlist );
    return Vector::shared_ptr( ptr );
}
template<typename TYPE, typename OPS, typename DATA>
void SimpleVector<TYPE, OPS, DATA>::resize( size_t N )
{
    d_DOFManager = std::make_shared<AMP::Discretization::DOFManager>( N, d_comm );
#if 1
    this->allocateVectorData(
        d_DOFManager->beginDOF(), d_DOFManager->numLocalDOF(), d_DOFManager->numGlobalDOF() );
#else
    this->allocate(
        d_DOFManager->beginDOF(), d_DOFManager->numLocalDOF(), d_DOFManager->numGlobalDOF() );
#endif
}


/****************************************************************
 * Vector operations                                             *
 ****************************************************************/
template<typename TYPE, typename OPS, typename DATA>
Vector::shared_ptr
SimpleVector<TYPE, OPS, DATA>::cloneVector( const Variable::shared_ptr name ) const
{
    return create( name, d_DOFManager, d_VectorData->getCommunicationList() );
}
template<typename TYPE, typename OPS, typename DATA>
void SimpleVector<TYPE, OPS, DATA>::swapVectors( Vector &rhs )
{
    auto x = dynamic_cast<SimpleVector *>( &rhs );
    AMP_INSIST( x != nullptr, "rhs is not a SimpleVector" );
    std::swap( d_comm, d_comm );
    d_VectorData->swapData( *( rhs.getVectorData() ) );
}

// temporary functions
template<typename TYPE, typename OPS, typename DATA>
void SimpleVector<TYPE, OPS, DATA>::allocateVectorData( size_t localSize,
                                                        size_t numLocal,
                                                        size_t numGlobal )
{
    d_VectorData = std::make_shared<DATA>( localSize, numLocal, numGlobal );
}

template<typename TYPE, typename OPS, typename DATA>
void SimpleVector<TYPE, OPS, DATA>::setDOFManager(
    std::shared_ptr<AMP::Discretization::DOFManager> dofManager )
{
    d_DOFManager = dofManager;
}

} // namespace LinearAlgebra
} // namespace AMP
