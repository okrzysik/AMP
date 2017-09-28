#include "math.h"

#include "SimpleVector.h"
#include "discretization/DOF_Manager.h"

namespace AMP {
namespace LinearAlgebra {


/****************************************************************
 * Constructors                                                  *
 ****************************************************************/
template<typename TYPE, typename OPS, typename DATA>
SimpleVector<TYPE, OPS, DATA>::SimpleVector() : Vector(), DATA()
{
}

template<typename TYPE, typename OPS, typename DATA>
Vector::shared_ptr SimpleVector<TYPE, OPS, DATA>::create( size_t localSize, const std::string &var )
{
    return create( localSize, AMP::make_shared<Variable>( var ) );
}

template<typename TYPE, typename OPS, typename DATA>
Vector::shared_ptr SimpleVector<TYPE, OPS, DATA>::create( size_t localSize,
                                                          Variable::shared_ptr var )
{
    AMP_MPI comm( AMP_COMM_SELF );
    auto DOFs = AMP::make_shared<AMP::Discretization::DOFManager>( localSize, comm );
    auto ptr  = new SimpleVector<TYPE, OPS, DATA>();
    ptr->setVariable( var );
    ptr->allocate( DOFs->beginDOF(), DOFs->numLocalDOF(), DOFs->numGlobalDOF() );
    ptr->d_comm       = comm;
    ptr->d_DOFManager = DOFs;
    ptr->setCommunicationList(
        AMP::LinearAlgebra::CommunicationList::createEmpty( localSize, comm ) );
    return Vector::shared_ptr( ptr );
}
template<typename TYPE, typename OPS, typename DATA>
Vector::shared_ptr
SimpleVector<TYPE, OPS, DATA>::create( size_t localSize, Variable::shared_ptr var, AMP_MPI comm )
{
    auto DOFs = AMP::make_shared<AMP::Discretization::DOFManager>( localSize, comm );
    auto ptr  = new SimpleVector<TYPE, OPS, DATA>();
    ptr->setVariable( var );
    ptr->allocate( DOFs->beginDOF(), DOFs->numLocalDOF(), DOFs->numGlobalDOF() );
    ptr->d_DOFManager = DOFs;
    ptr->d_comm       = comm;
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
    ptr->allocate( DOFs->beginDOF(), DOFs->numLocalDOF(), DOFs->numGlobalDOF() );
    ptr->d_DOFManager = DOFs;
    ptr->setCommunicationList( commlist );
    ptr->d_comm = DOFs->getComm();
    return Vector::shared_ptr( ptr );
}
template<typename TYPE, typename OPS, typename DATA>
void SimpleVector<TYPE, OPS, DATA>::resize( size_t N )
{
    d_DOFManager = AMP::make_shared<AMP::Discretization::DOFManager>( N, d_comm );
    this->allocate(
        d_DOFManager->beginDOF(), d_DOFManager->numLocalDOF(), d_DOFManager->numGlobalDOF() );
}


/****************************************************************
 * Vector operations                                             *
 ****************************************************************/
template<typename TYPE, typename OPS, typename DATA>
Vector::shared_ptr
SimpleVector<TYPE, OPS, DATA>::cloneVector( const Variable::shared_ptr name ) const
{
    return create( name, d_DOFManager, getCommunicationList() );
}
template<typename TYPE, typename OPS, typename DATA>
void SimpleVector<TYPE, OPS, DATA>::swapVectors( Vector &rhs )
{
    auto x = dynamic_cast<SimpleVector *>( &rhs );
    AMP_INSIST( x != nullptr, "rhs is not a SimpleVector" );
    std::swap( d_comm, d_comm );
    swapData( rhs );
}
template<typename TYPE, typename OPS, typename DATA>
void SimpleVector<TYPE, OPS, DATA>::aliasVector( Vector & )
{
    AMP_ERROR( "Not implemented" );
}
template<typename TYPE, typename OPS, typename DATA>
void SimpleVector<TYPE, OPS, DATA>::assemble()
{
    AMP_ERROR( "Not implemented" );
}
} // namespace LinearAlgebra
} // namespace AMP
