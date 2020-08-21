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
    ptr->setDOFManager(DOFs);
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
    ptr->setDOFManager(DOFs);
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
    ptr->setDOFManager(DOFs);
    auto comm = DOFs->getComm();
    ptr->setComm( comm );
#else
    ptr->allocate( DOFs->beginDOF(), DOFs->numLocalDOF(), DOFs->numGlobalDOF() );
    ptr->d_DOFManager = DOFs;
    ptr->d_comm = DOFs->getComm();
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

/****************************************************************
 * VectorData operations -- will move to Vector eventually      *
 ****************************************************************/
  
template<typename TYPE, typename OPS, typename DATA>
std::string SimpleVector<TYPE, OPS, DATA>::VectorDataName() const
{
  AMP_ASSERT(d_VectorData.get()!=this);
  return d_VectorData->VectorDataName();
}

template<typename TYPE, typename OPS, typename DATA>
size_t SimpleVector<TYPE, OPS, DATA>::numberOfDataBlocks() const
{
  AMP_ASSERT(d_VectorData.get()!=this);
  return d_VectorData->numberOfDataBlocks();
}
  
template<typename TYPE, typename OPS, typename DATA>
size_t SimpleVector<TYPE, OPS, DATA>::sizeOfDataBlock( size_t i ) const
{
  AMP_ASSERT(d_VectorData.get()!=this);
  return d_VectorData->sizeOfDataBlock(i);
}
 
template<typename TYPE, typename OPS, typename DATA>
void SimpleVector<TYPE, OPS, DATA>::putRawData( const double *buf )
{
  AMP_ASSERT(d_VectorData.get()!=this);
  d_VectorData->putRawData(buf);
}
  
template<typename TYPE, typename OPS, typename DATA>
void SimpleVector<TYPE, OPS, DATA>::copyOutRawData( double *buf ) const
{
  AMP_ASSERT(d_VectorData.get()!=this);
  d_VectorData->copyOutRawData(buf);
}
 
template<typename TYPE, typename OPS, typename DATA>
size_t SimpleVector<TYPE, OPS, DATA>::getLocalSize() const
{
  AMP_ASSERT(d_VectorData.get()!=this);
  return d_VectorData->getLocalSize();
}
  
template<typename TYPE, typename OPS, typename DATA>
size_t SimpleVector<TYPE, OPS, DATA>::getGlobalSize() const
{
  AMP_ASSERT(d_VectorData.get()!=this);
  return d_VectorData->getGlobalSize();
}
  
template<typename TYPE, typename OPS, typename DATA>
size_t SimpleVector<TYPE, OPS, DATA>::getLocalStartID() const
{
  AMP_ASSERT(d_VectorData.get()!=this);
  return d_VectorData->getLocalStartID();
}
  
template<typename TYPE, typename OPS, typename DATA>
void SimpleVector<TYPE, OPS, DATA>::setValuesByLocalID( int num, size_t *indices, const double *vals )
{
  AMP_ASSERT(d_VectorData.get()!=this);
  d_VectorData->setValuesByLocalID( num, indices, vals );
}
  
template<typename TYPE, typename OPS, typename DATA>
void SimpleVector<TYPE, OPS, DATA>::setLocalValuesByGlobalID( int num, size_t *indices, const double *vals )
{
  AMP_ASSERT(d_VectorData.get()!=this);
  d_VectorData->setLocalValuesByGlobalID( num, indices, vals );
}
  
template<typename TYPE, typename OPS, typename DATA>
void SimpleVector<TYPE, OPS, DATA>::addValuesByLocalID( int num, size_t *indices, const double *vals )
{
  AMP_ASSERT(d_VectorData.get()!=this);
  d_VectorData->addValuesByLocalID( num, indices, vals );
}
  
template<typename TYPE, typename OPS, typename DATA>
void SimpleVector<TYPE, OPS, DATA>::addLocalValuesByGlobalID( int num, size_t *indices, const double *vals )
{
  AMP_ASSERT(d_VectorData.get()!=this);
  d_VectorData->addLocalValuesByGlobalID( num, indices, vals );
}
  
template<typename TYPE, typename OPS, typename DATA>
void SimpleVector<TYPE, OPS, DATA>::getLocalValuesByGlobalID( int num, size_t *indices, double *vals ) const
{
  AMP_ASSERT(d_VectorData.get()!=this);
  d_VectorData->getLocalValuesByGlobalID( num, indices, vals );
}
  
template<typename TYPE, typename OPS, typename DATA>
uint64_t SimpleVector<TYPE, OPS, DATA>::getDataID() const
{
  AMP_ASSERT(d_VectorData.get()!=this);
  return d_VectorData->getDataID();
}
  
template<typename TYPE, typename OPS, typename DATA>
void *SimpleVector<TYPE, OPS, DATA>::getRawDataBlockAsVoid( size_t i )
{
  AMP_ASSERT(d_VectorData.get()!=this);
  return d_VectorData->getRawDataBlockAsVoid(i);
}
  
template<typename TYPE, typename OPS, typename DATA>
const void *SimpleVector<TYPE, OPS, DATA>::getRawDataBlockAsVoid( size_t i ) const
{
  AMP_ASSERT(d_VectorData.get()!=this);
  return d_VectorData->getRawDataBlockAsVoid(i);
}
  
template<typename TYPE, typename OPS, typename DATA>
size_t SimpleVector<TYPE, OPS, DATA>::sizeofDataBlockType( size_t i ) const
{
  AMP_ASSERT(d_VectorData.get()!=this);
  return d_VectorData->sizeofDataBlockType(i);
}
  
template<typename TYPE, typename OPS, typename DATA>
bool SimpleVector<TYPE, OPS, DATA>::isTypeId( size_t hash, size_t block ) const
{
  AMP_ASSERT(d_VectorData.get()!=this);
  return d_VectorData->isTypeId(hash, block);
}
  
template<typename TYPE, typename OPS, typename DATA>
void SimpleVector<TYPE, OPS, DATA>::swapData( VectorData &rhs )
{
  AMP_ASSERT(d_VectorData.get()!=this);
  d_VectorData->swapData(rhs);
}
  
template<typename TYPE, typename OPS, typename DATA>
std::shared_ptr<VectorData> SimpleVector<TYPE, OPS, DATA>::cloneData() const
{
  AMP_ASSERT(d_VectorData.get()!=this);
  return d_VectorData->cloneData();
}

// temporary functions
template<typename TYPE, typename OPS, typename DATA>
void SimpleVector<TYPE, OPS, DATA>::allocateVectorData( size_t localSize, size_t numLocal, size_t numGlobal )
{
  d_VectorData = std::make_shared<DATA>(localSize, numLocal, numGlobal);
  AMP_ASSERT(d_VectorData.get()!=this);
}

template<typename TYPE, typename OPS, typename DATA>
void SimpleVector<TYPE, OPS, DATA>::setDOFManager( std::shared_ptr<AMP::Discretization::DOFManager> dofManager)
{
  AMP_ASSERT(d_VectorData.get()!=this);
  d_DOFManager = dofManager;
}

template<typename TYPE, typename OPS, typename DATA>
void SimpleVector<TYPE, OPS, DATA>::setCommunicationList( CommunicationList::shared_ptr comm )
{
  AMP_ASSERT(d_VectorData.get()!=this);
  d_VectorData->setCommunicationList(comm);
}
 
} // namespace LinearAlgebra
} // namespace AMP
