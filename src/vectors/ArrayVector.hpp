#include "math.h"

#include "ArrayVector.h"
#include "discretization/DOF_Manager.h"

namespace AMP {
namespace LinearAlgebra {


/****************************************************************
* Constructors                                                  *
****************************************************************/
template <typename T>
ArrayVector<T>::ArrayVector (): 
    SimpleVector<T>()
{
}

template <typename T>
Vector::shared_ptr  ArrayVector<T>::create (const std::vector<size_t> &localSize , Variable::shared_ptr var )
{
    size_t N = 1;
    for(auto s: localSize) N*=s;

    AMP::shared_ptr<ArrayVector<T> > retVal( new ArrayVector<T> );
    retVal->setVariable ( var );
    retVal->resize ( N );
    auto internalArray = retVal->getArray();
    auto internalVec = retVal->getData();
    //    internalArray.viewRaw(localSize, (T*)&internalVec[0]);
    AMP_MPI comm(AMP_COMM_SELF);
    AMP::Discretization::DOFManager::shared_ptr DOFs( new AMP::Discretization::DOFManager( N, comm ) );
    retVal->d_DOFManager = DOFs;
    retVal->setCommunicationList( AMP::LinearAlgebra::CommunicationList::createEmpty( DOFs->numLocalDOF(), comm ) );
    retVal->d_comm = comm;
    retVal->d_globalSize = N;
    return retVal;
}

template <typename T>
Vector::shared_ptr  ArrayVector<T>::create ( const std::vector<size_t>& localSize, Variable::shared_ptr var, AMP_MPI comm )
{
    size_t N = 1;
    for(auto s: localSize) N*=s;

    AMP::shared_ptr<ArrayVector<T> > retVal( new ArrayVector<T> );
    retVal->setVariable ( var );
    retVal->resize ( N );
    AMP::Discretization::DOFManager::shared_ptr DOFs( new AMP::Discretization::DOFManager( N, comm ) );
    retVal->d_DOFManager = DOFs;
    retVal->setCommunicationList( AMP::LinearAlgebra::CommunicationList::createEmpty( DOFs->numLocalDOF(), comm ) );
    retVal->d_comm = comm;
    retVal->d_globalSize = comm.sumReduce(N);
    return retVal;
}

template <typename T>
Vector::shared_ptr  ArrayVector<T>::create ( Variable::shared_ptr var,
    AMP::Discretization::DOFManager::shared_ptr DOFs, 
    AMP::LinearAlgebra::CommunicationList::shared_ptr commlist )
{
    AMP::shared_ptr<ArrayVector<T> > retVal( new ArrayVector<T> );
    retVal->d_startIndex = DOFs->beginDOF();
    retVal->setVariable ( var );
    retVal->resize (  DOFs->numLocalDOF() );
    retVal->d_DOFManager = DOFs;
    retVal->setCommunicationList( commlist );
    retVal->d_comm = DOFs->getComm();
    retVal->d_globalSize = DOFs->numGlobalDOF();
    return retVal;
}
       
/****************************************************************
* Copy vector                                                   *
****************************************************************/
template <typename T>
void ArrayVector<T>::copyVector( Vector::const_shared_ptr src_vec )
{
    SimpleVector<T>::copyVector(src_vec);
    // add code for array here
}

template <typename T>
inline
Vector::shared_ptr ArrayVector<T>::cloneVector(const Variable::shared_ptr name) const
{
    return create ( name, this->getDOFManager(), this->getCommunicationList() );
}

template <typename T>
void ArrayVector<T>::swapVectors(Vector &rhs)
{
    this->getData().swap ( rhs.castTo<ArrayVector<T> >().getData() );
}

template <typename T>
void ArrayVector<T>::aliasVector(Vector &)
{
    AMP_ERROR( "Not implemented" );
}

}
}

