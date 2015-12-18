#include "math.h"

#include "ArrayVector.h"
#include "discretization/DOF_Manager.h"

namespace AMP {
namespace LinearAlgebra {


/****************************************************************
* Constructors                                                  *
****************************************************************/
template <typename T>
ArrayVector<T>::ArrayVector() : SimpleVector<T>()
{
}

template <typename T>
Vector::shared_ptr ArrayVector<T>::create( const std::vector<size_t> &localSize,
                                           Variable::shared_ptr var )
{
    AMP::shared_ptr<ArrayVector<T>> retVal( new ArrayVector<T> );
    retVal->setVariable( var );

    // allocate space for the vector
    size_t N = 1;
    for ( auto s : localSize )
        N *= s;
    retVal->resize( N );
    // extract pointers to the internal vector and array
    // do not use 'auto' in place of AMP::Array<T> &
    // and std::vector<T> & as these result in the
    // copy constructor being called!!
    AMP::Array<T> &internalArray = retVal->getArray();
    std::vector<T> &internalVec   = const_cast<std::vector<T> &>(retVal->getData());
    // set the data pointer for the array to point to the std:vector data
    internalArray.viewRaw( localSize, internalVec.data() );
    AMP_ASSERT(internalArray.size()==localSize);
    AMP_ASSERT(internalArray.ndim()==localSize.size());

    AMP_MPI comm( AMP_COMM_SELF );
    AMP::Discretization::DOFManager::shared_ptr DOFs(
        new AMP::Discretization::DOFManager( N, comm ) );
    retVal->d_DOFManager = DOFs;
    retVal->setCommunicationList(
        AMP::LinearAlgebra::CommunicationList::createEmpty( DOFs->numLocalDOF(), comm ) );
    retVal->d_comm       = comm;
    retVal->d_globalSize = N;
    return retVal;
}

template <typename T>
Vector::shared_ptr ArrayVector<T>::create( const std::vector<size_t> &localSize,
                                           Variable::shared_ptr var,
                                           AMP_MPI comm )
{
    AMP::shared_ptr<ArrayVector<T>> retVal( new ArrayVector<T> );
    retVal->setVariable( var );

    size_t N = 1;
    for ( auto s : localSize )
        N *= s;
    retVal->resize( N );
    // extract pointers to the internal vector and array
    auto internalArray = retVal->getArray();
    auto internalVec   = retVal->getData();
    // set the data pointer for the array to point to the std:vector data
    internalArray.viewRaw( localSize, internalVec.data() );

    AMP::Discretization::DOFManager::shared_ptr DOFs(
        new AMP::Discretization::DOFManager( N, comm ) );
    retVal->d_DOFManager = DOFs;
    retVal->setCommunicationList(
        AMP::LinearAlgebra::CommunicationList::createEmpty( DOFs->numLocalDOF(), comm ) );
    retVal->d_comm       = comm;
    retVal->d_globalSize = comm.sumReduce( N );
    return retVal;
}

template <typename T>
Vector::shared_ptr
ArrayVector<T>::create( Variable::shared_ptr var,
                        AMP::Discretization::DOFManager::shared_ptr DOFs,
                        AMP::LinearAlgebra::CommunicationList::shared_ptr commlist )
{
    AMP::shared_ptr<ArrayVector<T>> retVal( new ArrayVector<T> );
    retVal->d_startIndex = DOFs->beginDOF();
    retVal->setVariable( var );
    retVal->resize( DOFs->numLocalDOF() );
    retVal->d_DOFManager = DOFs;
    retVal->setCommunicationList( commlist );
    retVal->d_comm       = DOFs->getComm();
    retVal->d_globalSize = DOFs->numGlobalDOF();
    return retVal;
}

/****************************************************************
* Copy vector                                                   *
****************************************************************/
template <typename T>
void ArrayVector<T>::copyVector( Vector::const_shared_ptr src_vec )
{
    SimpleVector<T>::copyVector( src_vec );
}

template <typename T>
inline Vector::shared_ptr ArrayVector<T>::cloneVector( const Variable::shared_ptr name ) const
{
    return create( name, this->getDOFManager(), this->getCommunicationList() );
}

template <typename T>
void ArrayVector<T>::swapVectors( Vector &rhs )
{
    // get information on array dimensions
    auto internalArray = this->getArray();
    auto internalVec   = this->getData();
    auto mySize        = internalArray.size();

    auto otherArray = rhs.castTo<ArrayVector<T>>().getArray();
    auto otherVec   = rhs.castTo<ArrayVector<T>>().getData();
    auto otherSize  = otherArray.size();
    // for now assume arrays are same size (do we need/want to?)
    AMP_ASSERT( mySize == otherSize );

    // swap data vectors
    internalVec.swap( otherVec );
    // reset views
    internalArray.viewRaw( otherSize, internalVec.data() );
    otherArray.viewRaw( mySize, internalVec.data() );
}

template <typename T>
void ArrayVector<T>::aliasVector( Vector & )
{
    AMP_ERROR( "Not implemented" );
}
}
}
