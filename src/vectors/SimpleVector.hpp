#include "math.h"

#include "SimpleVector.h"
#include "discretization/DOF_Manager.h"

namespace AMP {
namespace LinearAlgebra {


/****************************************************************
* Constructors                                                  *
****************************************************************/
template <typename T>
SimpleVector<T>::SimpleVector() : Vector(), d_startIndex( 0 ), d_globalSize( 0 )
{
}

template <typename T>
Vector::shared_ptr SimpleVector<T>::create( size_t localSize, const std::string& var )
{
    return create( localSize, AMP::make_shared<Variable>( var ) );
}

template <typename T>
Vector::shared_ptr SimpleVector<T>::create( size_t localSize, Variable::shared_ptr var )
{
    AMP::shared_ptr<SimpleVector<T>> retVal( new SimpleVector<T> );
    retVal->d_startIndex = 0;
    retVal->setVariable( var );
    retVal->d_Data.resize( localSize );
    AMP_MPI comm( AMP_COMM_SELF );
    AMP::Discretization::DOFManager::shared_ptr DOFs(
        new AMP::Discretization::DOFManager( localSize, comm ) );
    retVal->d_DOFManager = DOFs;
    retVal->setCommunicationList(
        AMP::LinearAlgebra::CommunicationList::createEmpty( DOFs->numLocalDOF(), comm ) );
    retVal->d_comm       = comm;
    retVal->d_globalSize = localSize;
    return retVal;
}

template <typename T>
Vector::shared_ptr
SimpleVector<T>::create( size_t localSize, Variable::shared_ptr var, AMP_MPI comm )
{
    AMP::shared_ptr<SimpleVector<T>> retVal( new SimpleVector<T> );
    retVal->d_startIndex = 0;
    retVal->setVariable( var );
    retVal->d_Data.resize( localSize );
    AMP::Discretization::DOFManager::shared_ptr DOFs(
        new AMP::Discretization::DOFManager( localSize, comm ) );
    retVal->d_DOFManager = DOFs;
    retVal->setCommunicationList(
        AMP::LinearAlgebra::CommunicationList::createEmpty( DOFs->numLocalDOF(), comm ) );
    retVal->d_comm       = comm;
    retVal->d_globalSize = comm.sumReduce( localSize );
    return retVal;
}

template <typename T>
Vector::shared_ptr
SimpleVector<T>::create( Variable::shared_ptr var,
                         AMP::Discretization::DOFManager::shared_ptr DOFs,
                         AMP::LinearAlgebra::CommunicationList::shared_ptr commlist )
{
    AMP::shared_ptr<SimpleVector<T>> retVal( new SimpleVector<T> );
    retVal->d_startIndex = DOFs->beginDOF();
    retVal->setVariable( var );
    retVal->d_Data.resize( DOFs->numLocalDOF() );
    retVal->d_DOFManager = DOFs;
    retVal->setCommunicationList( commlist );
    retVal->d_comm       = DOFs->getComm();
    retVal->d_globalSize = DOFs->numGlobalDOF();
    return retVal;
}


/****************************************************************
* Copy raw data                                                 *
****************************************************************/
template <typename T>
void SimpleVector<T>::putRawData( const double *in )
{
    for ( size_t i = 0; i < d_Data.size(); ++i ) {
        d_Data[i] = static_cast<T>( in[i] );
    }
}

template <typename T>
void SimpleVector<T>::copyOutRawData( double *out ) const
{
    for ( size_t i = 0; i < d_Data.size(); ++i ) {
        out[i] = static_cast<double>( d_Data[i] );
    }
}


}
}
