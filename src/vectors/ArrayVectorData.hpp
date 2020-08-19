#ifndef included_ArrayVectorData_hpp
#define included_ArrayVectorData_hpp
#include "AMP/discretization/DOF_Manager.h"

#include "math.h"


namespace AMP {
namespace LinearAlgebra {


/****************************************************************
 * Constructors                                                  *
 ****************************************************************/
template<typename T, typename FUN, typename Allocator>
ArrayVectorData<T, FUN, Allocator>::ArrayVectorData() : VectorData(), d_globalSize( 0 )
{
}

template<typename T, typename FUN, typename Allocator>
std::shared_ptr<VectorData> ArrayVectorData<T, FUN, Allocator>::create( const std::vector<size_t> &localSize )
{
    std::shared_ptr<ArrayVectorData<T, FUN, Allocator>> retVal( new ArrayVectorData<T, FUN, Allocator>() );
    retVal->resize( localSize );
    const auto N = retVal->getArray().length();
    AMP_MPI comm( AMP_COMM_SELF );
    auto DOFs            = std::make_shared<AMP::Discretization::DOFManager>( N, comm );
    retVal->setCommunicationList(
        AMP::LinearAlgebra::CommunicationList::createEmpty( DOFs->numLocalDOF(), comm ) );
    retVal->d_globalSize = N;
    retVal->d_comm       = comm;
    return retVal;
}

template<typename T, typename FUN, typename Allocator>
std::shared_ptr<VectorData> ArrayVectorData<T, FUN, Allocator>::create( const std::vector<size_t> &localSize,
									AMP_MPI comm )
{
    std::shared_ptr<ArrayVectorData<T, FUN, Allocator>> retVal( new ArrayVectorData<T, FUN, Allocator>() );
    retVal->resize( localSize );
    const auto N         = retVal->getArray().length();
    auto DOFs            = std::make_shared<AMP::Discretization::DOFManager>( N, comm );
    retVal->setCommunicationList(
        AMP::LinearAlgebra::CommunicationList::createEmpty( DOFs->numLocalDOF(), comm ) );
    retVal->d_comm       = comm;
    retVal->d_globalSize = comm.sumReduce( N );
    return retVal;
}

template<typename T, typename FUN, typename Allocator>
std::shared_ptr<VectorData>
ArrayVectorData<T, FUN, Allocator>::create( AMP::LinearAlgebra::CommunicationList::shared_ptr commlist )
{
    std::shared_ptr<ArrayVectorData<T, FUN, Allocator>> retVal( new ArrayVectorData<T, FUN, Allocator>() );
    retVal->setCommunicationList( commlist );
    retVal->d_comm = commlist->getComm();
    AMP_ERROR( "This routine is not complete" );
    return retVal;
}

template<typename T, typename FUN, typename Allocator>
inline std::shared_ptr<VectorData>
ArrayVectorData<T, FUN, Allocator>::cloneData( void ) const
{
    const auto &array = this->getArray();
    std::vector<size_t> size( array.size().begin(), array.size().end() );
    return create( size, this->getComm() );
}

template<typename T, typename FUN, typename Allocator>
void ArrayVectorData<T, FUN, Allocator>::swapData( VectorData &rhs )
{
    // get internal arrays
    auto &internalArray = this->getArray();
    auto &otherArray    = dynamic_cast<ArrayVectorData<T, FUN, Allocator> &>( rhs ).getArray();
    // reset views
    internalArray.swap( otherArray );
}

template<typename T, typename FUN, typename Allocator>
void ArrayVectorData<T, FUN, Allocator>::resize( const std::vector<size_t> &localDims )
{
    d_array.resize( localDims );
}


template<typename T, typename FUN, typename Allocator>
void ArrayVectorData<T, FUN, Allocator>::putRawData( const double *buf )
{
    auto &array = this->getArray();
    array.copy( buf );
}

template<typename T, typename FUN, typename Allocator>
void ArrayVectorData<T, FUN, Allocator>::copyOutRawData( double *buf ) const
{
    auto &array = this->getArray();
    array.copyTo( buf );
}
} // namespace LinearAlgebra
} // namespace AMP
#endif
