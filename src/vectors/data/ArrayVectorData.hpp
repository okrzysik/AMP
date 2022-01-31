#ifndef included_AMP_ArrayVectorData_hpp
#define included_AMP_ArrayVectorData_hpp
#include "AMP/discretization/DOF_Manager.h"

#include "math.h"


namespace AMP::LinearAlgebra {


/****************************************************************
 * Constructors                                                  *
 ****************************************************************/
template<typename T, typename FUN, typename Allocator>
std::shared_ptr<ArrayVectorData<T, FUN, Allocator>>
ArrayVectorData<T, FUN, Allocator>::create( const ArraySize &localSize )
{
    return create( localSize, { 0, 0, 0 }, AMP_COMM_SELF );
}
template<typename T, typename FUN, typename Allocator>
std::shared_ptr<ArrayVectorData<T, FUN, Allocator>> ArrayVectorData<T, FUN, Allocator>::create(
    const ArraySize &localSize, const ArraySize &blockIndex, AMP_MPI comm )
{
    size_t N_blocks[3];
    for ( int d = 0; d < 3; d++ ) {
        N_blocks[d] = comm.maxReduce( blockIndex[d] + 1 );
        AMP_ASSERT( comm.maxReduce( localSize[d] ) == localSize[d] );
    }
    AMP_ASSERT( N_blocks[0] * N_blocks[1] * N_blocks[2] == (size_t) comm.getSize() );
    size_t blockOffset =
        blockIndex[0] + N_blocks[0] * ( blockIndex[1] + blockIndex[2] * N_blocks[1] );
    std::shared_ptr<ArrayVectorData<T, FUN, Allocator>> retVal(
        new ArrayVectorData<T, FUN, Allocator>() );
    retVal->d_array.resize( localSize );
    retVal->d_comm            = comm;
    retVal->d_blockIndex      = blockIndex;
    retVal->d_globalArraySize = { localSize[0] * N_blocks[0],
                                  localSize[1] * N_blocks[1],
                                  localSize[2] * N_blocks[2] };
    retVal->d_offset          = blockOffset * localSize.length();
    retVal->setCommunicationList(
        AMP::LinearAlgebra::CommunicationList::createEmpty( localSize.length(), comm ) );
    retVal->d_localSize  = localSize.length();
    retVal->d_globalSize = retVal->d_globalArraySize.length();
    retVal->d_localStart = retVal->d_CommList->getStartGID();
    return retVal;
}


/****************************************************************
 * Clone/Swap data                                               *
 ****************************************************************/
template<typename T, typename FUN, typename Allocator>
inline std::shared_ptr<VectorData> ArrayVectorData<T, FUN, Allocator>::cloneData() const
{
    auto retVal               = std::make_shared<ArrayVectorData<T, FUN, Allocator>>();
    retVal->d_array           = d_array;
    retVal->d_comm            = d_comm;
    retVal->d_blockIndex      = d_blockIndex;
    retVal->d_globalArraySize = d_globalSize;
    retVal->d_localSize       = d_localSize;
    retVal->d_globalSize      = d_globalSize;
    retVal->d_localStart      = d_localStart;
    retVal->setCommunicationList( getCommunicationList() );
    return retVal;
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


/****************************************************************
 * Resize data                                                   *
 ****************************************************************/
template<typename T, typename FUN, typename Allocator>
void ArrayVectorData<T, FUN, Allocator>::resize( const ArraySize &localDims )
{
    AMP_ASSERT( getComm().getSize() == 1 );
    d_array.resize( localDims );
    d_globalArraySize = localDims;
    setCommunicationList(
        AMP::LinearAlgebra::CommunicationList::createEmpty( localDims.length(), getComm() ) );
    d_localSize  = localDims.length();
    d_globalSize = d_globalArraySize.length();
    d_localStart = d_CommList->getStartGID();
}


/****************************************************************
 * Get/Set raw data                                              *
 ****************************************************************/
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


/****************************************************************
 * Get/Set values                                                *
 ****************************************************************/
template<typename T, typename FUN, typename Allocator>
void ArrayVectorData<T, FUN, Allocator>::setValuesByLocalID( int num,
                                                             size_t *indices,
                                                             const double *vals )
{
    for ( int i = 0; i < num; i++ )
        d_array( indices[i] ) = vals[i];
}
template<typename T, typename FUN, typename Allocator>
void ArrayVectorData<T, FUN, Allocator>::addValuesByLocalID( int num,
                                                             size_t *indices,
                                                             const double *vals )
{
    for ( int i = 0; i < num; i++ )
        d_array( indices[i] ) += vals[i];
}
template<typename T, typename FUN, typename Allocator>
void ArrayVectorData<T, FUN, Allocator>::getLocalValuesByGlobalID( int num,
                                                                   size_t *indices,
                                                                   double *vals ) const
{
    for ( int i = 0; i < num; i++ )
        vals[i] = d_array( indices[i] - d_offset );
}
template<typename T, typename FUN, typename Allocator>
void ArrayVectorData<T, FUN, Allocator>::setLocalValuesByGlobalID( int num,
                                                                   size_t *indices,
                                                                   const double *vals )
{
    for ( int i = 0; i < num; i++ )
        d_array( indices[i] - d_offset ) = vals[i];
}
template<typename T, typename FUN, typename Allocator>
void ArrayVectorData<T, FUN, Allocator>::addLocalValuesByGlobalID( int num,
                                                                   size_t *indices,
                                                                   const double *vals )
{
    for ( int i = 0; i < num; i++ )
        d_array( indices[i] - d_offset ) += vals[i];
}


} // namespace AMP::LinearAlgebra
#endif
