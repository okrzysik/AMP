#ifndef included_AMP_VectorDataCPU_hpp
#define included_AMP_VectorDataCPU_hpp

#include "AMP/vectors/data/VectorDataCPU.h"


namespace AMP::LinearAlgebra {


// Define some specializations
template<>
std::string VectorDataCPU<double>::VectorDataName() const;
template<>
std::string VectorDataCPU<float>::VectorDataName() const;


// Suppresses implicit instantiation below
extern template class VectorDataCPU<double>;
extern template class VectorDataCPU<float>;


/****************************************************************
 * Get the class type                                            *
 ****************************************************************/
template<typename TYPE>
std::string VectorDataCPU<TYPE>::VectorDataName() const
{
    return "VectorDataCPU<" + std::string( typeid( TYPE ).name() ) + ">";
}


/****************************************************************
 * Allocate the data                                             *
 ****************************************************************/
template<typename TYPE>
VectorDataCPU<TYPE>::VectorDataCPU( size_t start, size_t localSize, size_t globalSize )
{
    allocate( start, localSize, globalSize );
}
template<typename TYPE>
void VectorDataCPU<TYPE>::allocate( size_t start, size_t localSize, size_t globalSize )
{
    d_Data.resize( localSize );
    d_startIndex = start;
    d_globalSize = globalSize;
}


/****************************************************************
 * Return basic properties                                       *
 ****************************************************************/
template<typename TYPE>
inline size_t VectorDataCPU<TYPE>::numberOfDataBlocks() const
{
    return 1;
}
template<typename TYPE>
inline size_t VectorDataCPU<TYPE>::sizeOfDataBlock( size_t i ) const
{
    if ( i > 0 )
        return 0;
    return d_Data.size();
}
template<typename TYPE>
inline size_t VectorDataCPU<TYPE>::getLocalSize() const
{
    return d_Data.size();
}
template<typename TYPE>
inline size_t VectorDataCPU<TYPE>::getGlobalSize() const
{
    return d_globalSize;
}
template<typename TYPE>
inline size_t VectorDataCPU<TYPE>::getLocalStartID() const
{
    return d_startIndex;
}
template<typename TYPE>
uint64_t VectorDataCPU<TYPE>::getDataID() const
{
    return reinterpret_cast<uint64_t>( d_Data.data() );
}
template<typename TYPE>
bool VectorDataCPU<TYPE>::isTypeId( size_t hash, size_t ) const
{
    return hash == typeid( TYPE ).hash_code();
}
template<typename TYPE>
size_t VectorDataCPU<TYPE>::sizeofDataBlockType( size_t ) const
{
    return sizeof( TYPE );
}


/****************************************************************
 * Access the raw data blocks                                    *
 ****************************************************************/
template<typename TYPE>
inline void *VectorDataCPU<TYPE>::getRawDataBlockAsVoid( size_t i )
{
    if ( i != 0 ) {
        return 0;
    }
    return d_Data.data();
}
template<typename TYPE>
inline const void *VectorDataCPU<TYPE>::getRawDataBlockAsVoid( size_t i ) const
{
    if ( i != 0 ) {
        return 0;
    }
    return d_Data.data();
}


/****************************************************************
 * Access individual values                                      *
 ****************************************************************/
template<typename TYPE>
inline TYPE &VectorDataCPU<TYPE>::operator[]( size_t i )
{
    return d_Data[i];
}

template<typename TYPE>
inline const TYPE &VectorDataCPU<TYPE>::operator[]( size_t i ) const
{
    return d_Data[i];
}
template<typename TYPE>
inline void VectorDataCPU<TYPE>::setValuesByLocalID( int num, size_t *indices, const double *vals )
{
    for ( int i = 0; i != num; i++ )
        d_Data[indices[i]] = static_cast<TYPE>( vals[i] );
    if ( *d_UpdateState == UpdateState::UNCHANGED )
        *d_UpdateState = UpdateState::LOCAL_CHANGED;
}

template<typename TYPE>
inline void
VectorDataCPU<TYPE>::setLocalValuesByGlobalID( int num, size_t *indices, const double *vals )
{
    for ( int i = 0; i != num; i++ ) {
        AMP_ASSERT( indices[i] >= d_startIndex && indices[i] < d_startIndex + d_Data.size() );
        d_Data[indices[i] - d_startIndex] = static_cast<TYPE>( vals[i] );
    }
    if ( *d_UpdateState == UpdateState::UNCHANGED )
        *d_UpdateState = UpdateState::LOCAL_CHANGED;
}

template<typename TYPE>
inline void VectorDataCPU<TYPE>::addValuesByLocalID( int num, size_t *indices, const double *vals )
{
    for ( int i = 0; i != num; i++ )
        d_Data[indices[i]] += static_cast<TYPE>( vals[i] );
    if ( *d_UpdateState == UpdateState::UNCHANGED )
        *d_UpdateState = UpdateState::LOCAL_CHANGED;
}

template<typename TYPE>
inline void
VectorDataCPU<TYPE>::addLocalValuesByGlobalID( int num, size_t *indices, const double *vals )
{
    for ( int i = 0; i != num; i++ ) {
        AMP_ASSERT( indices[i] >= d_startIndex && indices[i] < d_startIndex + d_Data.size() );
        d_Data[indices[i] - d_startIndex] += static_cast<TYPE>( vals[i] );
    }
    if ( *d_UpdateState == UpdateState::UNCHANGED )
        *d_UpdateState = UpdateState::LOCAL_CHANGED;
}

template<typename TYPE>
inline void
VectorDataCPU<TYPE>::getLocalValuesByGlobalID( int num, size_t *indices, double *vals ) const
{
    for ( int i = 0; i != num; i++ ) {
        AMP_ASSERT( indices[i] >= d_startIndex && indices[i] < d_startIndex + d_Data.size() );
        vals[i] = static_cast<double>( d_Data[indices[i] - d_startIndex] );
    }
}


/****************************************************************
 * Copy raw data                                                 *
 ****************************************************************/
template<typename TYPE>
void VectorDataCPU<TYPE>::putRawData( const double *in )
{
    for ( size_t i = 0; i < d_Data.size(); ++i ) {
        d_Data[i] = static_cast<TYPE>( in[i] );
    }
}

template<typename TYPE>
void VectorDataCPU<TYPE>::copyOutRawData( double *out ) const
{
    for ( size_t i = 0; i < d_Data.size(); ++i ) {
        out[i] = static_cast<double>( d_Data[i] );
    }
}


/****************************************************************
 * Swap raw data                                                 *
 ****************************************************************/
template<typename TYPE>
void VectorDataCPU<TYPE>::swapData( VectorData &rhs )
{
    auto rhs2 = dynamic_cast<VectorDataCPU<TYPE> *>( &rhs );
    AMP_INSIST( rhs2, "Cannot swap with arbitrary VectorData" );
    std::swap( d_CommList, rhs2->d_CommList );
    std::swap( d_UpdateState, rhs2->d_UpdateState );
    std::swap( d_Ghosts, rhs2->d_Ghosts );
    std::swap( d_AddBuffer, rhs2->d_AddBuffer );
    std::swap( d_Data, rhs2->d_Data );
    std::swap( d_startIndex, rhs2->d_startIndex );
    std::swap( d_globalSize, rhs2->d_globalSize );
}


} // namespace AMP::LinearAlgebra

#endif
