#ifndef included_AMP_VectorDataGPU_hpp
#define included_AMP_VectorDataGPU_hpp

#include "vectors/data/cuda/VectorDataGPU.h"


#include <cuda.h>
#include <cuda_runtime_api.h>


namespace AMP {
namespace LinearAlgebra {


extern template class VectorDataGPU<double>; // Suppresses implicit instantiation below --
extern template class VectorDataGPU<float>;  // Suppresses implicit instantiation below --


/****************************************************************
 * Allocate the data                                             *
 ****************************************************************/
template<typename TYPE>
void VectorDataGPU<TYPE>::allocate( size_t start, size_t localSize, size_t globalSize )
{
    cudaMallocManaged( (void **) &d_Data, localSize * sizeof( TYPE ), cudaMemAttachGlobal );
    AMP_INSIST( d_Data, "Failed to allocate memory on device" );
    d_startIndex = start;
    d_localSize  = localSize;
    d_globalSize = globalSize;
}
template<typename TYPE>
VectorDataGPU<TYPE>::~VectorDataGPU()
{
    cudaFree( d_Data );
}


/****************************************************************
 * Return basic properties                                       *
 ****************************************************************/
template<typename TYPE>
size_t VectorDataGPU<TYPE>::numberOfDataBlocks() const
{
    return 1;
}
template<typename TYPE>
size_t VectorDataGPU<TYPE>::sizeOfDataBlock( size_t i ) const
{
    if ( i > 0 )
        return 0;
    return d_localSize;
}
template<typename TYPE>
size_t VectorDataGPU<TYPE>::getLocalSize() const
{
    return d_localSize;
}
template<typename TYPE>
size_t VectorDataGPU<TYPE>::getGlobalSize() const
{
    return d_globalSize;
}
template<typename TYPE>
uint64_t VectorDataGPU<TYPE>::getDataID() const
{
    return reinterpret_cast<uint64_t>( d_Data );
}
template<typename TYPE>
bool VectorDataGPU<TYPE>::isTypeId( size_t hash, size_t ) const
{
    return hash == typeid( TYPE ).hash_code();
}
template<typename TYPE>
size_t VectorDataGPU<TYPE>::sizeofDataBlockType( size_t ) const
{
    return sizeof( TYPE );
}


/****************************************************************
 * Access the raw data blocks                                    *
 ****************************************************************/
template<typename TYPE>
inline void *VectorDataGPU<TYPE>::getRawDataBlockAsVoid( size_t i )
{
    if ( i != 0 ) {
        return 0;
    }
    return d_Data;
}
template<typename TYPE>
inline const void *VectorDataGPU<TYPE>::getRawDataBlockAsVoid( size_t i ) const
{
    if ( i != 0 ) {
        return 0;
    }
    return d_Data;
}


/****************************************************************
 * Access individual values                                      *
 ****************************************************************/
template<typename TYPE>
inline TYPE &VectorDataGPU<TYPE>::operator[]( size_t i )
{
    return d_Data[i];
}

template<typename TYPE>
inline const TYPE &VectorDataGPU<TYPE>::operator[]( size_t i ) const
{
    return d_Data[i];
}
template<typename TYPE>
inline void VectorDataGPU<TYPE>::setValuesByLocalID( int num, size_t *indices, const double *vals )
{
    for ( int i = 0; i != num; i++ )
        d_Data[indices[i]] = static_cast<TYPE>( vals[i] );
    if ( *d_UpdateState == UpdateState::UNCHANGED )
        *d_UpdateState = UpdateState::LOCAL_CHANGED;
}

template<typename TYPE>
inline void
VectorDataGPU<TYPE>::setLocalValuesByGlobalID( int num, size_t *indices, const double *vals )
{
    for ( int i = 0; i != num; i++ ) {
        AMP_ASSERT( indices[i] >= d_startIndex && indices[i] < d_startIndex + d_localSize );
        d_Data[indices[i] - d_startIndex] = static_cast<TYPE>( vals[i] );
    }
    if ( *d_UpdateState == UpdateState::UNCHANGED )
        *d_UpdateState = UpdateState::LOCAL_CHANGED;
}

template<typename TYPE>
inline void VectorDataGPU<TYPE>::addValuesByLocalID( int num, size_t *indices, const double *vals )
{
    for ( int i = 0; i != num; i++ )
        d_Data[indices[i]] += static_cast<TYPE>( vals[i] );
    if ( *d_UpdateState == UpdateState::UNCHANGED )
        *d_UpdateState = UpdateState::LOCAL_CHANGED;
}

template<typename TYPE>
inline void
VectorDataGPU<TYPE>::addLocalValuesByGlobalID( int num, size_t *indices, const double *vals )
{
    for ( int i = 0; i != num; i++ ) {
        AMP_ASSERT( indices[i] >= d_startIndex && indices[i] < d_startIndex + d_localSize );
        d_Data[indices[i] - d_startIndex] += static_cast<TYPE>( vals[i] );
    }
    if ( *d_UpdateState == UpdateState::UNCHANGED )
        *d_UpdateState = UpdateState::LOCAL_CHANGED;
}

template<typename TYPE>
inline void
VectorDataGPU<TYPE>::getLocalValuesByGlobalID( int num, size_t *indices, double *vals ) const
{
    for ( int i = 0; i != num; i++ ) {
        AMP_ASSERT( indices[i] >= d_startIndex && indices[i] < d_startIndex + d_localSize );
        vals[i] = static_cast<double>( d_Data[indices[i] - d_startIndex] );
    }
}


/****************************************************************
 * Copy raw data                                                 *
 ****************************************************************/
template<typename TYPE>
void VectorDataGPU<TYPE>::putRawData( const double *in )
{
    for ( size_t i = 0; i < d_localSize; ++i ) {
        d_Data[i] = static_cast<TYPE>( in[i] );
    }
}

template<typename TYPE>
void VectorDataGPU<TYPE>::copyOutRawData( double *out ) const
{
    for ( size_t i = 0; i < d_localSize; ++i ) {
        out[i] = static_cast<double>( d_Data[i] );
    }
}


} // namespace LinearAlgebra
} // namespace AMP

#endif
