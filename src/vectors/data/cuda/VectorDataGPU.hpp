#ifndef included_AMP_VectorDataGPU_hpp
#define included_AMP_VectorDataGPU_hpp

#include "AMP/vectors/data/cuda/VectorDataGPU.h"


#include <cstring>
#include <cuda.h>
#include <cuda_runtime_api.h>
#include <string>


namespace AMP::LinearAlgebra {


/****************************************************************
 * Get the class type                                            *
 ****************************************************************/
template<typename TYPE>
std::string VectorDataGPU<TYPE>::VectorDataName() const
{
    return "VectorDataGPU<" + std::string( typeid( TYPE ).name() ) + ">";
}


/****************************************************************
 * Allocate the data                                             *
 ****************************************************************/
template<typename TYPE>
VectorDataGPU<TYPE>::VectorDataGPU( size_t start, size_t localSize, size_t globalSize )
{
    allocate( start, localSize, globalSize );
}
template<typename TYPE>
void VectorDataGPU<TYPE>::allocate( size_t start, size_t localSize, size_t globalSize )
{
    d_localSize  = localSize;
    d_globalSize = globalSize;
    d_localStart = start;
    if ( d_localSize > 0 ) {
        cudaMallocManaged( (void **) &d_Data, localSize * sizeof( TYPE ), cudaMemAttachGlobal );
        AMP_INSIST( d_Data, "Failed to allocate memory on device" );
    }
}
template<typename TYPE>
VectorDataGPU<TYPE>::~VectorDataGPU()
{
    cudaFree( d_Data );
}


/****************************************************************
 * Clone the data                                                *
 ****************************************************************/
template<typename TYPE>
std::shared_ptr<VectorData> VectorDataGPU<TYPE>::cloneData() const
{
    auto retVal = std::make_shared<VectorDataGPU<TYPE>>( d_localStart, d_localSize, d_globalSize );
    retVal->setCommunicationList( getCommunicationList() );
    return retVal;
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
uint64_t VectorDataGPU<TYPE>::getDataID() const
{
    return reinterpret_cast<uint64_t>( d_Data );
}
template<typename TYPE>
typeID VectorDataGPU<TYPE>::getType( size_t i ) const
{
    constexpr auto type = getTypeID<TYPE>();
    return type;
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
inline void VectorDataGPU<TYPE>::setValuesByLocalID( size_t num,
                                                     const size_t *indices,
                                                     const void *vals,
                                                     const typeID &id )
{
    if ( id == getTypeID<TYPE>() ) {
        auto data = reinterpret_cast<const TYPE *>( vals );
        for ( size_t i = 0; i < num; i++ )
            d_Data[indices[i]] = data[i];
    } else if ( id == getTypeID<double>() ) {
        auto data = reinterpret_cast<const double *>( vals );
        for ( size_t i = 0; i < num; ++i )
            d_Data[indices[i]] = static_cast<TYPE>( data[i] );
    } else {
        AMP_ERROR( "Conversion not supported yet" );
    }
    if ( *d_UpdateState == UpdateState::UNCHANGED )
        *d_UpdateState = UpdateState::LOCAL_CHANGED;
}
template<typename TYPE>
inline void VectorDataGPU<TYPE>::addValuesByLocalID( size_t num,
                                                     const size_t *indices,
                                                     const void *vals,
                                                     const typeID &id )
{
    if ( id == getTypeID<TYPE>() ) {
        auto data = reinterpret_cast<const TYPE *>( vals );
        for ( size_t i = 0; i < num; i++ )
            d_Data[indices[i]] += data[i];
    } else if ( id == getTypeID<double>() ) {
        auto data = reinterpret_cast<const double *>( vals );
        for ( size_t i = 0; i < num; ++i )
            d_Data[indices[i]] += static_cast<TYPE>( data[i] );
    } else {
        AMP_ERROR( "Conversion not supported yet" );
    }
    if ( *d_UpdateState == UpdateState::UNCHANGED )
        *d_UpdateState = UpdateState::LOCAL_CHANGED;
}
template<typename TYPE>
inline void VectorDataGPU<TYPE>::getValuesByLocalID( size_t num,
                                                     const size_t *indices,
                                                     void *vals,
                                                     const typeID &id ) const
{
    if ( id == getTypeID<TYPE>() ) {
        auto data = reinterpret_cast<TYPE *>( vals );
        for ( size_t i = 0; i < num; i++ )
            data[i] = d_Data[indices[i]];
    } else if ( id == getTypeID<double>() ) {
        auto data = reinterpret_cast<double *>( vals );
        for ( size_t i = 0; i < num; ++i )
            data[i] = d_Data[indices[i]];
    } else {
        AMP_ERROR( "Conversion not supported yet" );
    }
}


/****************************************************************
 * Copy raw data                                                 *
 ****************************************************************/
template<typename TYPE>
void VectorDataGPU<TYPE>::putRawData( const void *in, const typeID &id )
{
    if ( id == getTypeID<TYPE>() ) {
        memcpy( d_Data, in, d_localSize * sizeof( TYPE ) );
    } else if ( id == getTypeID<double>() ) {
        auto data = reinterpret_cast<const double *>( in );
        for ( size_t i = 0; i < d_localSize; ++i )
            d_Data[i] = static_cast<TYPE>( data[i] );
    } else {
        AMP_ERROR( "Conversion not supported yet" );
    }
}

template<typename TYPE>
void VectorDataGPU<TYPE>::getRawData( void *out, const typeID &id ) const
{
    if ( id == getTypeID<TYPE>() ) {
        memcpy( out, d_Data, d_localSize * sizeof( TYPE ) );
    } else if ( id == getTypeID<double>() ) {
        auto data = reinterpret_cast<double *>( out );
        for ( size_t i = 0; i < d_localSize; ++i )
            data[i] = static_cast<double>( d_Data[i] );
    } else {
        AMP_ERROR( "Conversion not supported yet" );
    }
}


/****************************************************************
 * Swap raw data                                                 *
 ****************************************************************/
template<typename TYPE>
void VectorDataGPU<TYPE>::swapData( VectorData &rhs )
{
    auto rhs2 = dynamic_cast<VectorDataGPU<TYPE> *>( &rhs );
    AMP_INSIST( rhs2, "Cannot swap with arbitrary VectorData" );
    std::swap( d_CommList, rhs2->d_CommList );
    std::swap( d_UpdateState, rhs2->d_UpdateState );
    std::swap( d_Ghosts, rhs2->d_Ghosts );
    std::swap( d_AddBuffer, rhs2->d_AddBuffer );
    std::swap( d_Data, rhs2->d_Data );
    std::swap( d_localSize, rhs2->d_localSize );
    std::swap( d_globalSize, rhs2->d_globalSize );
    std::swap( d_localStart, rhs2->d_localStart );
}


} // namespace AMP::LinearAlgebra

#endif
