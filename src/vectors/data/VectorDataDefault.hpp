#ifndef included_AMP_VectorDataDefault_hpp
#define included_AMP_VectorDataDefault_hpp

#include "AMP/IO/RestartManager.h"
#include "AMP/vectors/data/VectorDataDefault.h"

#include <cstring>


namespace AMP::LinearAlgebra {


// Suppresses implicit instantiation below
extern template class VectorDataDefault<double>;
extern template class VectorDataDefault<float>;


/****************************************************************
 * Get the class type                                            *
 ****************************************************************/
template<typename TYPE, class Allocator>
std::string VectorDataDefault<TYPE, Allocator>::VectorDataName() const
{
    constexpr typeID id = getTypeID<TYPE>();
    return "VectorDataDefault<" + std::string( id.name ) + ">";
}


/****************************************************************
 * Allocate the data                                             *
 ****************************************************************/
template<typename TYPE, class Allocator>
VectorDataDefault<TYPE, Allocator>::VectorDataDefault( size_t start,
                                                       size_t localSize,
                                                       size_t globalSize,
                                                       const Allocator &alloc )
    : d_alloc( alloc )
{
    static_assert( std::is_same_v<typename Allocator::value_type, TYPE> );
    d_localSize  = localSize;
    d_globalSize = globalSize;
    d_localStart = start;
    d_data       = d_alloc.allocate( localSize );
    for ( size_t i = 0; i < localSize; ++i )
        new ( d_data + i ) TYPE();
}
template<typename TYPE, class Allocator>
VectorDataDefault<TYPE, Allocator>::~VectorDataDefault()
{
    for ( size_t i = 0; i < d_localSize; ++i )
        d_data[i].~TYPE();
    d_alloc.deallocate( d_data, d_localSize );
}


/****************************************************************
 * Clone the data                                                *
 ****************************************************************/
template<typename TYPE, class Allocator>
std::shared_ptr<VectorData> VectorDataDefault<TYPE, Allocator>::cloneData() const
{
    auto retVal = std::make_shared<VectorDataDefault<TYPE, Allocator>>(
        d_localStart, d_localSize, d_globalSize );
    retVal->setCommunicationList( getCommunicationList() );
    return retVal;
}


/****************************************************************
 * Return basic properties                                       *
 ****************************************************************/
template<typename TYPE, class Allocator>
inline size_t VectorDataDefault<TYPE, Allocator>::numberOfDataBlocks() const
{
    return 1;
}
template<typename TYPE, class Allocator>
inline size_t VectorDataDefault<TYPE, Allocator>::sizeOfDataBlock( size_t i ) const
{
    if ( i > 0 )
        return 0;
    return d_localSize;
}
template<typename TYPE, class Allocator>
uint64_t VectorDataDefault<TYPE, Allocator>::getDataID() const
{
    return reinterpret_cast<uint64_t>( d_data );
}
template<typename TYPE, class Allocator>
typeID VectorDataDefault<TYPE, Allocator>::getType( size_t ) const
{
    constexpr auto type = getTypeID<TYPE>();
    return type;
}
template<typename TYPE, class Allocator>
size_t VectorDataDefault<TYPE, Allocator>::sizeofDataBlockType( size_t ) const
{
    return sizeof( TYPE );
}


/****************************************************************
 * Access the raw data blocks                                    *
 ****************************************************************/
template<typename TYPE, class Allocator>
inline void *VectorDataDefault<TYPE, Allocator>::getRawDataBlockAsVoid( size_t i )
{
    if ( i != 0 ) {
        return 0;
    }
    return d_data;
}
template<typename TYPE, class Allocator>
inline const void *VectorDataDefault<TYPE, Allocator>::getRawDataBlockAsVoid( size_t i ) const
{
    if ( i != 0 ) {
        return 0;
    }
    return d_data;
}


/****************************************************************
 * Access individual values                                      *
 ****************************************************************/
template<typename TYPE, class Allocator>
inline TYPE &VectorDataDefault<TYPE, Allocator>::operator[]( size_t i )
{
    return d_data[i];
}

template<typename TYPE, class Allocator>
inline const TYPE &VectorDataDefault<TYPE, Allocator>::operator[]( size_t i ) const
{
    return d_data[i];
}
template<typename TYPE, class Allocator>
inline void VectorDataDefault<TYPE, Allocator>::setValuesByLocalID( size_t num,
                                                                    const size_t *indices,
                                                                    const void *vals,
                                                                    const typeID &id )
{
    if ( id == getTypeID<TYPE>() ) {
        auto data = reinterpret_cast<const TYPE *>( vals );
        for ( size_t i = 0; i < num; i++ )
            d_data[indices[i]] = data[i];
    } else if ( id == getTypeID<double>() ) {
        auto data = reinterpret_cast<const double *>( vals );
        for ( size_t i = 0; i < num; ++i )
            d_data[indices[i]] = static_cast<TYPE>( data[i] );
    } else {
        AMP_ERROR( "Conversion not supported yet" );
    }
    if ( *d_UpdateState == UpdateState::UNCHANGED )
        *d_UpdateState = UpdateState::LOCAL_CHANGED;
}
template<typename TYPE, class Allocator>
inline void VectorDataDefault<TYPE, Allocator>::addValuesByLocalID( size_t num,
                                                                    const size_t *indices,
                                                                    const void *vals,
                                                                    const typeID &id )
{
    if ( id == getTypeID<TYPE>() ) {
        auto data = reinterpret_cast<const TYPE *>( vals );
        for ( size_t i = 0; i < num; i++ )
            d_data[indices[i]] += data[i];
    } else if ( id == getTypeID<double>() ) {
        auto data = reinterpret_cast<const double *>( vals );
        for ( size_t i = 0; i < num; ++i )
            d_data[indices[i]] += static_cast<TYPE>( data[i] );
    } else {
        AMP_ERROR( "Conversion not supported yet" );
    }
    if ( *d_UpdateState == UpdateState::UNCHANGED )
        *d_UpdateState = UpdateState::LOCAL_CHANGED;
}
template<typename TYPE, class Allocator>
inline void VectorDataDefault<TYPE, Allocator>::getValuesByLocalID( size_t num,
                                                                    const size_t *indices,
                                                                    void *vals,
                                                                    const typeID &id ) const
{
    if ( id == getTypeID<TYPE>() ) {
        auto data = reinterpret_cast<TYPE *>( vals );
        for ( size_t i = 0; i < num; i++ )
            data[i] = d_data[indices[i]];
    } else if ( id == getTypeID<double>() ) {
        auto data = reinterpret_cast<double *>( vals );
        for ( size_t i = 0; i < num; ++i )
            data[i] = d_data[indices[i]];
    } else {
        AMP_ERROR( "Conversion not supported yet" );
    }
}


/****************************************************************
 * Copy raw data                                                 *
 ****************************************************************/
template<typename TYPE, class Allocator>
void VectorDataDefault<TYPE, Allocator>::putRawData( const void *in, const typeID &id )
{
    if ( id == getTypeID<TYPE>() ) {
        memcpy( d_data, in, d_localSize * sizeof( TYPE ) );
    } else if ( id == getTypeID<double>() ) {
        auto data = reinterpret_cast<const double *>( in );
        for ( size_t i = 0; i < d_localSize; ++i )
            d_data[i] = static_cast<TYPE>( data[i] );
    } else {
        AMP_ERROR( "Conversion not supported yet" );
    }
}

template<typename TYPE, class Allocator>
void VectorDataDefault<TYPE, Allocator>::getRawData( void *out, const typeID &id ) const
{
    if ( id == getTypeID<TYPE>() ) {
        memcpy( out, d_data, d_localSize * sizeof( TYPE ) );
    } else if ( id == getTypeID<double>() ) {
        auto data = reinterpret_cast<double *>( out );
        for ( size_t i = 0; i < d_localSize; ++i )
            data[i] = static_cast<double>( d_data[i] );
    } else {
        AMP_ERROR( "Conversion not supported yet" );
    }
}


/****************************************************************
 * Swap raw data                                                 *
 ****************************************************************/
template<typename TYPE, class Allocator>
void VectorDataDefault<TYPE, Allocator>::swapData( VectorData &rhs )
{
    auto rhs2 = dynamic_cast<VectorDataDefault<TYPE, Allocator> *>( &rhs );
    AMP_INSIST( rhs2, "Cannot swap with arbitrary VectorData" );
    std::swap( d_CommList, rhs2->d_CommList );
    std::swap( d_UpdateState, rhs2->d_UpdateState );
    std::swap( d_Ghosts, rhs2->d_Ghosts );
    std::swap( d_AddBuffer, rhs2->d_AddBuffer );
    std::swap( d_data, rhs2->d_data );
    std::swap( d_localSize, rhs2->d_localSize );
    std::swap( d_globalSize, rhs2->d_globalSize );
    std::swap( d_localStart, rhs2->d_localStart );
}


/****************************************************************
 * Write/Read restart data                                       *
 ****************************************************************/
template<typename TYPE, class Allocator>
void VectorDataDefault<TYPE, Allocator>::registerChildObjects(
    AMP::IO::RestartManager *manager ) const
{
    VectorData::registerChildObjects( manager );
}
template<typename TYPE, class Allocator>
void VectorDataDefault<TYPE, Allocator>::writeRestart( int64_t fid ) const
{
    AMP::Array<TYPE> data( d_localSize );
    getRawData( data.data(), getTypeID<TYPE>() );
    writeHDF5( fid, "data", data );
    writeHDF5( fid, "localSize", d_localSize );
    writeHDF5( fid, "globalSize", d_globalSize );
    writeHDF5( fid, "localStart", d_localStart );
}
template<typename TYPE, class Allocator>
VectorDataDefault<TYPE, Allocator>::VectorDataDefault( int64_t fid, AMP::IO::RestartManager * )
{
    AMP::Array<TYPE> data;
    readHDF5( fid, "data", data );
    readHDF5( fid, "localSize", d_localSize );
    readHDF5( fid, "globalSize", d_globalSize );
    readHDF5( fid, "localStart", d_localStart );
    d_data = d_alloc.allocate( d_localSize );
    putRawData( data.data(), getTypeID<TYPE>() );
}


} // namespace AMP::LinearAlgebra

#endif
