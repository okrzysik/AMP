#ifndef included_AMP_MPI_pack
#define included_AMP_MPI_pack

#include "AMP/utils/UtilityMacros.h"

#include <cstring>
#include <tuple>
#include <vector>


namespace AMP {


/************************************************************************
 *  Functions for packing/unpacking data                                 *
 ************************************************************************/
template<class TYPE>
typename std::enable_if<!std::is_trivially_copyable<TYPE>::value, size_t>::type
packSize( const TYPE & );
template<class TYPE>
typename std::enable_if<!std::is_trivially_copyable<TYPE>::value, size_t>::type
pack( const TYPE &x, std::byte *buf );
template<class TYPE>
typename std::enable_if<!std::is_trivially_copyable<TYPE>::value, size_t>::type
unpack( TYPE &x, const std::byte *buf );
template<class TYPE>
typename std::enable_if<std::is_trivially_copyable<TYPE>::value, size_t>::type
packSize( const TYPE & )
{
    return sizeof( TYPE );
}
template<class TYPE>
typename std::enable_if<std::is_trivially_copyable<TYPE>::value, size_t>::type
pack( const TYPE &x, std::byte *buf )
{
    memcpy( buf, &x, sizeof( TYPE ) );
    return sizeof( TYPE );
}
template<class TYPE>
typename std::enable_if<std::is_trivially_copyable<TYPE>::value, size_t>::type
unpack( TYPE &x, const std::byte *buf )
{
    DISABLE_WARNINGS
    memcpy( &x, buf, sizeof( TYPE ) );
    ENABLE_WARNINGS
    return sizeof( TYPE );
}
template<class TYPE>
std::tuple<size_t, std::byte *> packArray( const TYPE *data, size_t length )
{
    size_t tot_bytes = 0;
    std::vector<size_t> bytes( length, 0 );
    for ( size_t i = 0; i < length; i++ ) {
        bytes[i] = packSize<TYPE>( data[i] );
        tot_bytes += bytes[i];
    }
    AMP_ASSERT( tot_bytes < static_cast<size_t>( std::numeric_limits<int>::max() ) );
    std::byte *buf = new std::byte[tot_bytes];
    memset( buf, 0, tot_bytes );
    size_t N = 0;
    for ( size_t i = 0; i < length; i++ ) {
        auto b = pack<TYPE>( data[i], &buf[N] );
        AMP_ASSERT( b == bytes[i] );
        N += b;
    }
    AMP_ASSERT( N == tot_bytes );
    return std::tuple<size_t, std::byte *>( tot_bytes, buf );
}
template<class TYPE>
size_t unpackArray( TYPE *data, size_t length, const std::byte *buf )
{
    size_t N = 0;
    for ( size_t i = 0; i < length; i++ ) {
        auto b = unpack<TYPE>( data[i], &buf[N] );
        N += b;
    }
    return N;
}

} // namespace AMP


#endif
