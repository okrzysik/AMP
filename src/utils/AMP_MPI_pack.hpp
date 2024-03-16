#ifndef included_AMP_MPI_pack
#define included_AMP_MPI_pack

#include "AMP/utils/Array.h"
#include "AMP/utils/TypeTraits.h"
#include "AMP/utils/UtilityMacros.h"

#include <cstddef>
#include <cstring>
#include <limits>
#include <tuple>
#include <vector>


namespace AMP {


/************************************************************************
 *  Default functions for packing/unpacking data                         *
 ************************************************************************/
// clang-format off
template<class TYPE>
typename std::enable_if_t<!std::is_trivially_copyable_v<TYPE> &&
                        !AMP::is_shared_ptr_v<TYPE> &&
                        !AMP::is_Array_v<TYPE> &&
                        !AMP::is_vector_v<TYPE>,
                        size_t>
packSize( const TYPE & );
template<class TYPE>
typename std::enable_if_t<!std::is_trivially_copyable_v<TYPE> &&
                        !AMP::is_shared_ptr_v<TYPE> &&
                        !AMP::is_Array_v<TYPE> &&
                        !AMP::is_vector_v<TYPE>,
                        size_t>
pack( const TYPE &x, std::byte *buf );
template<class TYPE>
typename std::enable_if_t<!std::is_trivially_copyable_v<TYPE> &&
                        !AMP::is_shared_ptr_v<TYPE> &&
                        !AMP::is_Array_v<TYPE> &&
                        !AMP::is_vector_v<TYPE>,
                        size_t>
unpack( TYPE &x, const std::byte *buf );
template<class TYPE>
typename std::enable_if_t<AMP::is_shared_ptr_v<TYPE>, size_t> packSize( const TYPE &x );
template<class TYPE>
typename std::enable_if_t<AMP::is_shared_ptr_v<TYPE>, size_t> pack( const TYPE &x,
                                                                             std::byte *buf );
template<class TYPE>
typename std::enable_if_t<AMP::is_shared_ptr_v<TYPE>, size_t>
unpack( TYPE &x, const std::byte *buf );
// clang-format on


/************************************************************************
 *  Functions for packing/unpacking trivial data                         *
 ************************************************************************/
template<class TYPE>
typename std::enable_if_t<std::is_trivially_copyable_v<TYPE>, size_t> packSize( const TYPE & )
{
    return sizeof( TYPE );
}
template<class TYPE>
typename std::enable_if_t<std::is_trivially_copyable_v<TYPE>, size_t> pack( const TYPE &x,
                                                                            std::byte *buf )
{
    memcpy( buf, &x, sizeof( TYPE ) );
    return sizeof( TYPE );
}
template<class TYPE>
typename std::enable_if_t<std::is_trivially_copyable_v<TYPE>, size_t> unpack( TYPE &x,
                                                                              const std::byte *buf )
{
    DISABLE_WARNINGS
    memcpy( &x, buf, sizeof( TYPE ) );
    ENABLE_WARNINGS
    return sizeof( TYPE );
}


/************************************************************************
 *  Functions for packing/unpacking AMP::Array                           *
 ************************************************************************/
template<class TYPE>
typename std::enable_if_t<AMP::is_Array_v<TYPE>, size_t> packSize( const TYPE &x )
{
    return x.packSize();
}
template<class TYPE>
typename std::enable_if_t<AMP::is_Array_v<TYPE>, size_t> pack( const TYPE &x, std::byte *buf )
{
    return x.pack( buf );
}
template<class TYPE>
typename std::enable_if_t<AMP::is_Array_v<TYPE>, size_t> unpack( TYPE &x, const std::byte *buf )
{
    return x.unpack( buf );
}


/************************************************************************
 *  Functions for packing/unpacking std::vector                          *
 ************************************************************************/
template<class TYPE>
typename std::enable_if_t<AMP::is_vector_v<TYPE>, size_t> packSize( const TYPE &x )
{
    size_t N = sizeof( size_t );
    for ( size_t i = 0; i < x.size(); i++ )
        N += packSize( x[i] );
    return N;
}
template<class TYPE>
typename std::enable_if_t<AMP::is_vector_v<TYPE>, size_t> pack( const TYPE &x, std::byte *buf )
{
    size_t N = pack<size_t>( x.size(), buf );
    for ( size_t i = 0; i < x.size(); i++ )
        N += pack( x[i], &buf[N] );
    return N;
}
template<class TYPE>
typename std::enable_if_t<AMP::is_vector_v<TYPE>, size_t> unpack( TYPE &x, const std::byte *buf )
{
    size_t length;
    size_t N = unpack<size_t>( length, buf );
    x.resize( length );
    for ( size_t i = 0; i < x.size(); i++ )
        N += unpack( x[i], &buf[N] );
    return N;
}


/************************************************************************
 *  Functions for packing/unpacking shared_pointers                      *
 ************************************************************************/
template<class TYPE>
typename std::enable_if_t<AMP::is_shared_ptr_v<TYPE>, size_t> packSize( const TYPE &x )
{
    return packSize( *x );
}
template<class TYPE>
typename std::enable_if_t<AMP::is_shared_ptr_v<TYPE>, size_t> pack( const TYPE &x, std::byte *buf )
{
    return pack( *x, buf );
}
template<class TYPE>
typename std::enable_if_t<AMP::is_shared_ptr_v<TYPE>, size_t> unpack( TYPE &x,
                                                                      const std::byte *buf )
{
    return unpack( *x, buf );
}


/************************************************************************
 *  Functions for packing/unpacking arrays of data                       *
 ************************************************************************/
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
