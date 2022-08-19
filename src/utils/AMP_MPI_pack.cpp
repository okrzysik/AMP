#include "AMP/utils/AMP_MPI_pack.hpp"
#include "AMP/utils/UtilityMacros.h"

#include <complex>
#include <vector>


#define PACK_VECTOR( TYPE )                                              \
    template<>                                                           \
    size_t AMP::packSize( const std::vector<TYPE> &x )                   \
    {                                                                    \
        return sizeof( size_t ) + sizeof( TYPE ) * x.size();             \
    }                                                                    \
    template<>                                                           \
    size_t AMP::pack( const std::vector<TYPE> &x, std::byte *buf )       \
    {                                                                    \
        pack<size_t>( x.size(), buf );                                   \
        buf += sizeof( size_t );                                         \
        if ( !x.empty() )                                                \
            memcpy( buf, (void *) x.data(), sizeof( TYPE ) * x.size() ); \
        return sizeof( size_t ) + sizeof( TYPE ) * x.size();             \
    }                                                                    \
    template<>                                                           \
    size_t AMP::unpack( std::vector<TYPE> &x, const std::byte *buf )     \
    {                                                                    \
        size_t N;                                                        \
        unpack<size_t>( N, buf );                                        \
        x.resize( N );                                                   \
        buf += sizeof( size_t );                                         \
        if ( !x.empty() )                                                \
            memcpy( (void *) x.data(), buf, sizeof( TYPE ) * x.size() ); \
        return sizeof( size_t ) + sizeof( TYPE ) * x.size();             \
    }

//clang-format off
PACK_VECTOR( char )
PACK_VECTOR( int8_t )
PACK_VECTOR( int16_t )
PACK_VECTOR( int32_t )
PACK_VECTOR( int64_t )
PACK_VECTOR( uint8_t )
PACK_VECTOR( uint16_t )
PACK_VECTOR( uint32_t )
PACK_VECTOR( uint64_t )
PACK_VECTOR( float )
PACK_VECTOR( double )
PACK_VECTOR( std::complex<float> )
PACK_VECTOR( std::complex<double> )
// clang-format onn


// std::vector<string>
template<>
size_t AMP::packSize( const std::vector<std::string> &x )
{
    size_t N = sizeof( size_t );
    for ( size_t i = 0; i < x.size(); i++ )
        N += packSize( x[i] );
    return N;
}
template<>
size_t AMP::pack( const std::vector<std::string> &x, std::byte *buf )
{
    size_t N = pack<size_t>( x.size(), buf );
    for ( size_t i = 0; i < x.size(); i++ )
        N += pack( x[i], &buf[N] );
    return N;
}
template<>
size_t AMP::unpack( std::vector<std::string> &x, const std::byte *buf )
{
    size_t length;
    size_t N = unpack<size_t>( length, buf );
    x.resize( length );
    for ( size_t i = 0; i < x.size(); i++ )
        N += unpack( x[i], &buf[N] );
    return N;
}
