#include "AMP/utils/AMP_MPI_pack.hpp"
#include "AMP/utils/UtilityMacros.h"

#include <complex>
#include <vector>


// Pack some common vectors
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
// clang-format off
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
// clang-format on


// std::string
template<>
size_t AMP::packSize<std::string>( const std::string &s )
{
    return s.size() + 1;
}
template<>
size_t AMP::pack<std::string>( const std::string &s, std::byte *buf )
{
    memcpy( buf, s.data(), s.size() + 1 );
    return s.size() + 1;
}
template<>
size_t AMP::unpack<std::string>( std::string &s, const std::byte *buf )
{
    s = std::string( reinterpret_cast<const char *>( buf ) );
    return s.size() + 1;
}


// std::vector<bool>
/*template<>
size_t AMP::packSize( const std::vector<bool>::reference &s )
{
    return s.size() + 1;
}
template<>
size_t AMP::pack( const std::vector<bool>::reference &s, std::byte *buf )
{
    memcpy( buf, s.data(), s.size() + 1 );
    return s.size() + 1;
}
template<>
size_t AMP::unpack( std::vector<bool>::reference &s, const std::byte *buf )
{
    s = std::vector<bool>::reference( reinterpret_cast<const char *>( buf ) );
    return s.size() + 1;
}*/
