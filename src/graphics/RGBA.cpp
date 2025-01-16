#include "AMP/graphics/RGBA.h"
#include "AMP/utils/Array.hpp"


namespace AMP {


/********************************************************
 *  RGBA32 / ARGB32                                      *
 ********************************************************/
static_assert( sizeof( ARGB32 ) == 4 );
static_assert( sizeof( RGBA32 ) == 4 );
static constexpr bool runTests()
{
    constexpr ARGB32 argb( (uint32_t) 0x01020304 );
    constexpr RGBA32 rgba( (uint32_t) 0x02030401 );
    constexpr ARGB32 argb2 = argb;
    static_assert( argb.red() == 2 );
    static_assert( argb.green() == 3 );
    static_assert( argb.blue() == 4 );
    static_assert( argb.alpha() == 1 );
    static_assert( argb.red() == rgba.red() );
    static_assert( argb.green() == rgba.green() );
    static_assert( argb.blue() == rgba.blue() );
    static_assert( argb.alpha() == rgba.alpha() );
    static_assert( argb2 == argb );
    return true;
}
static_assert( runTests() );


} // namespace AMP


/********************************************************
 *  Explicit instantiations of Array<RGB>                *
 ********************************************************/
instantiateArrayConstructors( AMP::RGBA32 );
instantiateArrayConstructors( AMP::ARGB32 );
typedef AMP::Array<AMP::RGBA32> RGBA_Array;
typedef AMP::Array<AMP::ARGB32> ARGB_Array;
template RGBA_Array RGBA_Array::repmat( const std::vector<size_t> & ) const;
template ARGB_Array ARGB_Array::repmat( const std::vector<size_t> & ) const;
template RGBA_Array RGBA_Array::subset( const std::vector<size_t> & ) const;
template ARGB_Array ARGB_Array::subset( const std::vector<size_t> & ) const;
template void RGBA_Array::copySubset( const std::vector<size_t> &, const RGBA_Array & );
template void ARGB_Array::copySubset( const std::vector<size_t> &, const ARGB_Array & );
