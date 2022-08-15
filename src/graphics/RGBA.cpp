#include "AMP/graphics/RGBA.h"
#include "AMP/utils/Array.hpp"


namespace AMP {


/********************************************************
 *  Explicit instantiations of Array<RGB>                *
 ********************************************************/
instantiateArrayConstructors( RGBA );
template Array<RGBA>::Array( const Range<RGBA> & );
template Array<std::complex<double>> &
Array<std::complex<double>>::operator=( const std::vector<std::complex<double>> & );
template Array<RGBA> Array<RGBA>::repmat( const std::vector<size_t> & ) const;
template void Array<RGBA>::copySubset( const std::vector<size_t> &, const Array<RGBA> & );
template Array<RGBA> Array<RGBA>::subset( const std::vector<size_t> & ) const;


} // namespace AMP
