// clang-format off
#include "AMP/utils/Array.h"
#include "AMP/utils/Array.hpp"
#if defined( USING_ICC )
#include "AMP/utils/UtilityMacros.h"
#endif

#include <complex>


namespace AMP {


/********************************************************
 *  Explicit instantiations of Array                     *
 ********************************************************/
#if defined( USING_ICC )
DISABLE_WARNINGS
#endif
template class Array<char,FunctionTable>;
template class Array<uint8_t,FunctionTable>;
template class Array<uint16_t,FunctionTable>;
template class Array<uint32_t,FunctionTable>;
template class Array<uint64_t,FunctionTable>;
template class Array<int8_t,FunctionTable>;
template class Array<int16_t,FunctionTable>;
template class Array<int32_t,FunctionTable>;
template class Array<int64_t,FunctionTable>;
template class Array<float,FunctionTable>;
template class Array<double,FunctionTable>;


/********************************************************
 *  Explicit instantiations of Array<bool>               *
 ********************************************************/
template Array<bool,FunctionTable>::Array();
template Array<bool,FunctionTable>::~Array();
template Array<bool,FunctionTable>::Array( size_t );
template Array<bool,FunctionTable>::Array( size_t, size_t );
template Array<bool,FunctionTable>::Array( size_t, size_t, size_t );
template Array<bool,FunctionTable>::Array( size_t, size_t, size_t, size_t );
template Array<bool,FunctionTable>::Array( size_t, size_t, size_t, size_t, size_t );
template Array<bool,FunctionTable>::Array( const std::vector<size_t>&, const bool* );
template Array<bool,FunctionTable>::Array( std::string );
template Array<bool,FunctionTable>::Array( std::initializer_list<bool> );
template Array<bool,FunctionTable>::Array( const Array<bool,FunctionTable>& );
template Array<bool,FunctionTable>::Array( Array<bool,FunctionTable>&& );
template Array<bool,FunctionTable>& Array<bool,FunctionTable>::operator=( const Array<bool,FunctionTable>& );
template Array<bool,FunctionTable>& Array<bool,FunctionTable>::operator=( Array<bool,FunctionTable>&& );
template Array<bool,FunctionTable>& Array<bool,FunctionTable>::operator=( const std::vector<bool>& );
template void Array<bool,FunctionTable>::clear();
template bool Array<bool,FunctionTable>::operator==(Array<bool,FunctionTable> const&) const;
template void Array<bool,FunctionTable>::resize( ArraySize const& );


/********************************************************
 *  Explicit instantiations of Array<std::complex>       *
 ********************************************************/
instantiateArrayConstructors( std::complex<float> )
instantiateArrayConstructors( std::complex<double> )
template Array<std::complex<float>,FunctionTable>::Array( const Range<std::complex<float>>& range );
template Array<std::complex<double>,FunctionTable>::Array( const Range<std::complex<double>>& range );
template void Array<std::complex<float>,FunctionTable>::resize( ArraySize const& );
template void Array<std::complex<double>,FunctionTable>::resize( ArraySize const& );
template void Array<std::complex<float>,FunctionTable>::clear();
template void Array<std::complex<double>,FunctionTable>::clear();
template bool Array<std::complex<float>,FunctionTable>::operator==(Array<std::complex<float>,FunctionTable> const&) const;
template bool Array<std::complex<double>,FunctionTable>::operator==(Array<std::complex<double>,FunctionTable> const&) const;
template Array<std::complex<float>,FunctionTable> Array<std::complex<float>,FunctionTable>::repmat(std::vector<unsigned long> const&) const;
template Array<std::complex<double>,FunctionTable> Array<std::complex<double>,FunctionTable>::repmat(std::vector<unsigned long> const&) const;
template void Array<std::complex<float>,FunctionTable>::copySubset(std::vector<unsigned long> const&, Array<std::complex<float>,FunctionTable> const&);
template void Array<std::complex<double>,FunctionTable>::copySubset(std::vector<unsigned long> const&, Array<std::complex<double>,FunctionTable> const&);
template Array<std::complex<float>,FunctionTable> Array<std::complex<float>,FunctionTable>::subset(std::vector<unsigned long> const&) const;
template Array<std::complex<double>,FunctionTable> Array<std::complex<double>,FunctionTable>::subset(std::vector<unsigned long> const&) const;
template bool Array<std::complex<float>,FunctionTable>::NaNs() const;
template bool Array<std::complex<double>,FunctionTable>::NaNs() const;


/********************************************************
 *  Explicit instantiations of Array<std::string>        *
 ********************************************************/
instantiateArrayConstructors( std::string )
template void Array<std::string,FunctionTable>::resize( ArraySize const& );

#if defined( USING_ICC )
ENABLE_WARNINGS
#endif

} // AMP namespace
