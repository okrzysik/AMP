#include "AMP/utils/Array.h"
#include "AMP/utils/Array.hpp"
#include "AMP/utils/UtilityMacros.h"

#include <complex>


namespace AMP {


/********************************************************
 *  Explicit instantiations of Array                     *
 ********************************************************/
#if defined( USING_ICC )
DISABLE_WARNINGS
#endif
template class Array<char>;
template class Array<uint8_t>;
template class Array<uint16_t>;
template class Array<uint32_t>;
template class Array<uint64_t>;
template class Array<int8_t>;
template class Array<int16_t>;
template class Array<int32_t>;
template class Array<int64_t>;
template class Array<float>;
template class Array<double>;
template class Array<long double>;
instantiateArrayConstructors( bool );
instantiateArrayConstructors( std::string );
instantiateArrayConstructors( std::string_view );


/********************************************************
 *  Explicit instantiations of Array<std::complex>       *
 ********************************************************/
instantiateArrayConstructors( std::complex<float> );
instantiateArrayConstructors( std::complex<double> );
template Array<std::complex<float>>::Array( const Range<std::complex<float>> &range );
template Array<std::complex<double>>::Array( const Range<std::complex<double>> &range );
template Array<std::complex<float>>
Array<std::complex<float>>::repmat( std::vector<unsigned long> const & ) const;
template Array<std::complex<double>>
Array<std::complex<double>>::repmat( std::vector<unsigned long> const & ) const;
template void Array<std::complex<float>>::copySubset( std::vector<unsigned long> const &,
                                                      Array<std::complex<float>> const & );
template void Array<std::complex<double>>::copySubset( std::vector<unsigned long> const &,
                                                       Array<std::complex<double>> const & );
template Array<std::complex<float>>
Array<std::complex<float>>::subset( std::vector<unsigned long> const & ) const;
template Array<std::complex<double>>
Array<std::complex<double>>::subset( std::vector<unsigned long> const & ) const;
template bool Array<std::complex<float>>::NaNs() const;
template bool Array<std::complex<double>>::NaNs() const;
template void AMP::Array<std::complex<double>, AMP::FunctionTable>::rand();


#if defined( USING_ICC )
ENABLE_WARNINGS
#endif

} // namespace AMP
