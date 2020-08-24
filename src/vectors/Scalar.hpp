#ifndef included_AMP_Scalar_hpp
#define included_AMP_Scalar_hpp


#include "AMP/utils/Utilities.h"
#include "AMP/vectors/Scalar.h"

#include <complex>
#include <limits>


namespace AMP {


template<class TYPE>
inline char get_type()
{
    if constexpr ( std::is_integral<TYPE>::value ) {
        return 'i';
    } else if constexpr ( std::is_floating_point<TYPE>::value ) {
        return 'f';
    } else if constexpr ( std::is_same<TYPE, std::complex<float>>::value ||
                          std::is_same<TYPE, std::complex<double>>::value ) {
        return 'c';
    } else {
        return 0;
    }
}


/********************************************************************
 * Contructor                                                        *
 ********************************************************************/
template<class TYPE>
Scalar::Scalar( TYPE x ) : d_type( get_type<TYPE>() ), d_data( nullptr )
{
    // Store the scalar
    if constexpr ( std::is_integral<TYPE>::value ) {
        if ( std::numeric_limits<int64_t>::max() && x < std::numeric_limits<int64_t>::max() ) {
            d_data = std::any( static_cast<int64_t>( x ) );
        } else {
            std::any( x );
        }
    } else if constexpr ( std::is_same<TYPE, float>::value ) {
        d_data = std::any( static_cast<double>( x ) );
    } else if constexpr ( std::is_same<TYPE, double>::value ) {
        d_data = std::any( x );
    } else if constexpr ( std::is_same<TYPE, long double>::value ) {
        d_data = std::any( static_cast<long double>( x ) );
    } else if constexpr ( std::is_same<TYPE, std::complex<float>>::value ) {
        d_data = std::any( std::complex<double>( x.real(), x.imag() ) );
    } else if constexpr ( std::is_same<TYPE, std::complex<double>>::value ) {
        d_data = std::any( x );
    } else {
        std::any( x );
    }
    // Check that we can get the data back
    auto y = get<TYPE>();
    if constexpr ( std::is_integral<TYPE>::value ) {
        AMP_ASSERT( x == y );
    } else {
        AMP_ASSERT( std::abs( x - y ) <= 1e-14 * std::abs( x ) );
    }
}


/********************************************************************
 * Get                                                               *
 ********************************************************************/
template<class T1, class T2>
std::tuple<T1, double> convert( const std::any &x0 )
{
    T2 x     = std::any_cast<T2>( x0 );
    T1 y     = static_cast<T1>( x );
    double e = std::abs<double>( x - static_cast<T2>( std::real( y ) ) );
    return std::tie( y, e );
}
template<class TYPE>
TYPE Scalar::get( double tol ) const
{
    if ( !d_data.has_value() )
        AMP_ERROR( "Scalar has no data" );
    auto type = d_data.type().hash_code();
    TYPE y;
    double e;
    if ( type == typeid( TYPE ).hash_code() ) {
        return std::any_cast<TYPE>( d_data );
    } else if ( type == typeid( int64_t ).hash_code() ) {
        std::tie( y, e ) = convert<TYPE, int64_t>( d_data );
    } else if ( type == typeid( double ).hash_code() ) {
        std::tie( y, e ) = convert<TYPE, double>( d_data );
    } else if ( type == typeid( long double ).hash_code() ) {
        std::tie( y, e ) = convert<TYPE, long double>( d_data );
    } else if ( type == typeid( std::complex<double> ).hash_code() ) {
        auto x = std::any_cast<std::complex<double>>( d_data );
        if constexpr ( std::is_same<TYPE, std::complex<float>>::value ) {
            return std::complex<float>( x.real(), x.imag() );
        } else if constexpr ( std::is_same<TYPE, std::complex<double>>::value ) {
            return x;
        } else {
            y = x.real();
            e = std::abs<double>( x - static_cast<std::complex<double>>( std::real( y ) ) );
        }
    } else {
        std::string t1 = d_data.type().name();
        std::string t2 = typeid( TYPE ).name();
        AMP_ERROR( "Unable to convert data: " + t1 + " -> " + t2 );
    }
    if ( e > tol )
        AMP_ERROR( "Error exceeds tolerance converting data" );
    return y;
}


} // namespace AMP

#endif
