#ifndef included_AMP_Scalar_hpp
#define included_AMP_Scalar_hpp


#include "AMP/utils/UtilityMacros.h"
#include "AMP/vectors/Scalar.h"

#include <complex>
#include <limits>


// is_complex
// clang-format off
namespace AMP {
template<class T> struct is_complex : public std::false_type {};
template<class T> struct is_complex<const T> : public is_complex<T> {};
template<class T> struct is_complex<volatile const T> : public is_complex<T> {};
template<class T> struct is_complex<volatile T> : public is_complex<T> {};
template<class T> struct is_complex<std::complex<T>> : public std::true_type {};
}
// clang-format on


namespace AMP {


/********************************************************************
 * Helper functions                                                  *
 ********************************************************************/
template<class TYPE>
constexpr char Scalar::get_type()
{
    if constexpr ( std::is_integral<TYPE>::value ) {
        return 'i';
    } else if constexpr ( std::is_floating_point<TYPE>::value ) {
        return 'f';
    } else if constexpr ( AMP::is_complex<TYPE>::value ) {
        return 'c';
    } else {
        return 0;
    }
}
template<class TYPE>
constexpr double Scalar::getTol()
{
    if constexpr ( std::is_integral<TYPE>::value ) {
        return 0;
    } else if constexpr ( std::is_floating_point<TYPE>::value ) {
        return 10 * std::abs( std::numeric_limits<TYPE>::epsilon() );
    } else {
        return 10 * std::abs( std::numeric_limits<TYPE>::epsilon() );
    }
}
template<class T1, class T2>
std::tuple<T1, double> Scalar::convert( const std::any &x0 )
{
    T2 x     = std::any_cast<T2>( x0 );
    T1 y     = static_cast<T1>( x );
    double e = std::abs<double>( x - static_cast<T2>( std::real( y ) ) );
    return std::tie( y, e );
}


/********************************************************************
 * Contructor                                                        *
 ********************************************************************/
template<class TYPE>
Scalar::Scalar( TYPE x ) : d_type( get_type<TYPE>() ), d_data( nullptr )
{
    // Store the scalar
    if constexpr ( std::is_integral<TYPE>::value ) {
        if constexpr ( std::is_signed<TYPE>::value ) {
            if ( x >= std::numeric_limits<int64_t>::min() &&
                 x <= std::numeric_limits<int64_t>::max() ) {
                d_data = std::any( static_cast<int64_t>( x ) );
            } else {
                d_data = std::any( x );
            }
        } else {
            if ( x <= std::numeric_limits<int64_t>::max() ) {
                d_data = std::any( static_cast<int64_t>( x ) );
            } else {
                d_data = std::any( x );
            }
        }
    } else if constexpr ( std::is_same<TYPE, float>::value ) {
        d_data = std::any( static_cast<double>( x ) );
    } else if constexpr ( std::is_same<TYPE, double>::value ) {
        d_data = std::any( x );
    } else if constexpr ( std::is_same<TYPE, long double>::value ) {
        d_data = std::any( static_cast<long double>( x ) );
    } else if constexpr ( AMP::is_complex<TYPE>::value ) {
        d_data = std::any( std::complex<double>( x.real(), x.imag() ) );
    } else {
        d_data = std::any( x );
    }
    // Check that we can get the data back
#if ( defined( DEBUG ) || defined( _DEBUG ) ) && !defined( NDEBUG )
    auto y = get<TYPE>();
    if constexpr ( std::is_integral<TYPE>::value ) {
        AMP_ASSERT( x == y );
    } else {
        auto tol = 10 * std::abs( std::numeric_limits<TYPE>::epsilon() );
        AMP_ASSERT( std::abs( x - y ) <= tol * std::abs( x ) );
    }
#endif
}


/********************************************************************
 * Get                                                               *
 ********************************************************************/
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
    } else if ( type == typeid( double ).hash_code() ) {
        std::tie( y, e ) = convert<TYPE, double>( d_data );
    } else if ( type == typeid( long double ).hash_code() ) {
        std::tie( y, e ) = convert<TYPE, long double>( d_data );
    } else if ( type == typeid( int64_t ).hash_code() ) {
        std::tie( y, e ) = convert<TYPE, int64_t>( d_data );
    } else if ( type == typeid( std::complex<double> ).hash_code() ) {
        auto x = std::any_cast<std::complex<double>>( d_data );
        if constexpr ( AMP::is_complex<TYPE>::value ) {
            return TYPE( x.real(), x.imag() );
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
