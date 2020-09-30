#ifndef included_AMP_Scalar_hpp
#define included_AMP_Scalar_hpp


#include "AMP/utils/UtilityMacros.h"

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
inline std::tuple<T1, double> Scalar::convert( const std::any &x0 )
{
    T2 x     = std::any_cast<T2>( x0 );
    T1 y     = static_cast<T1>( x );
    double e = std::abs<double>( x - static_cast<T2>( std::real( y ) ) );
    return std::tie( y, e );
}
template<class TYPE>
inline size_t Scalar::get_hash()
{
    static auto hash = typeid( TYPE ).hash_code();
    return hash;
}
template<class TYPE>
inline void Scalar::store( const TYPE &x )
{
    d_hash = get_hash<TYPE>();
    d_data = std::make_any<TYPE>( x );
}


/********************************************************************
 * Contructor                                                        *
 ********************************************************************/
inline Scalar::Scalar() : d_type( 0 ), d_hash( 0 ) {}
template<class TYPE>
Scalar::Scalar( TYPE x ) : d_type( get_type<TYPE>() ), d_hash( 0 )
{
    // Store the scalar
    if constexpr ( std::is_integral<TYPE>::value ) {
        if constexpr ( std::is_signed<TYPE>::value ) {
            if ( x >= std::numeric_limits<int64_t>::min() &&
                 x <= std::numeric_limits<int64_t>::max() ) {
                store( static_cast<int64_t>( x ) );
            } else {
                store( x );
            }
        } else {
            if ( x <= std::numeric_limits<int64_t>::max() ) {
                store( static_cast<int64_t>( x ) );
            } else {
                store( x );
            }
        }
    } else if constexpr ( std::is_same<TYPE, float>::value ) {
        store( static_cast<double>( x ) );
    } else if constexpr ( std::is_same<TYPE, double>::value ) {
        store( x );
    } else if constexpr ( std::is_same<TYPE, long double>::value ) {
        store( static_cast<long double>( x ) );
    } else if constexpr ( AMP::is_complex<TYPE>::value ) {
        store( std::complex<double>( x.real(), x.imag() ) );
    } else {
        store( x );
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
inline TYPE Scalar::get( double tol ) const
{
    // Special cases for performance
    if ( d_hash == get_hash<TYPE>() )
        return std::any_cast<TYPE>( d_data );
    // Check that the data exists
    if ( !d_data.has_value() )
        AMP_ERROR( "Scalar has no data" );
    // Convert the data
    TYPE y   = 0;
    double e = 0;
    if ( d_hash == get_hash<double>() ) {
        std::tie( y, e ) = convert<TYPE, double>( d_data );
    } else if ( d_hash == get_hash<long double>() ) {
        std::tie( y, e ) = convert<TYPE, long double>( d_data );
    } else if ( d_hash == get_hash<int64_t>() ) {
        std::tie( y, e ) = convert<TYPE, int64_t>( d_data );
    } else if ( d_hash == get_hash<std::complex<double>>() ) {
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
