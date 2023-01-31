#ifndef included_AMP_Scalar_hpp
#define included_AMP_Scalar_hpp


#include "AMP/utils/AMP_MPI.h"
#include "AMP/utils/UtilityMacros.h"

#include <complex>
#include <limits>
#include <math.h>


// is_complex
// clang-format off
namespace AMP {
template<class T> struct is_complex : public std::false_type {};
template<class T> struct is_complex<const T> : public is_complex<T> {};
template<class T> struct is_complex<volatile const T> : public is_complex<T> {};
template<class T> struct is_complex<volatile T> : public is_complex<T> {};
template<class T> struct is_complex<std::complex<T>> : public std::true_type {};
template<class T> inline constexpr bool is_complex_v = is_complex<T>::value;
}
// clang-format on


namespace AMP {


/********************************************************************
 * Helper functions                                                  *
 ********************************************************************/
#if defined( USING_ICC )
DISABLE_WARNINGS
#endif
template<class TYPE>
constexpr double Scalar::getTol()
{
    if constexpr ( std::is_integral_v<TYPE> ) {
        return 0;
    } else if constexpr ( std::is_floating_point_v<TYPE> ) {
        constexpr double tol = 10 * std::numeric_limits<TYPE>::epsilon();
        return tol;
    } else if constexpr ( AMP::is_complex_v<TYPE> ) {
        constexpr double tol = std::numeric_limits<TYPE>::epsilon().real();
        return tol;
    } else {
        constexpr double tol = 10 * std::numeric_limits<TYPE>::epsilon();
        return tol;
    }
}
template<class TYPE>
inline void Scalar::store( const TYPE &x )
{
    constexpr auto hash = getTypeID<TYPE>().hash;
    d_hash              = hash;
    d_data              = std::make_any<TYPE>( x );
    if constexpr ( std::is_integral_v<TYPE> ) {
        d_type = 'i';
    } else if constexpr ( std::is_floating_point_v<TYPE> ) {
        d_type = 'f';
    } else if constexpr ( AMP::is_complex_v<TYPE> ) {
        d_type = 'c';
    } else {
        d_type = 0;
    }
}
#if defined( USING_ICC )
ENABLE_WARNINGS
#endif


/********************************************************************
 * Contructor                                                        *
 ********************************************************************/
inline Scalar::Scalar() : d_type( 0 ), d_hash( 0 ) {}
template<class TYPE>
Scalar::Scalar( const TYPE &x )
{
    // Special case if we are dealing with a Scalar
    if constexpr ( std::is_same_v<TYPE, Scalar> ) {
        d_type = x.d_type;
        d_hash = x.d_hash;
        d_data = x.d_data;
    } else {
        // Store the data
        if constexpr ( std::is_integral_v<TYPE> ) {
            if constexpr ( std::is_signed_v<TYPE> ) {
                if ( x >= std::numeric_limits<int64_t>::min() &&
                     x <= std::numeric_limits<int64_t>::max() ) {
                    store<int64_t>( x );
                } else {
                    store<TYPE>( x );
                }
            } else {
                if ( x <= static_cast<uint64_t>( std::numeric_limits<int64_t>::max() ) ) {
                    store<int64_t>( x );
                } else {
                    store<TYPE>( x );
                }
            }
        } else if constexpr ( std::is_same_v<TYPE, float> ) {
            store<double>( x );
        } else if constexpr ( std::is_same_v<TYPE, double> ) {
            store<double>( x );
        } else if constexpr ( std::is_same_v<TYPE, long double> ) {
            store<long double>( x );
        } else if constexpr ( AMP::is_complex_v<TYPE> ) {
            store( std::complex<double>( x.real(), x.imag() ) );
        } else {
            store<TYPE>( x );
        }
        // Check that we can get the data back
#if ( defined( DEBUG ) || defined( _DEBUG ) ) && !defined( NDEBUG )
        if constexpr ( std::is_integral_v<TYPE> ) {
            auto z = get<TYPE>();
            AMP_ASSERT( x == z );
        } else if constexpr ( std::is_floating_point_v<TYPE> ) {
            constexpr TYPE inf = std::numeric_limits<TYPE>::infinity();
            constexpr TYPE tol = 10 * std::abs( std::numeric_limits<TYPE>::epsilon() );
            if ( x != inf && x != -inf && x == x ) {
                auto z = get<TYPE>();
                AMP_ASSERT( std::abs( x - z ) <= tol * std::abs( x ) );
            }
        }
#endif
    }
}


/********************************************************************
 * Get                                                               *
 ********************************************************************/
template<class TYPE>
inline TYPE Scalar::get( double tol ) const
{
    // Special cases for performance
    constexpr auto typeHash = getTypeID<TYPE>().hash;
    if ( d_hash == typeHash )
        return std::any_cast<TYPE>( d_data );
    // Convert the data
    TYPE y;
    if constexpr ( AMP::is_complex_v<TYPE> ) {
        // Return complex data
        constexpr auto complexHash = getTypeID<std::complex<double>>().hash;
        if ( d_hash == complexHash ) {
            auto x = std::any_cast<std::complex<double>>( d_data );
            return TYPE( x.real(), x.imag() );
        } else {
            double x = get<double>( tol );
            return TYPE( x );
        }
    } else {
        double e;
        constexpr auto doubleHash  = getTypeID<double>().hash;
        constexpr auto longHash    = getTypeID<long double>().hash;
        constexpr auto int64Hash   = getTypeID<int64_t>().hash;
        constexpr auto complexHash = getTypeID<std::complex<double>>().hash;
        if ( d_hash == doubleHash ) {
            auto x = std::any_cast<double>( d_data );
            y      = static_cast<TYPE>( x );
            e      = std::abs<double>( x - static_cast<double>( y ) );
        } else if ( d_hash == longHash ) {
            auto x = std::any_cast<long double>( d_data );
            y      = static_cast<TYPE>( x );
            e      = std::abs<double>( x - static_cast<long double>( y ) );
        } else if ( d_hash == int64Hash ) {
            auto x = std::any_cast<int64_t>( d_data );
            y      = static_cast<TYPE>( x );
            e      = std::abs<double>( x - static_cast<int64_t>( y ) );
        } else if ( d_hash == complexHash ) {
            auto x = std::any_cast<std::complex<double>>( d_data );
            y      = x.real();
            e      = std::abs<double>( x - static_cast<std::complex<double>>( y ) );
        } else if ( !d_data.has_value() ) {
            // Data does not exist
            AMP_ERROR( "Scalar has no data" );
        } else {
            std::string t1 = d_data.type().name();
            std::string t2 = typeid( TYPE ).name();
            AMP_ERROR( "Unable to convert data: " + t1 + " -> " + t2 );
        }
        if ( e > tol )
            AMP_ERROR( "Error exceeds tolerance converting data" );
    }
    return y;
}


/********************************************************************
 * Operator overloading                                              *
 ********************************************************************/
Scalar operator-( const Scalar &x );
Scalar operator+( const Scalar &x, const Scalar &y );
Scalar operator-( const Scalar &x, const Scalar &y );
Scalar operator*( const Scalar &x, const Scalar &y );
Scalar operator/( const Scalar &x, const Scalar &y );
inline bool operator==( double x, const Scalar &y ) { return y.operator==( x ); }


/********************************************************************
 * Special functions                                                 *
 ********************************************************************/
Scalar minReduce( const AMP::AMP_MPI &comm, const Scalar &x );
Scalar maxReduce( const AMP::AMP_MPI &comm, const Scalar &x );
Scalar sumReduce( const AMP::AMP_MPI &comm, const Scalar &x );


/********************************************************
 *  ostream operator                                     *
 ********************************************************/
template<class TYPE>
typename std::enable_if_t<std::is_same_v<TYPE, Scalar>, std::ostream &>
operator<<( std::ostream &out, const TYPE &x );


} // namespace AMP

#endif
