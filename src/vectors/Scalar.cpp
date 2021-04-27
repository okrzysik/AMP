#include "AMP/vectors/Scalar.h"

#include <math.h>


namespace AMP {


/********************************************************************
 * boolean operator overloading                                      *
 ********************************************************************/
bool Scalar::operator<( const Scalar &rhs ) const
{
    if ( !has_value() || !rhs.has_value() )
        AMP_ERROR( "Comparing empty scalar" );
    if ( is_complex() || rhs.is_complex() ) {
        AMP_ERROR( "No operator < for std::complex" );
    } else if ( is_floating_point() || rhs.is_floating_point() ) {
        return get<double>() < rhs.get<double>();
    } else if ( is_integral() && rhs.is_integral() ) {
        return get<int64_t>() < rhs.get<int64_t>();
    } else {
        AMP_ERROR( "Unable to get types for Scalar" );
    }
    return false;
}
bool Scalar::operator==( const Scalar &rhs ) const
{
    try {
        if ( has_value() != rhs.has_value() )
            return false;
        if ( is_floating_point() ) {
            return get<double>() == rhs.get<double>();
        } else if ( is_integral() ) {
            return get<int64_t>() == rhs.get<int64_t>();
        } else if ( is_complex() ) {
            return get<std::complex<double>>() == rhs.get<std::complex<double>>();
        } else {
            AMP_ERROR( "Unable to get types for Scalar" );
        }
    } catch ( ... ) {
    }
    return false;
}
bool Scalar::operator>( const Scalar &rhs ) const { return rhs < *this; }
bool Scalar::operator<=( const Scalar &rhs ) const { return !( *this > rhs ); }
bool Scalar::operator>=( const Scalar &rhs ) const { return !( *this < rhs ); }
bool Scalar::operator!=( const Scalar &rhs ) const { return !( *this == rhs ); }


/********************************************************************
 * arithmetic operator overloading                                   *
 ********************************************************************/
Scalar operator+( const Scalar &x, const Scalar &y )
{
    if ( !x.has_value() )
        return y;
    if ( !y.has_value() )
        return x;
    if ( x.is_floating_point() ) {
        return Scalar::create( x.get<double>() + y.get<double>() );
    } else if ( x.is_integral() ) {
        return Scalar::create( x.get<int64_t>() + y.get<int64_t>() );
    } else if ( x.is_complex() ) {
        return Scalar::create( x.get<std::complex<double>>() + y.get<std::complex<double>>() );
    } else {
        AMP_ERROR( "Unable to get types for Scalar" );
    }
    return Scalar();
}
Scalar operator-( const Scalar &x, const Scalar &y )
{
    if ( !x.has_value() )
        return y;
    if ( !y.has_value() )
        return x;
    if ( x.is_floating_point() ) {
        return Scalar::create( x.get<double>() - y.get<double>() );
    } else if ( x.is_integral() ) {
        return Scalar::create( x.get<int64_t>() - y.get<int64_t>() );
    } else if ( x.is_complex() ) {
        return Scalar::create( x.get<std::complex<double>>() - y.get<std::complex<double>>() );
    } else {
        AMP_ERROR( "Unable to get types for Scalar" );
    }
    return Scalar();
}
Scalar operator*( const Scalar &x, const Scalar &y )
{
    if ( !x.has_value() || !y.has_value() )
        return 0.0;
    if ( x.is_floating_point() ) {
        return Scalar::create( x.get<double>() * y.get<double>() );
    } else if ( x.is_integral() ) {
        return Scalar::create( x.get<int64_t>() * y.get<int64_t>() );
    } else if ( x.is_complex() ) {
        return Scalar::create( x.get<std::complex<double>>() * y.get<std::complex<double>>() );
    } else {
        AMP_ERROR( "Unable to get types for Scalar" );
    }
    return Scalar();
}


/********************************************************************
 * Special functions                                                 *
 ********************************************************************/
Scalar minReduce( const AMP::AMP_MPI &comm, const Scalar &x )
{
    if ( comm.getSize() <= 1 )
        return x;
    if ( x.is_floating_point() ) {
        return Scalar::create( comm.minReduce( x.get<double>() ) );
    } else if ( x.is_integral() ) {
        return Scalar::create( comm.minReduce( x.get<int64_t>() ) );
    } else if ( x.is_complex() ) {
        return Scalar::create( comm.minReduce( x.get<std::complex<double>>() ) );
    } else {
        AMP_ERROR( "Unable to get types for Scalar" );
    }
    return Scalar();
}
Scalar maxReduce( const AMP::AMP_MPI &comm, const Scalar &x )
{
    if ( comm.getSize() <= 1 )
        return x;
    if ( x.is_floating_point() ) {
        return Scalar::create( comm.maxReduce( x.get<double>() ) );
    } else if ( x.is_integral() ) {
        return Scalar::create( comm.maxReduce( x.get<int64_t>() ) );
    } else if ( x.is_complex() ) {
        return Scalar::create( comm.maxReduce( x.get<std::complex<double>>() ) );
    } else {
        AMP_ERROR( "Unable to get types for Scalar" );
    }
    return Scalar();
}
Scalar sumReduce( const AMP::AMP_MPI &comm, const Scalar &x )
{
    if ( comm.getSize() <= 1 )
        return x;
    if ( x.is_floating_point() ) {
        return Scalar::create( comm.sumReduce( x.get<double>() ) );
    } else if ( x.is_integral() ) {
        return Scalar::create( comm.sumReduce( x.get<int64_t>() ) );
    } else if ( x.is_complex() ) {
        return Scalar::create( comm.sumReduce( x.get<std::complex<double>>() ) );
    } else {
        AMP_ERROR( "Unable to get types for Scalar" );
    }
    return Scalar();
}
Scalar Scalar::sqrt() const
{
    if ( !has_value() )
        return 0.0;
    if ( is_floating_point() ) {
        return Scalar::create( ::sqrt( get<double>() ) );
    } else if ( is_integral() ) {
        return Scalar::create( ::sqrt( get<int64_t>() ) );
    } else if ( is_complex() ) {
        return Scalar::create( ::sqrt( get<std::complex<double>>() ) );
    } else {
        AMP_ERROR( "Unable to get types for Scalar" );
    }
    return Scalar();
}
Scalar Scalar::abs() const
{
    if ( !has_value() )
        return Scalar();
    if ( is_floating_point() ) {
        return Scalar::create( std::abs( get<double>() ) );
    } else if ( is_integral() ) {
        return Scalar::create( std::abs( get<int64_t>() ) );
    } else if ( is_complex() ) {
        return Scalar::create( std::abs( get<std::complex<double>>() ) );
    } else {
        AMP_ERROR( "Unable to get types for Scalar" );
    }
    return Scalar();
}

/********************************************************
 *  ostream operator                                     *
 ********************************************************/
template<>
std::enable_if<std::is_same<AMP::Scalar, AMP::Scalar>::value, std::ostream &>::type
operator<<<AMP::Scalar>( std::ostream &out, const AMP::Scalar &x )
{
    if ( !x.has_value() )
        return out;
    if ( x.is_floating_point() ) {
        out << x.get<double>();
    } else if ( x.is_integral() ) {
        out << x.get<int64_t>();
    } else if ( x.is_complex() ) {
        out << x.get<std::complex<double>>();
    } else {
        AMP_ERROR( "Unable to get types for Scalar" );
    }
    return out;
}

} // namespace AMP
