#include "AMP/vectors/Scalar.h"


/********************************************************************
 * boolean operator overloading                                      *
 ********************************************************************/
bool AMP::Scalar::operator<( const AMP::Scalar &rhs ) const
{
    if ( has_value() || rhs.has_value() )
        AMP_ERROR( "Comparing empty scalar" );
    if ( is_floating_point() ) {
        return get<double>() < rhs.get<double>();
    } else if ( is_integral() ) {
        return get<int64_t>() < rhs.get<int64_t>();
    } else if ( is_complex() ) {
        AMP_ERROR( "No operator < for std::complex" );
    } else {
        AMP_ERROR( "Unable to get types for Scalar" );
    }
    return false;
}
bool AMP::Scalar::operator==( const AMP::Scalar &rhs ) const
{
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
    return false;
}
bool AMP::Scalar::operator>( const AMP::Scalar &rhs ) const { return rhs < *this; }
bool AMP::Scalar::operator<=( const AMP::Scalar &rhs ) const { return !( *this > rhs ); }
bool AMP::Scalar::operator>=( const AMP::Scalar &rhs ) const { return !( *this < rhs ); }
bool AMP::Scalar::operator!=( const AMP::Scalar &rhs ) const { return !( *this == rhs ); }


/********************************************************************
 * arithmetic operator overloading                                   *
 ********************************************************************/
AMP::Scalar operator+( const AMP::Scalar &x, const AMP::Scalar &y )
{
    if ( !x.has_value() )
        return y;
    if ( !y.has_value() )
        return x;
    if ( x.is_floating_point() ) {
        return x.get<double>() + y.get<double>();
    } else if ( x.is_integral() ) {
        return x.get<int64_t>() + y.get<int64_t>();
    } else if ( x.is_complex() ) {
        return x.get<std::complex<double>>() + y.get<std::complex<double>>();
    } else {
        AMP_ERROR( "Unable to get types for Scalar" );
    }
    return AMP::Scalar();
}
AMP::Scalar operator-( const AMP::Scalar &x, const AMP::Scalar &y )
{
    if ( !x.has_value() )
        return y;
    if ( !y.has_value() )
        return x;
    if ( x.is_floating_point() ) {
        return x.get<double>() - y.get<double>();
    } else if ( x.is_integral() ) {
        return x.get<int64_t>() - y.get<int64_t>();
    } else if ( x.is_complex() ) {
        return x.get<std::complex<double>>() - y.get<std::complex<double>>();
    } else {
        AMP_ERROR( "Unable to get types for Scalar" );
    }
    return AMP::Scalar();
}
AMP::Scalar operator*( const AMP::Scalar &x, const AMP::Scalar &y )
{
    if ( !x.has_value() || !y.has_value() )
        return 0.0;
    if ( x.is_floating_point() ) {
        return x.get<double>() * y.get<double>();
    } else if ( x.is_integral() ) {
        return x.get<int64_t>() * y.get<int64_t>();
    } else if ( x.is_complex() ) {
        return x.get<std::complex<double>>() * y.get<std::complex<double>>();
    } else {
        AMP_ERROR( "Unable to get types for Scalar" );
    }
    return AMP::Scalar();
}


/********************************************************************
 * Special functions                                                 *
 ********************************************************************/
AMP::Scalar minReduce( const AMP::AMP_MPI &comm, const AMP::Scalar &x )
{
    if ( comm.getSize() <= 1 )
        return x;
    if ( x.is_floating_point() ) {
        return comm.minReduce( x.get<double>() );
    } else if ( x.is_integral() ) {
        return comm.minReduce( x.get<int64_t>() );
    } else if ( x.is_complex() ) {
        return comm.minReduce( x.get<std::complex<double>>() );
    } else {
        AMP_ERROR( "Unable to get types for Scalar" );
    }
    return AMP::Scalar();
}
AMP::Scalar maxReduce( const AMP::AMP_MPI &comm, const AMP::Scalar &x )
{
    if ( comm.getSize() <= 1 )
        return x;
    if ( x.is_floating_point() ) {
        return comm.maxReduce( x.get<double>() );
    } else if ( x.is_integral() ) {
        return comm.maxReduce( x.get<int64_t>() );
    } else if ( x.is_complex() ) {
        return comm.maxReduce( x.get<std::complex<double>>() );
    } else {
        AMP_ERROR( "Unable to get types for Scalar" );
    }
    return AMP::Scalar();
}
AMP::Scalar sumReduce( const AMP::AMP_MPI &comm, const AMP::Scalar &x )
{
    if ( comm.getSize() <= 1 )
        return x;
    if ( x.is_floating_point() ) {
        return comm.sumReduce( x.get<double>() );
    } else if ( x.is_integral() ) {
        return comm.sumReduce( x.get<int64_t>() );
    } else if ( x.is_complex() ) {
        return comm.sumReduce( x.get<std::complex<double>>() );
    } else {
        AMP_ERROR( "Unable to get types for Scalar" );
    }
    return AMP::Scalar();
}
AMP::Scalar sqrt( const AMP::Scalar &x )
{
    if ( !x.has_value() )
        return 0.0;
    if ( x.is_floating_point() ) {
        return sqrt( x.get<double>() );
    } else if ( x.is_integral() ) {
        return sqrt( x.get<int64_t>() );
    } else if ( x.is_complex() ) {
        return sqrt( x.get<std::complex<double>>() );
    } else {
        AMP_ERROR( "Unable to get types for Scalar" );
    }
    return AMP::Scalar();
}
