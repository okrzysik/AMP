#include "AMP/utils/AMP_MPI_pack.hpp"
#include "AMP/utils/Array.h"
#include "AMP/utils/Database.h"
#include "AMP/utils/Database.hpp"
#include "AMP/utils/MathExpr.h"
#include "AMP/utils/Utilities.h"

#include <cmath>
#include <complex>
#include <cstdlib>
#include <cstring>
#include <iomanip>
#include <limits>
#include <sstream>
#include <string>
#include <typeinfo>


namespace AMP {


/********************************************************************
 * KeyData factory functions                                         *
 ********************************************************************/
void registerMaterial( const std::string &name, std::function<std::unique_ptr<KeyData>()> fun )
{
    AMP::FactoryStrategy<KeyData>::registerFactory( name, fun );
}


/********************************************************************
 * Helper function to perform data conversion                        *
 ********************************************************************/
template<class TYPE>
static constexpr TYPE getTol()
{
    if constexpr ( std::is_same<TYPE, float>::value ) {
        return 1e-6;
    } else if constexpr ( std::is_same<TYPE, double>::value ) {
        return 1e-12;
    } else {
        return 0;
    }
}
template<class TYPE>
static constexpr TYPE abs( const TYPE &x )
{
    if constexpr ( std::is_signed<TYPE>::value ) {
        return x < 0 ? -x : x;
    } else {
        return x;
    }
}
template<class TYPE1, class TYPE2>
Array<TYPE2> convert( const Array<TYPE1> &x )
{
    if constexpr ( std::is_same<TYPE1, TYPE2>::value ) {
        return x;
    } else if constexpr ( std::is_arithmetic<TYPE1>::value && std::is_arithmetic<TYPE2>::value ) {
        Array<TYPE2> y( x.size() );
        y.fill( 0 );
        bool pass = true;
        auto tol  = getTol<TYPE1>();
        for ( size_t i = 0; i < x.length(); i++ ) {
            y( i ) = static_cast<TYPE2>( x( i ) );
            pass   = pass && abs( static_cast<TYPE1>( y( i ) - x( i ) ) ) <= tol * x( i );
        }
        if ( !pass ) {
            std::string type1 = typeid( TYPE1 ).name();
            std::string type2 = typeid( TYPE2 ).name();
            std::string msg = "Converting " + type1 + "-" + type2 + " results in loss of precision";
            AMP_WARNING( msg );
        }
        return y;
    } else {
        std::string type1 = typeid( TYPE1 ).name();
        std::string type2 = typeid( TYPE2 ).name();
        std::string msg =
            "Invalid conversion: " + type1 + "-" + type2 + " results in loss of precision";
        throw std::logic_error( msg );
    }
}


/********************************************************************
 * Scale data                                                        *
 ********************************************************************/
template<class TYPE>
void scaleData( TYPE &data, double factor )
{
    if constexpr ( std::is_same<TYPE, bool>::value ) {
        throw std::logic_error( "Unable to scale bool" );
    } else if constexpr ( std::is_same<TYPE, std::complex<float>>::value ) {
        data = static_cast<float>( factor ) * data;
    } else if constexpr ( std::is_same<TYPE, std::complex<double>>::value ) {
        data = factor * data;
    } else if constexpr ( std::is_arithmetic<TYPE>::value ) {
        data = static_cast<TYPE>( factor ) * data;
    } else {
        NULL_USE( factor );
        std::string type = typeid( TYPE ).name();
        throw std::logic_error( "Unable to scale " + type );
    }
}
template<class TYPE>
void scaleData( Array<TYPE> &data, double factor )
{
    if constexpr ( std::is_same<TYPE, bool>::value ) {
        throw std::logic_error( "Unable to scale bool" );
    } else if constexpr ( std::is_arithmetic<TYPE>::value ) {
        data.scale( factor );
    } else {
        NULL_USE( factor );
        std::string type = typeid( TYPE ).name();
        throw std::logic_error( "Unable to scale " + type );
    }
}


/********************************************************************
 * KeyDataScalar::operator==                                         *
 ********************************************************************/
template<class TYPE>
static inline bool compare( const TYPE &x, const TYPE &y );
template<>
inline bool compare( const double &x, const double &y )
{
    bool test = x == y;
    test      = test || fabs( x - y ) <= 1e-12 * fabs( x + y );
    test      = test || ( ( x != x ) && ( y != y ) );
    return test;
}
template<>
inline bool compare( const float &x, const float &y )
{
    bool test = x == y;
    test      = test || fabs( x - y ) <= 1e-7 * fabs( x + y );
    test      = test || ( ( x != x ) && ( y != y ) );
    return test;
}
template<>
inline bool compare( const std::complex<double> &x, const std::complex<double> &y )
{
    bool test = x == y;
    test      = test || std::abs( x - y ) <= 1e-12 * std::abs( x + y );
    test      = test || ( ( x != x ) && ( y != y ) );
    return test;
}
template<>
inline bool compare( const std::complex<float> &x, const std::complex<float> &y )
{
    bool test = x == y;
    test      = test || std::abs( x - y ) <= 1e-7 * std::abs( x + y );
    test      = test || ( ( x != x ) && ( y != y ) );
    return test;
}
template<class TYPE>
static inline bool compare( const TYPE &x, const TYPE &y )
{
    return x == y;
}
template<class TYPE>
bool KeyDataScalar<TYPE>::operator==( const KeyData &rhs ) const
{
    auto tmp1 = dynamic_cast<const KeyDataScalar<TYPE> *>( &rhs );
    auto tmp2 = dynamic_cast<const KeyDataArray<TYPE> *>( &rhs );
    if ( tmp1 ) {
        return compare( d_data, tmp1->d_data );
    } else if ( tmp2 ) {
        if ( tmp2->get().size() != 1 )
            return false;
        return compare( d_data, tmp2->get()( 0 ) );
    } else if ( ( is_floating_point() || is_integral() ) &&
                ( rhs.is_floating_point() || rhs.is_integral() ) ) {
        auto data1 = convertToDouble();
        auto data2 = rhs.convertToDouble();
        if ( data1.size() != data2.size() )
            return false;
        bool test = true;
        for ( size_t i = 0; i < data1.length(); i++ )
            test = test && compare( data1( i ), data2( i ) );
        return test;
    }
    return false;
}
template<class TYPE>
bool KeyDataArray<TYPE>::operator==( const KeyData &rhs ) const
{
    auto tmp1 = dynamic_cast<const KeyDataScalar<TYPE> *>( &rhs );
    auto tmp2 = dynamic_cast<const KeyDataArray<TYPE> *>( &rhs );
    if ( tmp1 ) {
        if ( d_data.size() != 1 )
            return false;
        return compare( d_data( 0 ), tmp1->get() );
    } else if ( tmp2 ) {
        return compare( d_data, tmp2->get() );
    } else if ( ( is_floating_point() || is_integral() ) &&
                ( rhs.is_floating_point() || rhs.is_integral() ) ) {
        auto data1 = convertToDouble();
        auto data2 = rhs.convertToDouble();
        if ( data1.size() != data2.size() )
            return false;
        bool test = true;
        for ( size_t i = 0; i < data1.length(); i++ )
            test = test && compare( data1( i ), data2( i ) );
        return test;
    }
    return false;
}


/********************************************************************
 * EquationKeyData                                                   *
 ********************************************************************/
EquationKeyData::EquationKeyData( std::string_view eq, const Units &unit ) : KeyData( unit )
{
    eq = deblank( eq );
    AMP_ASSERT( eq.size() > 4 );
    size_t i = eq.find( ')' );
    if ( eq[0] != '@' || eq[1] != '(' || eq.back() != ';' || i == std::string::npos )
        AMP_ERROR( "Equation appears invalid: " + std::string( eq ) );
    auto eq_var   = deblank( eq.substr( 2, i - 2 ) );
    auto equation = std::string( deblank( eq.substr( i + 1, eq.size() - i - 2 ) ) );
    std::vector<std::string> vars;
    if ( !eq_var.empty() ) {
        for ( size_t i = 0; i < eq_var.size(); ) {
            size_t j = std::min( eq_var.find( ',' ), eq_var.size() );
            vars.emplace_back( deblank( eq_var.substr( i, j - i ) ) );
            i = j + 1;
        }
    }
    d_eq = std::make_shared<MathExpr>( equation, vars );
}
EquationKeyData::EquationKeyData( std::shared_ptr<const MathExpr> eq, const Units &unit )
    : KeyData( unit ), d_eq( eq )
{
}
std::unique_ptr<KeyData> EquationKeyData::clone() const
{
    return std::make_unique<EquationKeyData>( d_eq, d_unit );
}
void EquationKeyData::print( std::ostream &os, std::string_view indent, bool, bool printType ) const
{
    os << indent << "@(";
    if ( d_eq ) {
        auto vars = d_eq->getVars();
        auto expr = d_eq->getExpr();
        if ( !vars.empty() ) {
            os << vars[0];
            for ( size_t i = 1; i < vars.size(); i++ )
                os << "," << vars[i];
        }
        os << ") " << expr << ";";
    }
    if ( printType )
        os << "  // " << getDataType().name;
    os << std::endl;
}
typeID EquationKeyData::getDataType() const { return AMP::getTypeID<MathExpr>(); }
bool EquationKeyData::is_floating_point() const { return d_eq->getVars().empty(); }
bool EquationKeyData::is_integral() const
{
    if ( d_eq->getVars().empty() ) {
        double v = ( *d_eq )();
        if ( static_cast<double>( static_cast<int64_t>( v ) ) == v )
            return true;
    }
    return false;
}
ArraySize EquationKeyData::arraySize() const
{
    return d_eq->getVars().empty() ? ArraySize() : ArraySize( 1 );
}
Array<double> EquationKeyData::convertToDouble() const
{
    if ( d_eq->getVars().empty() ) {
        double v = ( *d_eq )();
        return Array<double>( { 1 }, &v );
    }
    return Array<double>();
}
Array<int64_t> EquationKeyData::convertToInt64() const
{
    if ( d_eq->getVars().empty() ) {
        double v   = ( *d_eq )();
        int64_t v2 = v;
        if ( v2 == v )
            return Array<int64_t>( { 1 }, &v2 );
    }
    return Array<int64_t>();
}
bool EquationKeyData::operator==( const KeyData &rhs ) const
{
    auto rhs2 = dynamic_cast<const EquationKeyData *>( &rhs );
    if ( rhs2 == nullptr )
        return false;
    if ( d_eq == rhs2->d_eq )
        return true;
    if ( !d_eq || !rhs2->d_eq )
        return false;
    return *d_eq == *rhs2->d_eq;
}
size_t EquationKeyData::packSize() const
{
    std::string expr;
    std::vector<std::string> vars;
    if ( d_eq ) {
        expr = d_eq->getExpr();
        vars = d_eq->getVars();
    }
    return AMP::packSize( expr ) + AMP::packSize( vars );
}
size_t EquationKeyData::pack( std::byte *buf ) const
{
    std::string expr;
    std::vector<std::string> vars;
    if ( d_eq ) {
        expr = d_eq->getExpr();
        vars = d_eq->getVars();
    }
    size_t N = 0;
    N += AMP::pack( expr, &buf[N] );
    N += AMP::pack( vars, &buf[N] );
    return N;
}
size_t EquationKeyData::unpack( const std::byte *buf )
{
    std::string expr;
    std::vector<std::string> vars;
    size_t N = 0;
    N += AMP::unpack( expr, &buf[N] );
    N += AMP::unpack( vars, &buf[N] );
    d_eq.reset();
    if ( !expr.empty() )
        d_eq = std::make_shared<MathExpr>( expr, vars );
    return N;
}


/********************************************************************
 * DatabaseBox                                                       *
 ********************************************************************/
DatabaseBox::DatabaseBox() : d_dim( 0 )
{
    d_lower.fill( 0 );
    d_upper.fill( 0 );
}
DatabaseBox::DatabaseBox( int dim, const int *lower, const int *upper ) : d_dim( dim )
{
    AMP_ASSERT( dim <= 5 );
    d_lower.fill( 0 );
    d_upper.fill( 0 );
    for ( int d = 0; d < dim; d++ ) {
        d_lower[d] = lower[d];
        d_upper[d] = upper[d];
    }
}
DatabaseBox::DatabaseBox( std::string_view str ) : d_dim( 0 )
{
    d_lower.fill( 0 );
    d_upper.fill( 0 );
    size_t j     = 0;
    size_t k     = 0;
    char tmp[10] = { 0 };
    for ( size_t i = 0; i < str.length(); i++ ) {
        if ( str[i] == ' ' || str[i] == '(' || str[i] == '[' || str[i] == ']' ) {
            continue;
        } else if ( str[i] == ',' || str[i] == ')' ) {
            if ( str[i] == ')' && d_dim == 0 ) {
                d_dim = k + 1;
                while ( str[i] != ',' && i + 1 < str.length() ) {
                    ++i;
                }
            }
            tmp[j] = 0;
            j      = 0;
            int v  = atoi( tmp );
            if ( k < d_dim || d_dim == 0 )
                d_lower[k++] = v;
            else
                d_upper[( k++ ) - d_dim] = v;
        } else {
            tmp[j++] = str[i];
        }
    }
}
uint8_t &DatabaseBox::dim() { return d_dim; }
uint8_t DatabaseBox::dim() const { return d_dim; }
bool DatabaseBox::empty() const
{
    if ( d_dim == 0 )
        return true;
    for ( int d = 0; d < d_dim; d++ ) {
        if ( d_upper[d] < d_lower[d] )
            return true;
    }
    return false;
}
int &DatabaseBox::lower( uint8_t d )
{
    AMP_ASSERT( d_dim <= 5 && d < d_dim );
    return d_lower[d];
}
int DatabaseBox::lower( uint8_t d ) const
{
    AMP_ASSERT( d_dim <= 5 && d < d_dim );
    return d_lower[d];
}
int &DatabaseBox::upper( uint8_t d )
{
    AMP_ASSERT( d_dim <= 5 && d < d_dim );
    return d_upper[d];
}
int DatabaseBox::upper( uint8_t d ) const
{
    AMP_ASSERT( d_dim <= 5 && d < d_dim );
    return d_upper[d];
}
bool DatabaseBox::operator==( const DatabaseBox &box ) const
{
    bool equal = d_dim == box.d_dim;
    for ( int d = 0; d < d_dim; d++ ) {
        equal = equal && d_lower[d] == box.d_lower[d];
        equal = equal && d_upper[d] == box.d_upper[d];
    }
    return equal;
}
std::ostream &operator<<( std::ostream &out, const DatabaseBox &box )
{
    if ( box.empty() )
        out << "[(),()]";
    out << "[(" << box.lower( 0 );
    for ( int d = 1; d < box.dim(); d++ )
        out << "," << box.lower( d );
    out << "),(" << box.upper( 0 );
    for ( int d = 1; d < box.dim(); d++ )
        out << "," << box.upper( d );
    out << ")";
    return out;
}


/********************************************************************
 * Explicit instantiations                                           *
 ********************************************************************/
#define instantiate( FUN )       \
    FUN( bool );                 \
    FUN( char );                 \
    FUN( int8_t );               \
    FUN( int16_t );              \
    FUN( int32_t );              \
    FUN( int64_t );              \
    FUN( uint8_t );              \
    FUN( uint16_t );             \
    FUN( uint32_t );             \
    FUN( uint64_t );             \
    FUN( float );                \
    FUN( double );               \
    FUN( long double );          \
    FUN( std::complex<float> );  \
    FUN( std::complex<double> ); \
    FUN( std::string );          \
    FUN( DatabaseBox )
#define instantiateConvert( TYPE )                                       \
    template Array<TYPE> convert( const Array<bool> & );                 \
    template Array<TYPE> convert( const Array<char> & );                 \
    template Array<TYPE> convert( const Array<int8_t> & );               \
    template Array<TYPE> convert( const Array<int16_t> & );              \
    template Array<TYPE> convert( const Array<int32_t> & );              \
    template Array<TYPE> convert( const Array<int64_t> & );              \
    template Array<TYPE> convert( const Array<uint8_t> & );              \
    template Array<TYPE> convert( const Array<uint16_t> & );             \
    template Array<TYPE> convert( const Array<uint32_t> & );             \
    template Array<TYPE> convert( const Array<uint64_t> & );             \
    template Array<TYPE> convert( const Array<float> & );                \
    template Array<TYPE> convert( const Array<double> & );               \
    template Array<TYPE> convert( const Array<long double> & );          \
    template Array<TYPE> convert( const Array<std::complex<float>> & );  \
    template Array<TYPE> convert( const Array<std::complex<double>> & ); \
    template Array<TYPE> convert( const Array<std::string> & );          \
    template Array<TYPE> convert( const Array<DatabaseBox> & )
#define instantiateScaleData( TYPE )                             \
    template void scaleData<TYPE>( TYPE & data, double factor ); \
    template void scaleData<TYPE>( Array<TYPE> & data, double factor )
#define instantiateKeyDataScalar( TYPE ) template class KeyDataScalar<TYPE>
#define instantiateKeyDataArray( TYPE ) template class KeyDataArray<TYPE>
#define instantiateIsType( TYPE ) \
    template bool Database::isType<TYPE>( std::string_view, source_location ) const
#define instantiatePutScalar( TYPE ) \
    template void Database::putScalar<TYPE>( std::string_view, TYPE, Units, Database::Check )
#define instantiateGetScalar( TYPE )                                                            \
    template TYPE Database::getScalar<TYPE>( std::string_view, const Units &, source_location ) \
        const
#define instantiateGetArray( TYPE )                                                             \
    template Array<TYPE> Database::getArray( std::string_view, const Units &, source_location ) \
        const
#define instantiateGetVector( TYPE )                \
    template std::vector<TYPE> Database::getVector( \
        std::string_view, const Units &, source_location ) const
#define instantiatePutVector( TYPE )         \
    template void Database::putVector<TYPE>( \
        std::string_view, const std::vector<TYPE> &, Units, Database::Check )
#define instantiatePutArray( TYPE ) \
    template void Database::putArray<TYPE>( std::string_view, Array<TYPE>, Units, Database::Check )
#define instantiateGetWithDefault( TYPE )                                                          \
    template TYPE Database::getWithDefault<TYPE>(                                                  \
        std::string_view, Database::IdentityType<TYPE const &>::type, const Units & ) const;       \
    template std::vector<TYPE> Database::getWithDefault<std::vector<TYPE>>(                        \
        std::string_view, Database::IdentityType<std::vector<TYPE> const &>::type, const Units & ) \
        const;                                                                                     \
    template Array<TYPE> Database::getWithDefault<Array<TYPE>>(                                    \
        std::string_view, Database::IdentityType<Array<TYPE> const &>::type, const Units & ) const

instantiate( instantiateConvert );        // convert
instantiate( instantiateScaleData );      // scaleData
instantiate( instantiateKeyDataScalar );  // KeyDataScalar
instantiate( instantiateKeyDataArray );   // KeyDataArray
instantiate( instantiateIsType );         // Database::isType
instantiate( instantiateGetScalar );      // Database::getScalar
instantiate( instantiateGetVector );      // Database::getVector
instantiate( instantiateGetArray );       // Database::getArray
instantiate( instantiatePutScalar );      // Database::putScalar
instantiate( instantiatePutVector );      // Database::putVector
instantiate( instantiatePutArray );       // Database::putArray
instantiate( instantiateGetWithDefault ); // Database::getWithDefault
template void
Database::putScalar<const char *>( std::string_view, const char *, Units, Database::Check );
template void
Database::putScalar<std::_Bit_reference>( std::string_view key, std::_Bit_reference, Units, Check );


/********************************************************
 *  Register KeyData with factory                        *
 ********************************************************/
#define registerkeyData2( TYPE, TYPENAME )                             \
    REGISTER_KEYDATA( KeyDataScalar<TYPE>, KeyDataScalar_##TYPENAME ); \
    REGISTER_KEYDATA( KeyDataArray<TYPE>, KeyDataArray##TYPENAME )
#define registerkeyData( TYPE ) registerkeyData2( TYPE, TYPE )
registerkeyData( bool );
registerkeyData( char );
registerkeyData( int8_t );
registerkeyData( int16_t );
registerkeyData( int32_t );
registerkeyData( int64_t );
registerkeyData( uint8_t );
registerkeyData( uint16_t );
registerkeyData( uint32_t );
registerkeyData( uint64_t );
registerkeyData( float );
registerkeyData( double );
registerkeyData2( long double, long_double );
registerkeyData2( std::complex<float>, complex_float );
registerkeyData2( std::complex<double>, complex_double );
registerkeyData2( std::string, string );
registerkeyData( DatabaseBox );
REGISTER_KEYDATA( EmptyKeyData, EmptyKeyData );
REGISTER_KEYDATA( DatabaseVector, DatabaseVector );
REGISTER_KEYDATA( EquationKeyData, EquationKeyData );

} // namespace AMP


/********************************************************
 *  Explicit instantiations of Array<DatabaseBox>        *
 ********************************************************/
#include "AMP/utils/Array.hpp"
instantiateArrayConstructors( AMP::DatabaseBox );
PACK_UNPACK_ARRAY( AMP::DatabaseBox )
