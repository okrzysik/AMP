#include "AMP/materials/Property.h"
#include "AMP/utils/Array.hpp"
#include "AMP/utils/Utilities.h"
#include "AMP/vectors/MultiVector.h"
#include "AMP/vectors/Vector.h"

#include <algorithm>


namespace AMP::Materials {


/************************************************************************
 *  Constructor                                                          *
 ************************************************************************/
Property::Property( std::string name,
                    const ArraySize &size,
                    const Units &unit,
                    std::string source,
                    std::vector<std::string> args,
                    std::vector<std::array<double, 2>> ranges,
                    std::vector<Units> units )
    : d_name( std::move( name ) ),
      d_dim( size ),
      d_units( unit ),
      d_source( std::move( source ) ),
      d_arguments( std::move( args ) ),
      d_argUnits( std::move( units ) ),
      d_defaults( {} ),
      d_ranges( std::move( ranges ) )
{
    if ( d_argUnits.empty() )
        d_argUnits.resize( d_arguments.size() );
    AMP_ASSERT( d_arguments.size() == d_ranges.size() );
    AMP_ASSERT( d_arguments.size() == d_argUnits.size() );
    for ( size_t i = 0; i < d_arguments.size(); i++ )
        d_argToIndexMap[d_arguments[i]] = i;
    d_defaults.resize( d_arguments.size() );
    for ( size_t i = 0; i < d_arguments.size(); i++ )
        d_defaults[i] = d_ranges[i][0];
}


/************************************************************************
 *  make_map                                                             *
 ************************************************************************/
static std::string getArgName( const std::string &vecname,
                               const std::map<std::string, std::string> &translator )
{
    if ( translator.empty() )
        return vecname;
    for ( auto &elem : translator ) {
        if ( elem.second == vecname )
            return elem.first;
    }
    return "";
}
std::map<std::string, std::shared_ptr<AMP::LinearAlgebra::Vector>>
Property::make_map( const std::shared_ptr<AMP::LinearAlgebra::MultiVector> &args,
                    const std::map<std::string, std::string> &translator ) const
{
    std::map<std::string, std::shared_ptr<AMP::LinearAlgebra::Vector>> result;
    if ( !d_arguments.empty() ) {
        size_t xls = translator.size();
        AMP_INSIST( xls > 0, "attempt to make MultiVector map without setting translator" );
        for ( auto vec : *args ) {
            auto key = getArgName( vec->getName(), translator );
            auto it  = std::find( d_arguments.begin(), d_arguments.end(), key );
            if ( it != d_arguments.end() )
                result.insert( std::make_pair( key, vec ) );
        }
    }
    return result;
}


/************************************************************************
 *  eval                                                                *
 ************************************************************************/
double Property::eval( const Units &unit,
                       const std::vector<double> &args,
                       const std::vector<std::string> &names,
                       const std::vector<Units> &argUnits ) const
{
    double value = 0;
    auto result  = AMP::Array<double>::staticView( { 1 }, &value );
    // Get the arguments
    AMP::Array<double> args2;
    if ( names.empty() && argUnits.empty() && args.size() == d_arguments.size() ) {
        args2 =
            AMP::Array<double>::staticView( { args.size() }, const_cast<double *>( args.data() ) );
    } else {
        // Set the default values
        args2.resize( d_defaults.size() );
        for ( size_t i = 0; i < d_defaults.size(); i++ )
            args2( i ) = d_defaults[i];
        // Override default arguments
        AMP_ASSERT( names.size() == args.size() );
        for ( size_t i = 0; i < args.size(); i++ ) {
            for ( size_t j = 0; j < d_arguments.size(); j++ ) {
                if ( names[i] == d_arguments[j] ) {
                    args2( i ) = args[i];
                    if ( !argUnits.empty() ) {
                        if ( !argUnits[i].isNull() && !d_argUnits[j].isNull() )
                            args2( i ) *= argUnits[i].convert( d_argUnits[j] );
                    }
                }
            }
        }
    }
    // Evaluate the property
    eval( result, args2 );
    // Convert units if required
    if ( !unit.isNull() )
        value *= d_units.convert( unit );
    return value;
}


/************************************************************************
 *  evalv                                                                *
 ************************************************************************/
static inline size_t getSize( const std::vector<double> &x ) { return x.size(); }
static inline size_t getSize( const AMP::LinearAlgebra::Vector &x ) { return x.getLocalSize(); }
template<class VEC>
AMP::Array<double> Property::loadArgs( size_t N,
                                       const std::vector<argumentDataStruct<VEC>> &args ) const
{
    AMP::Array<double> args2( d_arguments.size(), N );
    // Load default values
    for ( size_t j = 0; j < N; j++ )
        for ( size_t i = 0; i < d_arguments.size(); i++ )
            args2( i, j ) = d_defaults[i];
    // Load the arguments
    for ( const auto &arg : args ) {
        AMP_ASSERT( getSize( arg.vec ) == N );
        int i = get_arg_index( std::string( arg.str ) );
        if ( i == -1 )
            continue;
        double scale = 1.0;
        if ( !arg.units.isNull() && !d_argUnits[i].isNull() )
            scale = arg.units.convert( d_argUnits[i] );
        auto range = d_ranges[i];
        auto it    = arg.vec.begin();
        for ( size_t j = 0; j < N; j++, ++it ) {
            double value  = *it * scale;
            args2( i, j ) = value;
            bool test     = value >= range[0] && value <= range[1];
            if ( !test ) {
                std::stringstream ss;
                ss << "Property '" + d_arguments[i] + "' out of range in function '" + d_name +
                          "'\n";
                ss << "Value is " << value << " ";
                ss << std::endl
                   << "Valid range is [" << range[0] << "," << range[1] << "]" << std::endl;
                AMP_ERROR( ss.str() );
            }
        }
    }
    return args2;
}
template<class OUT, class IN>
void Property::evalv( OUT &r,
                      const Units &unit,
                      const std::vector<argumentDataStruct<IN>> &args ) const
{
    AMP_ASSERT( d_dim.length() == 1 );
    size_t N = getSize( r );

    // Load the arguments into an array
    auto args2 = loadArgs( N, args );

    // Allocate temporary output array
    Array<double> r2( N );
    r2.fill( 0 );

    // Call eval
    eval( r2, args2 );

    // Convert units
    if ( !unit.isNull() )
        r2.scale( d_units.convert( unit ) );

    // Copy the results back
    auto it = r.begin();
    for ( size_t i = 0; i < N; i++, ++it )
        *it = r2( i );
}
template void Property::evalv( std::vector<double> &,
                               const Units &unit,
                               const std::vector<argumentDataStruct<std::vector<double>>> & ) const;
void Property::evalv(
    std::shared_ptr<AMP::LinearAlgebra::Vector> &r,
    const std::map<std::string, std::shared_ptr<AMP::LinearAlgebra::Vector>> &args ) const
{
    std::vector<argumentDataStruct<AMP::LinearAlgebra::Vector>> args2;
    for ( const auto &t : args )
        args2.emplace_back( std::get<0>( t ), *std::get<1>( t ), std::get<1>( t )->getUnits() );
    evalv( *r, r->getUnits(), args2 );
}
void Property::evalv( std::shared_ptr<AMP::LinearAlgebra::Vector> &r,
                      const std::shared_ptr<AMP::LinearAlgebra::MultiVector> &args,
                      const std::map<std::string, std::string> &translator ) const
{
    auto mapargs = make_map( args, translator );
    evalv( r, mapargs );
}
void Property::evalv(
    std::vector<double> &r,
    const std::map<std::string, std::shared_ptr<std::vector<double>>> &args ) const
{
    std::vector<argumentDataStruct<std::vector<double>>> args2;
    for ( const auto &tmp : args )
        args2.emplace_back( tmp.first, *tmp.second );
    evalv( r, Units(), args2 );
}


/************************************************************************
 *  Misc functions                                                       *
 ************************************************************************/
std::array<double, 2> Property::get_arg_range( const std::string &name ) const
{
    int index = get_arg_index( name );
    AMP_ASSERT( index >= 0 );
    return d_ranges[index];
}
bool Property::is_argument( const std::string &name ) const { return get_arg_index( name ) >= 0; }
void Property::checkArgs( const AMP::Array<double> &args ) const
{
    AMP_ASSERT( args.ndim() <= 2 );
    AMP_ASSERT( args.size( 0 ) == d_dim.length() );
    for ( size_t i = 0; i < args.size( 0 ); i++ ) {
        auto range = d_ranges[i];
        for ( size_t j = 0; j < args.size( 1 ); j++ ) {
            auto value = args( i, j );
            bool test  = value >= range[0] && value <= range[1];
            if ( !test ) {
                std::stringstream ss;
                ss << "Property '" + d_arguments[i] + "' out of range in function '" + d_name +
                          "'\n";
                ss << "Value is " << value << " ";
                ss << std::endl
                   << "Valid range is [" << range[0] << "," << range[1] << "]" << std::endl;
                AMP_ERROR( ss.str() );
            }
        }
    }
}


/************************************************************************
 *  evalv (vector)                                                       *
 ************************************************************************/
template<class OUT, class IN>
void Property::evalv( std::vector<std::shared_ptr<OUT>> &r,
                      const Units &unit,
                      const std::vector<argumentDataStruct<IN>> &args ) const
{
    AMP_ASSERT( d_dim == ArraySize( r.size() ) );
    size_t N = getSize( *r[0] );

    // Load the arguments into an array
    auto args2 = loadArgs( N, args );

    // Allocate temporary output array
    Array<double> r2( r.size(), N );
    r2.fill( 0 );

    // Call eval
    eval( r2, args2 );

    // Convert units
    if ( !unit.isNull() )
        r2.scale( d_units.convert( unit ) );

    // Copy the results back
    size_t N0 = r.size();
    for ( size_t i = 0; i < N0; i++ ) {
        AMP_ASSERT( getSize( *r[i] ) == N );
        auto it = r[i]->begin();
        for ( size_t j = 0; j < N; j++, ++it )
            *it = r2( i + j * N0 );
    }
}
void Property::evalv(
    std::vector<std::shared_ptr<AMP::LinearAlgebra::Vector>> &r,
    const std::map<std::string, std::shared_ptr<AMP::LinearAlgebra::Vector>> &args ) const
{
    std::vector<argumentDataStruct<AMP::LinearAlgebra::Vector>> args2;
    for ( const auto &tmp : args )
        args2.emplace_back( tmp.first, *tmp.second );
    evalv( r, {}, args2 );
}
void Property::evalv( std::vector<std::shared_ptr<AMP::LinearAlgebra::Vector>> &r,
                      const std::shared_ptr<AMP::LinearAlgebra::MultiVector> &args,
                      const std::map<std::string, std::string> &translator ) const
{
    auto mapargs = make_map( args, translator );
    evalv( r, mapargs );
}
void Property::evalv(
    std::vector<std::shared_ptr<std::vector<double>>> &r,
    const std::map<std::string, std::shared_ptr<std::vector<double>>> &args ) const
{
    std::vector<argumentDataStruct<std::vector<double>>> args2;
    for ( const auto &tmp : args )
        args2.emplace_back( tmp.first, *tmp.second );
    evalv( r, {}, args2 );
}


/************************************************************************
 *  evalv (vector)                                                       *
 ************************************************************************/
template<class OUT, class IN>
void Property::evalv( AMP::Array<std::shared_ptr<OUT>> &r,
                      const Units &unit,
                      const std::vector<argumentDataStruct<IN>> &args ) const
{
    AMP_ASSERT( d_dim == r.size() );
    size_t N = getSize( *r( 0 ) );

    // Load the arguments into an array
    auto args2 = loadArgs( N, args );

    // Allocate temporary output array
    ArraySize rs = d_dim;
    rs.setNdim( d_dim.ndim() + 1 );
    rs.resize( d_dim.ndim(), N );
    Array<double> r2( rs );
    r2.fill( 0 );

    // Call eval
    eval( r2, args2 );

    // Convert units
    if ( !unit.isNull() )
        r2.scale( d_units.convert( unit ) );

    // Copy the results back
    size_t N0 = r.length();
    for ( size_t i = 0; i < N0; i++ ) {
        AMP_ASSERT( getSize( *r( i ) ) == N );
        auto it = r( i )->begin();
        for ( size_t j = 0; j < N; j++, ++it )
            *it = r2( i + j * N0 );
    }
}
void Property::evalv(
    AMP::Array<std::shared_ptr<AMP::LinearAlgebra::Vector>> &r,
    const std::map<std::string, std::shared_ptr<AMP::LinearAlgebra::Vector>> &args ) const
{
    std::vector<argumentDataStruct<AMP::LinearAlgebra::Vector>> args2;
    for ( const auto &tmp : args )
        args2.emplace_back( tmp.first, *tmp.second );
    evalv( r, {}, args2 );
}
void Property::evalv( AMP::Array<std::shared_ptr<AMP::LinearAlgebra::Vector>> &r,
                      const std::shared_ptr<AMP::LinearAlgebra::MultiVector> &args,
                      const std::map<std::string, std::string> &translator ) const
{
    auto mapargs = make_map( args, translator );
    evalv( r, mapargs );
}
void Property::evalv(
    AMP::Array<std::shared_ptr<std::vector<double>>> &r,
    const std::map<std::string, std::shared_ptr<std::vector<double>>> &args ) const
{
    std::vector<argumentDataStruct<std::vector<double>>> args2;
    for ( const auto &tmp : args )
        args2.emplace_back( tmp.first, *tmp.second );
    evalv( r, {}, args2 );
}


} // namespace AMP::Materials


/********************************************************
 *  Explicit instantiations of Array                     *
 ********************************************************/
instantiateArrayConstructors( std::shared_ptr<std::vector<double>> );
instantiateArrayConstructors( std::shared_ptr<AMP::LinearAlgebra::Vector> );
