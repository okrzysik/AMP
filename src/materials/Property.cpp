#include "AMP/materials/Property.h"
#include "AMP/materials/ScalarProperty.h"
#include "AMP/utils/MathExpr.h"
#include "AMP/utils/Utilities.h"

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
    // Set the defaults
    d_defaults.resize( d_arguments.size() );
    for ( size_t i = 0; i < d_arguments.size(); i++ )
        d_defaults[i] = d_ranges[i][0] * ( 1.0000001 );
}


/************************************************************************
 *  Evaluate the input arguments                                         *
 ************************************************************************/
void Property::evalArg( AMP::Array<double> &args,
                        const std::string &name,
                        const Units &unit,
                        const std::vector<double> &v ) const
{
    size_t N = args.size( 1 );
    AMP_INSIST( v.size() == N, "Argument " + name + " size does not match input" );
    int i = get_arg_index( name );
    if ( i >= 0 ) {
        double scale = 1.0;
        if ( !unit.isNull() )
            scale = unit.convert( d_argUnits[i] );
        for ( size_t j = 0; j < N; j++ )
            args( i, j ) = scale * v[j];
    }
}
void Property::evalv( const AMP::Array<double> &args,
                      AMP::Array<std::vector<double> *> &r,
                      const Units &unit ) const
{
    // Allocate temporary output array
    size_t N     = args.size( 1 );
    ArraySize rs = d_dim;
    rs.setNdim( d_dim.ndim() + 1 );
    rs.resize( d_dim.ndim(), N );
    Array<double> r2( rs );
    r2.fill( 0 );

    // Call eval
    eval( r2, args );

    // Convert units
    if ( !unit.isNull() )
        r2.scale( d_units.convert( unit ) );

    // Copy the results back
    size_t N0 = r.length();
    for ( size_t i = 0; i < N0; i++ ) {
        AMP_ASSERT( r( i )->size() == N );
        for ( size_t j = 0; j < N; j++ )
            ( *r( i ) )[j] = r2( i + j * N0 );
    }
}
std::string Property::evalString() const
{
    AMP_ERROR( "Evaluating a property as a string is not supported for " + d_name );
}


/************************************************************************
 *  Create the default arguments                                         *
 ************************************************************************/
AMP::Array<double> Property::defaultArgs( size_t N ) const
{
    AMP::Array<double> args2( d_arguments.size(), N );
    if ( args2.empty() )
        return args2;
    for ( size_t j = 0; j < N; j++ )
        memcpy( &args2( 0, j ), d_defaults.data(), d_defaults.size() * sizeof( double ) );
    return args2;
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
    AMP_ASSERT( args.size( 0 ) == d_arguments.size() );
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
 *  Create a property from a database key data object                    *
 ************************************************************************/
std::unique_ptr<Property> createProperty( const std::string &key, const Database &db )
{
    auto keyData = db.getData( key );
    auto unit    = keyData->unit();
    if ( db.isDatabase( key ) ) {
        // We are dealing with a database
        AMP_ERROR( "Not finished (Database)" );
    } else if ( db.isEquation( key ) ) {
        // We are dealing with an equation
        return std::make_unique<EquationProperty>( key, db.getEquation( key ), unit );
    } else if ( keyData->is_floating_point() ) {
        // We are dealing with a scalar
        auto data = keyData->convertToDouble();
        return std::make_unique<ScalarProperty>( key, data, unit );
    } else {
        AMP_ERROR( "Unknown data type" );
    }
    return nullptr;
}


/************************************************************************
 *  Get default value                                                    *
 ************************************************************************/
double Property::get_default( const std::string &name ) const
{
    for ( size_t i = 0; i < d_arguments.size(); i++ ) {
        if ( name == d_arguments[i] )
            return d_defaults[i];
    }
    return std::numeric_limits<double>::quiet_NaN();
}


} // namespace AMP::Materials


/********************************************************
 *  Explicit instantiations of Array                     *
 ********************************************************/
#include "AMP/utils/Array.hpp"
instantiateArrayConstructors( std::vector<double> * );
instantiateArrayConstructors( std::shared_ptr<std::vector<double>> );
