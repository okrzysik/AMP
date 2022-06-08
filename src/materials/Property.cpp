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
 *  Evaluate the input arguments                                         *
 ************************************************************************/
void Property::evalArgs( AMP::Array<double> &args2,
                         const std::shared_ptr<AMP::LinearAlgebra::MultiVector> &args,
                         const std::map<std::string, std::string> &translator ) const
{
    auto mapargs = make_map( args, translator );
    for ( auto arg : mapargs )
        evalArgs( args2, arg.first, *arg.second );
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


} // namespace AMP::Materials


/********************************************************
 *  Explicit instantiations of Array                     *
 ********************************************************/
instantiateArrayConstructors( std::vector<double> * );
instantiateArrayConstructors( AMP::LinearAlgebra::Vector * );
instantiateArrayConstructors( std::shared_ptr<std::vector<double>> );
instantiateArrayConstructors( std::shared_ptr<AMP::LinearAlgebra::Vector> );
