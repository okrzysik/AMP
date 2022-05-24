#include "Property.h"
#include "AMP/utils/Utilities.h"
#include "AMP/vectors/MultiVector.h"
#include "AMP/vectors/Vector.h"

#include <algorithm>


namespace AMP::Materials {


/************************************************************************
 *  Constructor                                                          *
 ************************************************************************/
Property::Property( std::string name,
                    std::string source,
                    std::vector<double> params,
                    std::vector<std::string> args,
                    std::vector<std::array<double, 2>> ranges,
                    std::vector<Units> units )
    : d_name( std::move( name ) ),
      d_source( std::move( source ) ),
      d_params( std::move( params ) ),
      d_arguments( std::move( args ) ),
      d_argUnits( std::move( units ) ),
      d_defaults( {} ),
      d_ranges( std::move( ranges ) ),
      d_variableNumberParameters( false )
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
 *  evalvActual                                                          *
 ************************************************************************/
static inline size_t size( const std::vector<double> &x ) { return x.size(); }
static inline size_t size( const AMP::LinearAlgebra::Vector &x ) { return x.getLocalSize(); }
template<class OUT, class IN>
void Property::evalv( OUT &r, const std::vector<argumentDataStruct<IN>> &args )
{
    // Check that all vectors are the same size
    for ( size_t i = 0; i < args.size(); i++ )
        AMP_ASSERT( size( args[i].vec ) == size( r ) );

    // Check that all arguments are in range
    for ( const auto &arg : args )
        in_range( std::string( arg.str ), arg.vec, arg.units, true );

    // Set the default values
    auto eval_args = d_defaults;

    // Get the indices and iterators to override
    std::vector<int> index( args.size(), -1 );
    std::vector<double> scaleArgs( args.size(), 1.0 );
    std::vector<typename IN::const_iterator> arg_it( args.size() );
    for ( size_t i = 0; i < args.size(); i++ ) {
        auto it = d_argToIndexMap.find( args[i].str );
        if ( it != d_argToIndexMap.end() ) {
            size_t j  = it->second;
            index[j]  = i;
            arg_it[j] = args[i].vec.begin();
            if ( !args[i].units.isNull() && !d_argUnits[j].isNull() )
                scaleArgs[j] = args[i].units.convert( d_argUnits[j] );
        }
    }

    // Call eval for each entry
    for ( auto r_it = r.begin(); r_it != r.end(); ++r_it ) {
        // Update the arguments
        for ( size_t i = 0; i < args.size(); i++ ) {
            if ( index[i] >= 0 ) {
                eval_args[i] = scaleArgs[i] * ( *arg_it[i] );
                ++arg_it[i];
            }
        }
        // Call eval
        *r_it = eval( eval_args );
    }
}
double Property::eval( const std::vector<std::string> &names,
                       const std::vector<double> &args,
                       const std::vector<Units> &units )
{
    // Set the default values
    auto eval_args = d_defaults;
    // Override default arguments
    AMP_ASSERT( names.size() == args.size() );
    for ( size_t i = 0; i < args.size(); i++ ) {
        for ( size_t j = 0; j < d_arguments.size(); j++ ) {
            if ( names[i] == d_arguments[j] ) {
                eval_args[j] = args[i];
                if ( !units.empty() ) {
                    if ( !units[i].isNull() && !d_argUnits[j].isNull() )
                        eval_args[j] *= units[i].convert( d_argUnits[j] );
                }
            }
        }
    }
    return eval( eval_args );
}
template void Property::evalv( std::vector<double> &,
                               const std::vector<argumentDataStruct<std::vector<double>>> & );


/************************************************************************
 *  evalv                                                                *
 ************************************************************************/
double Property::eval( const std::vector<double> & )
{
    AMP_INSIST( false, "function is not implemented for this property" );
    return 0;
}
void Property::evalv(
    std::shared_ptr<AMP::LinearAlgebra::Vector> &r,
    const std::map<std::string, std::shared_ptr<AMP::LinearAlgebra::Vector>> &args )
{
    std::vector<argumentDataStruct<AMP::LinearAlgebra::Vector>> args2;
    for ( const auto &t : args )
        args2.emplace_back( std::get<0>( t ), *std::get<1>( t ), std::get<1>( t )->getUnits() );
    evalv( *r, args2 );
}
void Property::evalv( std::shared_ptr<AMP::LinearAlgebra::Vector> &r,
                      const std::shared_ptr<AMP::LinearAlgebra::MultiVector> &args,
                      const std::map<std::string, std::string> &translator )
{
    auto mapargs = make_map( args, translator );
    evalv( r, mapargs );
}
void Property::evalv( std::vector<double> &r,
                      const std::map<std::string, std::shared_ptr<std::vector<double>>> &args )
{
    std::vector<argumentDataStruct<std::vector<double>>> args2;
    for ( const auto &tmp : args )
        args2.emplace_back( tmp.first, *tmp.second );
    evalv( r, args2 );
}


/************************************************************************
 *  Get/Set auxiliary data                                               *
 ************************************************************************/
void Property::setAuxiliaryData( const std::string &key, const double val )
{
    auto loc = d_AuxiliaryDataDouble.find( key );
    AMP_ASSERT( loc != d_AuxiliaryDataDouble.end() );
    loc->second = val;
}
void Property::setAuxiliaryData( const std::string &key, const int val )
{
    auto loc = d_AuxiliaryDataInteger.find( key );
    AMP_ASSERT( loc != d_AuxiliaryDataInteger.end() );
    loc->second = val;
}
void Property::setAuxiliaryData( const std::string &key, const std::string &val )
{
    auto loc = d_AuxiliaryDataString.find( key );
    AMP_ASSERT( loc != d_AuxiliaryDataString.end() );
    loc->second = val;
}
void Property::getAuxiliaryData( const std::string &key, double &val )
{
    auto p = d_AuxiliaryDataDouble.find( key );
    AMP_ASSERT( p != d_AuxiliaryDataDouble.end() );
    val = p->second;
}
void Property::getAuxiliaryData( const std::string &key, int &val )
{
    auto p = d_AuxiliaryDataInteger.find( key );
    AMP_ASSERT( p != d_AuxiliaryDataInteger.end() );
    val = p->second;
}
void Property::getAuxiliaryData( const std::string &key, std::string &val )
{
    auto p = d_AuxiliaryDataString.find( key );
    AMP_ASSERT( p != d_AuxiliaryDataString.end() );
    val = p->second;
}

/************************************************************************
 *  Misc functions                                                       *
 ************************************************************************/
void Property::set_parameters_and_number( std::vector<double> params )
{
    AMP_INSIST( d_variableNumberParameters,
                "changing number of parameters for this property not allowed" );
    d_params = std::move( params );
}
std::array<double, 2> Property::get_arg_range( const std::string &name )
{
    int index = get_arg_index( name );
    AMP_ASSERT( index >= 0 );
    return d_ranges[index];
}
bool Property::is_argument( const std::string &name ) { return get_arg_index( name ) >= 0; }


} // namespace AMP::Materials
