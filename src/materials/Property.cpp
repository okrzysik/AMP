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
                    std::vector<std::array<double, 2>> ranges )
    : d_name( std::move( name ) ),
      d_source( std::move( source ) ),
      d_params( std::move( params ) ),
      d_arguments( std::move( args ) ),
      d_defaults( {} ),
      d_ranges( std::move( ranges ) ),
      d_variableNumberParameters( false )
{
    AMP_ASSERT( d_arguments.size() == d_ranges.size() );
    for ( size_t i = 0; i < d_arguments.size(); i++ ) {
        d_argToIndexMap.insert( std::pair<std::string, size_t>( d_arguments[i], i ) );
    }
    d_defaults.resize( d_arguments.size() );
    for ( size_t i = 0; i < d_arguments.size(); i++ )
        d_defaults[i] = d_ranges[i][0];
    d_defaultsAreSet = true;
}


/************************************************************************
 *  make_map                                                             *
 ************************************************************************/
std::map<std::string, std::shared_ptr<AMP::LinearAlgebra::Vector>>
Property::make_map( const std::shared_ptr<AMP::LinearAlgebra::MultiVector> &args )
{
    std::map<std::string, std::shared_ptr<AMP::LinearAlgebra::Vector>> result;
    if ( !d_arguments.empty() ) {
        size_t xls = d_translator.size();
        AMP_INSIST( xls > 0, "attempt to make MultiVector map without setting translator" );
        for ( auto vec : *args ) {
            std::string name = vec->getVariable()->getName();
            for ( auto &elem : d_translator ) {
                std::string key = elem.first;
                if ( elem.second == name ) {
                    result.insert( std::make_pair( key, vec ) );
                }
            }
        }
    }
    return result;
}


/************************************************************************
 *  evalvActual                                                          *
 ************************************************************************/
/**
 * \brief       Evaluate a vector of results
 * \tparam      VTYPE       A vector type container.
 *
 *              Template argument:
 *              VTYPE   A vector type that must support:
 *                      VTYPE::iterator
 *                      VTYPE::const_iterator
 *                      VTYPE::iterator begin()
 *                      VTYPE::iterator end()
 *              Constraints:
 *                      type of *VTYPE::iterator is a Number.
 *                      type of *VTYPE::const_iterator is a const Number.
 * \param[out]  r      vector of results
 * \param[in]   args  input arguments corresponding to d_sequence
 *                      Must be in the correct order: T, u, burn
 */
template<class INPUT_VTYPE, class RETURN_VTYPE>
void Property::evalvActual( RETURN_VTYPE &r,
                            const std::map<std::string, std::shared_ptr<INPUT_VTYPE>> &args )
{
    std::vector<double> eval_args( d_arguments.size() ); // list of arguments for each input type

    // First we make sure that all of the vectors have something in them
    AMP_ASSERT( r.begin() != r.end() );

    // Make a vector of iterators - one for each d_arguments
    std::vector<typename INPUT_VTYPE::iterator> parameter_iter;
    std::vector<size_t> parameter_indices;
    std::vector<typename std::map<std::string, std::shared_ptr<INPUT_VTYPE>>::const_iterator>
        parameter_map_iter;

    // Walk through d_arguments and set the iterator at the beginning of the map vector to which it
    // corresponds
    for ( size_t i = 0; i < d_arguments.size(); ++i ) {
        auto mapIter = args.find( d_arguments[i] );
        if ( mapIter == args.end() ) {
            eval_args[i] = d_defaults[i];
        } else {
            parameter_iter.push_back( mapIter->second->begin() );
            parameter_indices.push_back( i );
            parameter_map_iter.push_back( mapIter );
        }
    }
    const size_t npresent = parameter_iter.size();

    for ( auto r_iter = r.begin(); r_iter != r.end(); ++r_iter ) {
        // Loop through the list of actually present parameter iterators and assign their values to
        // the vector being sent to eval
        // Check that parameter iterators have not gone off the end - meaning result and input sizes
        // do not match
        if ( !d_arguments.empty() ) {
            for ( size_t ipresent = 0; ipresent < npresent; ipresent++ ) {
                AMP_INSIST( parameter_iter[ipresent] != parameter_map_iter[ipresent]->second->end(),
                            std::string( "size mismatch between results and arguments - too few "
                                         "argument values for argument " ) +
                                d_arguments[parameter_indices[ipresent]] + std::string( "\n" ) );
                eval_args[parameter_indices[ipresent]] = *( parameter_iter[ipresent] );
            }
        }
        *r_iter = eval( eval_args );

        // update the parameter iterators
        for ( size_t i = 0; i < npresent; i++ ) {
            parameter_iter[i]++;
        }
    }
    // Make sure the input value iterators all got to the end.
    if ( !d_arguments.empty() ) {
        for ( size_t ipresent = 0; ipresent < npresent; ipresent++ ) {
            AMP_INSIST( parameter_iter[ipresent] == parameter_map_iter[ipresent]->second->end(),
                        "size mismatch between results and arguments - too few results\n" );
        }
    }
}


/************************************************************************
 *  evalv                                                                *
 ************************************************************************/
double Property::eval( const std::vector<double> & )
{
    AMP_INSIST( false, "function is not implemented for this property" );
    return 0;
}
void Property::evalv( std::vector<double> &r,
                      const std::map<std::string, std::shared_ptr<std::vector<double>>> &args )
{
    if ( !in_range( args ) ) {
        for ( const auto &arg : args ) {
            if ( !in_range( arg.first, *( arg.second ) ) ) {
                auto range  = get_arg_range( arg.first );
                auto values = *( arg.second );
                std::stringstream ss;
                ss << "Property '" + arg.first + "' out of range in function '" + d_name + "'."
                   << std::endl;
                ss << "Values are ";
                for ( auto &value : values )
                    ss << value << " ";
                ss << std::endl;
                ss << "Valid range is [" << range[0] << "," << range[1] << "]" << std::endl;
                AMP_ERROR( ss.str() );
            }
        }
    }
    evalvActual( r, args );
}
void Property::evalv(
    std::shared_ptr<AMP::LinearAlgebra::Vector> &r,
    const std::map<std::string, std::shared_ptr<AMP::LinearAlgebra::Vector>> &args )
{
    bool in = in_range( args );
    AMP_INSIST( in, "Property out of range" );
    evalvActual( *r, args );
}
void Property::evalv( std::shared_ptr<AMP::LinearAlgebra::Vector> &r,
                      const std::shared_ptr<AMP::LinearAlgebra::MultiVector> &args )
{
    auto mapargs = make_map( args );
    evalv( r, mapargs );
}


/************************************************************************
 *  Get/Set auxiliary data                                               *
 ************************************************************************/
void Property::setAuxiliaryData( const std::string &key, const double val )
{
    auto loc = d_AuxiliaryDataDouble.find( key );

    // this guarantees that the host code can not set just any old data value.
    // the property constructor must set up the database for legitimate use.
    AMP_ASSERT( loc != d_AuxiliaryDataDouble.end() );

    loc->second = val;
}
void Property::setAuxiliaryData( const std::string &key, const int val )
{
    auto loc = d_AuxiliaryDataInteger.find( key );

    // this guarantees that the host code can not set just any old data value.
    // the property constructor must set up the database for legitimate use.
    AMP_ASSERT( loc != d_AuxiliaryDataInteger.end() );

    loc->second = val;
}
void Property::setAuxiliaryData( const std::string &key, const std::string &val )
{
    auto loc = d_AuxiliaryDataString.find( key );

    // this guarantees that the host code can not set just any old data value.
    // the property constructor must set up the database for legitimate use.
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
void Property::set_translator( const std::map<std::string, std::string> &xlator )
{
    // assure incoming map has correct keys
    for ( auto iter = xlator.begin(); iter != xlator.end(); ++iter ) {
        AMP_ASSERT( std::find( d_arguments.begin(), d_arguments.end(), iter->first ) !=
                    d_arguments.end() );
    }
    d_translator = xlator;
}
std::array<double, 2> Property::get_arg_range( const std::string &argname )
{
    auto it = d_argToIndexMap.find( argname );
    if ( it == d_argToIndexMap.end() )
        return { -std::numeric_limits<double>::max(), std::numeric_limits<double>::max() };
    size_t index = it->second;
    return d_ranges[index];
}


} // namespace AMP::Materials
