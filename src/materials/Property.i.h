#include <algorithm>
#include <sstream>

namespace AMP {
namespace Materials {

// ====================================================================================
// FUNCTIONS
// ====================================================================================

// set auxiliary data
template<class Number>
void Property<Number>::setAuxiliaryData( const std::string &key, const double val )
{
    std::map<std::string, double>::iterator loc = d_AuxiliaryDataDouble.find( key );

    // this guarantees that the host code can not set just any old data value.
    // the property constructor must set up the database for legitimate use.
    AMP_ASSERT( loc != d_AuxiliaryDataDouble.end() );

    loc->second = val;
}

// set auxiliary data
template<class Number>
void Property<Number>::setAuxiliaryData( const std::string &key, const int val )
{
    std::map<std::string, int>::iterator loc = d_AuxiliaryDataInteger.find( key );

    // this guarantees that the host code can not set just any old data value.
    // the property constructor must set up the database for legitimate use.
    AMP_ASSERT( loc != d_AuxiliaryDataInteger.end() );

    loc->second = val;
}

// set auxiliary data
template<class Number>
void Property<Number>::setAuxiliaryData( const std::string &key, const std::string &val )
{
    std::map<std::string, std::string>::iterator loc = d_AuxiliaryDataString.find( key );

    // this guarantees that the host code can not set just any old data value.
    // the property constructor must set up the database for legitimate use.
    AMP_ASSERT( loc != d_AuxiliaryDataString.end() );

    loc->second = val;
}

/// get auxiliary data
template<class Number>
void Property<Number>::getAuxiliaryData( const std::string &key, double &val )
{
    val = *d_AuxiliaryDataDouble.find( key );
}

/// get auxiliary data
template<class Number>
void Property<Number>::getAuxiliaryData( const std::string &key, int &val )
{
    val = *d_AuxiliaryDataInteger.find( key );
}

/// get auxiliary data
template<class Number>
void Property<Number>::getAuxiliaryData( const std::string &key, std::string &val )
{
    val = *d_AuxiliaryDataString.find( key );
}


/************************************************************************
 *  Determine if a set of values are all within range or not             *
 ************************************************************************/
template<class Number>
bool Property<Number>::in_range( const std::string &argname, const Number value )
{
    if ( !is_argument( argname ) )
        return true;
    std::vector<Number> range = get_arg_range( argname );
    return value >= range[0] && value <= range[1];
}
template<class Number>
template<class INPUT_VTYPE>
bool Property<Number>::in_range( const std::string &argname, const INPUT_VTYPE &values )
{
    if ( !is_argument( argname ) )
        return true;
    std::vector<Number> range                = get_arg_range( argname );
    bool result                              = true;
    typename INPUT_VTYPE::const_iterator pos = values.begin();
    typename INPUT_VTYPE::const_iterator end = values.end();
    while ( pos != end ) {
        result = result && *pos >= range[0] && *pos <= range[1];
        ++pos;
    }
    return result;
}
template<class Number>
template<class INPUT_VTYPE>
bool Property<Number>::in_range( const std::map<std::string, AMP::shared_ptr<INPUT_VTYPE>> &values )
{
    bool result = true;
    for ( const auto &value : values ) {
        result = result && in_range( value.first, *( value.second ) );
    }
    return result;
}

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
template<class Number>
template<class INPUT_VTYPE, class RETURN_VTYPE>
void Property<Number>::evalvActual(
    RETURN_VTYPE &r, const std::map<std::string, AMP::shared_ptr<INPUT_VTYPE>> &args )
{
    std::vector<Number> eval_args( d_n_arguments ); // list of arguments for each input type

    // First we make sure that all of the vectors have something in them
    AMP_ASSERT( r.begin() != r.end() );

    // Make a vector of iterators - one for each d_arguments
    std::vector<typename INPUT_VTYPE::iterator> parameter_iter;
    std::vector<size_t> parameter_indices;
    std::vector<typename std::map<std::string, AMP::shared_ptr<INPUT_VTYPE>>::const_iterator>
        parameter_map_iter;

    // Walk through d_arguments and set the iterator at the beginning of the map vector to which it
    // corresponds
    for ( size_t i = 0; i < d_arguments.size(); ++i ) {
        typename std::map<std::string, AMP::shared_ptr<INPUT_VTYPE>>::const_iterator mapIter;
        mapIter = args.find( d_arguments[i] );
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
        // the vector being
        // sent to eval
        // Check that parameter iterators have not gone off the end - meaning result and input sizes
        // do not match
        if ( d_n_arguments > 0 ) {
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
    if ( d_n_arguments > 0 ) {
        for ( size_t ipresent = 0; ipresent < npresent; ipresent++ ) {
            AMP_INSIST( parameter_iter[ipresent] == parameter_map_iter[ipresent]->second->end(),
                        "size mismatch between results and arguments - too few results\n" );
        }
    }
}

template<class Number>
Number Property<Number>::eval( std::vector<Number> & )
{
    AMP_INSIST( false, "function is not implemented for this property" );
    Number x = 0.;
    return x;
}

template<class Number>
void Property<Number>::evalv(
    std::vector<Number> &r,
    const std::map<std::string, AMP::shared_ptr<std::vector<Number>>> &args )
{
    if ( !in_range( args ) ) {
        for ( const auto &arg : args ) {
            if ( !in_range( arg.first, *( arg.second ) ) ) {
                std::vector<Number> range   = get_arg_range( arg.first );
                std::vector<Number> &values = *( arg.second );
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


} // namespace Materials
} // namespace AMP
