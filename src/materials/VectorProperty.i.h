#include <algorithm>

namespace AMP {
namespace Materials {

// ====================================================================================
// FUNCTIONS
// ====================================================================================

/**
 * \brief	   Evaluate a vector of vectors of results
 * \tparam	  VTYPE	   A vector type container.
 *
 *			  Template argument:
 *			  VTYPE   A vector type that must support:
 *					  VTYPE::iterator
 *					  VTYPE::const_iterator
 *					  VTYPE::iterator begin()
 *					  VTYPE::iterator end()
 *			  Constraints:
 *					  type of *VTYPE::iterator is a Number.
 *					  type of *VTYPE::const_iterator is a const Number.
 * \param[out]  r	  vector of vectors of results
 * \param[in]   args  input arguments corresponding to d_sequence
 *					  Must be in the correct order: T, u, burn
 */
template <class Number>
template <class INPUT_VTYPE, class RETURN_VTYPE>
void VectorProperty<Number>::evalvActual(
    std::vector<AMP::shared_ptr<RETURN_VTYPE>> &r,
    const std::map<std::string, AMP::shared_ptr<INPUT_VTYPE>> &args )
{
    size_t rdim0 = r.size(); // number of results vectors to return

    // First we make sure that all of the vectors have something in them
    for ( size_t i = 0; i < rdim0; i++ ) {
        AMP_ASSERT( r[i]->begin() != r[i]->end() );
    }
    AMP_INSIST( rdim0 == d_dimension, "incorrect dimensionality for output vector" );

    // Make a vector of iterators for the results vector
    typename RETURN_VTYPE::iterator dummy_iter = r[0]->begin(); // for idiosyncracy of AMP::Vector
    std::vector<typename RETURN_VTYPE::iterator> r_iter( rdim0, dummy_iter );

    // Make a vector of iterators - one for each d_arguments
    std::vector<Number> eval_args(
        Property<Number>::d_n_arguments ); // list of arguments for each input type
    std::vector<typename INPUT_VTYPE::iterator> parameter_iter;
    std::vector<size_t> parameter_indices;
    std::vector<typename std::map<std::string, AMP::shared_ptr<INPUT_VTYPE>>::const_iterator>
        parameter_map_iter;

    // Walk through d_arguments and set the iterator at the beginning of the map vector to which it
    // corresponds
    for ( size_t i = 0; i < Property<Number>::d_arguments.size(); ++i ) {
        typename std::map<std::string, AMP::shared_ptr<INPUT_VTYPE>>::const_iterator mapIter;
        mapIter = args.find( Property<Number>::d_arguments[i] );
        if ( mapIter == args.end() ) {
            eval_args[i] = Property<Number>::d_defaults[i];
        }
        else {
            parameter_iter.push_back( mapIter->second->begin() );
            parameter_indices.push_back( i );
            parameter_map_iter.push_back( mapIter );
        }
    }
    const size_t npresent = parameter_iter.size();

    for ( size_t i = 0; i < rdim0; i++ ) {
        r_iter[i] = r[i]->begin();
    }
    bool goAgain = true;
    while ( goAgain ) {
        // Loop through the list of actually present parameter iterators and assign their values to
        // the vector being
        // sent to eval
        // Check that parameter iterators have not gone off the end - meaning result and input sizes
        // do not match
        if ( Property<Number>::d_n_arguments > 0 ) {
            for ( size_t ipresent = 0; ipresent < npresent; ipresent++ ) {
                AMP_INSIST( parameter_iter[ipresent] != parameter_map_iter[ipresent]->second->end(),
                            std::string( "size mismatch between results and arguments - too few "
                                         "argument values for argument " ) +
                                Property<Number>::d_arguments[parameter_indices[ipresent]] +
                                std::string( "\n" ) );
                eval_args[parameter_indices[ipresent]] = *( parameter_iter[ipresent] );
            }
        }
        std::vector<Number> evalResult = evalVector( eval_args );
        for ( size_t i = 0; i < rdim0; i++ ) {
            *( r_iter[i] ) = evalResult[i];
        }

        // update the parameter iterators
        for ( size_t i = 0; i < npresent; i++ ) {
            parameter_iter[i]++;
        }

        // update result iterators;
        for ( size_t i = 0; i < rdim0; i++ ) {
            ++( r_iter[i] );
        }

        // check if any of result iterators reached end
        bool alldone = true;
        for ( size_t i = 0; i < rdim0; i++ ) {
            if ( r_iter[i] == r[i]->end() )
                goAgain = false;                           // if goAgain true, none reached the end
            alldone = alldone && r_iter[i] == r[i]->end(); // if alldone true, all reached the end
        }
        // if one reached the end, make sure all did
        if ( !goAgain ) AMP_INSIST( alldone, "vector result vectors have unequal sizes" );
    }
    // Make sure the input value iterators all got to the end.
    if ( Property<Number>::d_n_arguments > 0 ) {
        for ( size_t ipresent = 0; ipresent < npresent; ipresent++ ) {
            AMP_INSIST( parameter_iter[ipresent] == parameter_map_iter[ipresent]->second->end(),
                        "size mismatch between results and arguments - too few results\n" );
        }
    }
}

template <class Number>
void VectorProperty<Number>::evalv(
    std::vector<AMP::shared_ptr<std::vector<Number>>> &r,
    const std::map<std::string, AMP::shared_ptr<std::vector<Number>>> &args )
{
    AMP_ASSERT( this->in_range( args ) );
    evalvActual( r, args );
}

} // Materials namespace
} // AMP namespace
