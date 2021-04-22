#include "TensorProperty.h"
#include "AMP/utils/Utilities.h"

#include <algorithm>

namespace AMP {
namespace Materials {


void TensorProperty::evalv(
    std::vector<std::vector<std::shared_ptr<AMP::LinearAlgebra::Vector>>> &r,
    const std::map<std::string, std::shared_ptr<AMP::LinearAlgebra::Vector>> &args )
{
    bool in = this->in_range( args );
    AMP_INSIST( in, "Property out of range" );
    evalvActual( r, args );
}

void TensorProperty::evalv(
    std::vector<std::vector<std::shared_ptr<AMP::LinearAlgebra::Vector>>> &r,
    const std::shared_ptr<AMP::LinearAlgebra::MultiVector> &args )
{
    auto mapargs = make_map( args );
    evalv( r, mapargs );
}


/**
 * \brief	   Evaluate a tensor of vectors of results
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
template<class INPUT_VTYPE, class RETURN_VTYPE>
void TensorProperty::evalvActual( std::vector<std::vector<std::shared_ptr<RETURN_VTYPE>>> &r,
                                  const std::map<std::string, std::shared_ptr<INPUT_VTYPE>> &args )
{
    size_t rdim0 = r.size();    // dimension 0 of results tensor to return
    size_t rdim1 = r[0].size(); // dimension 1 of results tensor to return

    // Check rows of tensor are same size
    for ( size_t i = 1; i < rdim0; i++ ) {
        AMP_ASSERT( r[i].size() == rdim1 );
    }

    // First we make sure that all of the vectors have something in them
    for ( size_t i = 0; i < rdim0; i++ )
        for ( size_t j = 0; j < rdim1; j++ ) {
            AMP_ASSERT( r[i][j]->begin() != r[i][j]->end() );
        }
    AMP_INSIST( rdim0 == d_dimensions[0] && rdim1 == d_dimensions[1],
                "incorrect dimensionality of output tensor" );

    // Make a tensor of iterators for the results tensor
    typename RETURN_VTYPE::iterator dummy_iter =
        r[0][0]->begin(); // for idiosyncracy of AMP::Vector
    std::vector<std::vector<typename RETURN_VTYPE::iterator>> r_iter(
        rdim0,
        std::vector<typename RETURN_VTYPE::iterator>(
            rdim1, dummy_iter ) ); // Iterator for the results tensor

    // Make a vector of iterators - one for each d_arguments
    std::vector<double> eval_args(
        Property::d_n_arguments ); // list of arguments for each input type
    std::vector<typename INPUT_VTYPE::iterator> parameter_iter;
    std::vector<size_t> parameter_indices;
    std::vector<typename std::map<std::string, std::shared_ptr<INPUT_VTYPE>>::const_iterator>
        parameter_map_iter;

    // Walk through d_arguments and set the iterator at the beginning of the map vector to which it
    // corresponds
    for ( size_t i = 0; i < Property::d_arguments.size(); ++i ) {
        typename std::map<std::string, std::shared_ptr<INPUT_VTYPE>>::const_iterator mapIter;
        mapIter = args.find( Property::d_arguments[i] );
        if ( mapIter == args.end() ) {
            eval_args[i] = Property::d_defaults[i];
        } else {
            parameter_iter.push_back( mapIter->second->begin() );
            parameter_indices.push_back( i );
            parameter_map_iter.push_back( mapIter );
        }
    }
    const size_t npresent = parameter_iter.size();

    for ( size_t i = 0; i < rdim0; i++ )
        for ( size_t j = 0; j < rdim1; j++ ) {
            r_iter[i][j] = r[i][j]->begin();
        }
    bool goAgain = true;
    while ( goAgain ) {
        // Loop through the list of actually present parameter iterators and assign their values to
        // the vector being
        // sent to eval
        // Check that parameter iterators have not gone off the end - meaning result and input sizes
        // do not match
        if ( Property::d_n_arguments > 0 ) {
            for ( size_t ipresent = 0; ipresent < npresent; ipresent++ ) {
                AMP_INSIST( parameter_iter[ipresent] != parameter_map_iter[ipresent]->second->end(),
                            std::string( "size mismatch between results and arguments - too few "
                                         "argument values for argument " ) +
                                Property::d_arguments[parameter_indices[ipresent]] +
                                std::string( "\n" ) );
                eval_args[parameter_indices[ipresent]] = *( parameter_iter[ipresent] );
            }
        }
        std::vector<std::vector<double>> evalResult( evalTensor( eval_args ) );
        for ( size_t i = 0; i < rdim0; i++ )
            for ( size_t j = 0; j < rdim1; j++ ) {
                *r_iter[i][j] = evalResult[i][j];
            }

        // update the parameter iterators
        for ( size_t i = 0; i < npresent; i++ ) {
            parameter_iter[i]++;
        }

        // update result iterators;
        for ( size_t i = 0; i < rdim0; i++ )
            for ( size_t j = 0; j < rdim1; j++ ) {
                ++r_iter[i][j];
            }

        // check if any of result iterators reached end
        bool alldone = true;
        for ( size_t i = 0; i < rdim0; i++ )
            for ( size_t j = 0; j < rdim1; j++ ) {
                if ( r_iter[i][j] == r[i][j]->end() )
                    goAgain = false; // if goAgain true, none reached the end
                alldone = alldone &&
                          r_iter[i][j] == r[i][j]->end(); // if alldone true, all reached the end
            }
        // if one reached the end, make sure all did
        if ( !goAgain )
            AMP_INSIST( alldone, "tensor result vectors have unequal sizes" );
    }
    // Make sure the input value iterators all got to the end.
    if ( Property::d_n_arguments > 0 ) {
        for ( size_t ipresent = 0; ipresent < npresent; ipresent++ ) {
            AMP_INSIST( parameter_iter[ipresent] == parameter_map_iter[ipresent]->second->end(),
                        "size mismatch between results and arguments - too few results\n" );
        }
    }
}

void TensorProperty::evalv(
    std::vector<std::vector<std::shared_ptr<std::vector<double>>>> &r,
    const std::map<std::string, std::shared_ptr<std::vector<double>>> &args )
{
    bool in = this->in_range( args );
    AMP_INSIST( in, "Property out of range" );
    evalvActual( r, args );
}

} // namespace Materials
} // namespace AMP
