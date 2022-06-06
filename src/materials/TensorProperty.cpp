#include "TensorProperty.h"
#include "AMP/utils/Array.hpp"
#include "AMP/utils/Utilities.h"

#include <algorithm>

namespace AMP::Materials {


/************************************************************************
 *  Constructor                                                          *
 ************************************************************************/
TensorProperty::TensorProperty( std::string name,
                                std::string source,
                                std::vector<std::string> args,
                                std::vector<std::array<double, 2>> ranges,
                                const AMP::ArraySize &dimensions )
    : Property(
          std::move( name ), Units(), std::move( source ), std::move( args ), std::move( ranges ) ),
      d_dimensions( dimensions )
{
    AMP_INSIST( d_dimensions.size() == 2, "there must be two dimensions" );
    AMP_INSIST( d_dimensions[0] > 0, "must have first return tensor dimension > 0" );
    AMP_INSIST( d_dimensions[1] > 0, "must have second return tensor dimension > 0" );
}


/************************************************************************
 *  evalv (vector)                                                       *
 ************************************************************************/
static inline size_t size( const std::vector<double> &x ) { return x.size(); }
static inline size_t size( const AMP::LinearAlgebra::Vector &x ) { return x.getLocalSize(); }
template<class OUT, class IN>
void TensorProperty::evalv( AMP::Array<std::shared_ptr<OUT>> &r,
                            const std::vector<argumentDataStruct<IN>> &args ) const
{

    // Check that argument vectors are the same size
    for ( size_t i = 0; i < args.size(); i++ )
        AMP_ASSERT( size( args[i].vec ) == size( *r( 0 ) ) );

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

    // Check rows of tensor are same size
    size_t rdim0 = r.size( 0 );
    size_t rdim1 = r.size( 1 );

    // Get the iterators to the results
    size_t N   = size( *r( 0 ) );
    using R_IT = decltype( r( 0 )->begin() );
    AMP::Array<R_IT> r_it( rdim0, rdim1 );
    for ( size_t i = 0; i < rdim0; i++ ) {
        for ( size_t j = 0; j < rdim1; j++ ) {
            AMP_ASSERT( size( *r( i, j ) ) == N );
            r_it( i, j ) = r( i, j )->begin();
        }
    }

    // Call eval for each entry
    double scale = 1.0;
    // if ( !unit.isNull() )
    //    scale = d_units.convert( unit );
    for ( size_t k = 0; k < N; k++ ) {
        // Update the arguments
        for ( size_t i = 0; i < args.size(); i++ ) {
            if ( index[i] >= 0 ) {
                eval_args[i] = scaleArgs[i] * ( *arg_it[i] );
                ++arg_it[i];
            }
        }
        // Call eval
        auto result = evalTensor( eval_args );
        for ( size_t i = 0; i < rdim0; i++ ) {
            for ( size_t j = 0; j < rdim1; j++ ) {
                *r_it( i, j ) = scale * result( i, j );
                ++r_it( i, j );
            }
        }
    }
}
void TensorProperty::evalv(
    AMP::Array<std::shared_ptr<AMP::LinearAlgebra::Vector>> &r,
    const std::map<std::string, std::shared_ptr<AMP::LinearAlgebra::Vector>> &args ) const
{
    std::vector<argumentDataStruct<AMP::LinearAlgebra::Vector>> args2;
    for ( const auto &tmp : args )
        args2.emplace_back( tmp.first, *tmp.second );
    evalv( r, args2 );
}
void TensorProperty::evalv( AMP::Array<std::shared_ptr<AMP::LinearAlgebra::Vector>> &r,
                            const std::shared_ptr<AMP::LinearAlgebra::MultiVector> &args,
                            const std::map<std::string, std::string> &translator ) const
{
    auto mapargs = make_map( args, translator );
    evalv( r, mapargs );
}
void TensorProperty::evalv(
    AMP::Array<std::shared_ptr<std::vector<double>>> &r,
    const std::map<std::string, std::shared_ptr<std::vector<double>>> &args ) const
{
    std::vector<argumentDataStruct<std::vector<double>>> args2;
    for ( const auto &tmp : args )
        args2.emplace_back( tmp.first, *tmp.second );
    evalv( r, args2 );
}


} // namespace AMP::Materials


/********************************************************
 *  Explicit instantiations of Array                     *
 ********************************************************/
instantiateArrayConstructors( std::shared_ptr<std::vector<double>> );
instantiateArrayConstructors( std::shared_ptr<AMP::LinearAlgebra::Vector> );
