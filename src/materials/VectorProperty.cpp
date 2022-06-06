#include "VectorProperty.h"
#include "AMP/utils/Utilities.h"

#include <algorithm>

namespace AMP::Materials {


/************************************************************************
 *  Constructor                                                          *
 ************************************************************************/
VectorProperty::VectorProperty( std::string name,
                                std::string source,
                                std::vector<std::string> args,
                                std::vector<std::array<double, 2>> ranges,
                                size_t dimension )
    : Property(
          std::move( name ), Units(), std::move( source ), std::move( args ), std::move( ranges ) ),
      d_dimension( dimension )
{
    AMP_INSIST( d_dimension > 0, "must return at least one value" );
}


/************************************************************************
 *  evalv (vector)                                                       *
 ************************************************************************/
static inline size_t size( const std::vector<double> &x ) { return x.size(); }
static inline size_t size( const AMP::LinearAlgebra::Vector &x ) { return x.getLocalSize(); }
template<class OUT, class IN>
void VectorProperty::evalv( std::vector<std::shared_ptr<OUT>> &r,
                            const std::vector<argumentDataStruct<IN>> &args ) const
{
    // Check that argument vectors are the same size
    for ( size_t i = 0; i < args.size(); i++ )
        AMP_ASSERT( size( args[i].vec ) == size( *r[0] ) );

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

    // Get the iterators to the results
    size_t N = size( *r[0] );
    std::vector<decltype( r[0]->begin() )> r_it( r.size() );
    for ( size_t i = 0; i < r.size(); i++ ) {
        AMP_ASSERT( size( *r[i] ) == N );
        r_it[i] = r[i]->begin();
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
        auto result = evalVector( eval_args );
        for ( size_t i = 0; i < r.size(); i++ ) {
            *r_it[i] = scale * result[i];
            ++r_it[i];
        }
    }
}

void VectorProperty::evalv(
    std::vector<std::shared_ptr<AMP::LinearAlgebra::Vector>> &r,
    const std::map<std::string, std::shared_ptr<AMP::LinearAlgebra::Vector>> &args ) const
{
    std::vector<argumentDataStruct<AMP::LinearAlgebra::Vector>> args2;
    for ( const auto &tmp : args )
        args2.emplace_back( tmp.first, *tmp.second );
    evalv( r, args2 );
}
void VectorProperty::evalv( std::vector<std::shared_ptr<AMP::LinearAlgebra::Vector>> &r,
                            const std::shared_ptr<AMP::LinearAlgebra::MultiVector> &args,
                            const std::map<std::string, std::string> &translator ) const
{
    auto mapargs = make_map( args, translator );
    evalv( r, mapargs );
}
void VectorProperty::evalv(
    std::vector<std::shared_ptr<std::vector<double>>> &r,
    const std::map<std::string, std::shared_ptr<std::vector<double>>> &args ) const
{
    std::vector<argumentDataStruct<std::vector<double>>> args2;
    for ( const auto &tmp : args )
        args2.emplace_back( tmp.first, *tmp.second );
    evalv( r, args2 );
}

} // namespace AMP::Materials
