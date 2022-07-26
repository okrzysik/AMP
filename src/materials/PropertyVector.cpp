#include "AMP/materials/Property.h"
#include "AMP/vectors/MultiVector.h"
#include "AMP/vectors/Vector.h"


namespace AMP::Materials {


/************************************************************************
 *  Evaluate the input arguments                                         *
 ************************************************************************/
size_t Property::getSize( const AMP::LinearAlgebra::Vector &v ) { return v.getLocalSize(); }
Units Property::getUnits( const AMP::LinearAlgebra::Vector &v ) { return v.getUnits(); }
void Property::evalArgs( AMP::Array<double> &args2,
                         const std::shared_ptr<AMP::LinearAlgebra::MultiVector> &args,
                         const std::map<std::string, std::string> &translator ) const
{
    auto mapargs = make_map( args, translator );
    for ( auto arg : mapargs )
        evalArgs( args2, arg.first, *arg.second );
}
void Property::evalArg( AMP::Array<double> &args,
                        const std::string &name,
                        const Units &unit,
                        const AMP::LinearAlgebra::Vector &v ) const
{
    size_t N = args.size( 1 );
    AMP_INSIST( v.getLocalSize() == N, "Argument " + name + " size does not match input" );
    int i = get_arg_index( name );
    if ( i >= 0 ) {
        double scale = 1.0;
        if ( !unit.isNull() )
            scale = unit.convert( d_argUnits[i] );
        auto it = v.begin();
        for ( size_t j = 0; j < N; j++, ++it )
            args( i, j ) = scale * *it;
    }
}
void Property::evalv( const AMP::Array<double> &args,
                      AMP::Array<AMP::LinearAlgebra::Vector *> &r ) const
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

    // Copy the results back
    size_t N0 = r.length();
    for ( size_t i = 0; i < N0; i++ ) {
        double scale = 1.0;
        if ( !r( i )->getUnits().isNull() )
            scale = d_units.convert( r( i )->getUnits() );
        AMP_ASSERT( r( i )->getLocalSize() == N );
        auto it = r( i )->begin();
        for ( size_t j = 0; j < N; j++, ++it )
            *it = scale * r2( i + j * N0 );
    }
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


} // namespace AMP::Materials


/********************************************************
 *  Explicit instantiations of Array                     *
 ********************************************************/
#include "AMP/utils/Array.hpp"
instantiateArrayConstructors( AMP::LinearAlgebra::Vector * );
instantiateArrayConstructors( std::shared_ptr<AMP::LinearAlgebra::Vector> );
