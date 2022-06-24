#include <algorithm>
#include <sstream>

#include "AMP/vectors/Vector.h"
#include "AMP/vectors/data/VectorDataIterator.h"

namespace AMP::Materials {


/************************************************************************
 *  evalv interfaces                                                     *
 ************************************************************************/
template<class... Args>
double Property::eval( const Units &unit, Args... args ) const
{
    double value = 0;
    auto result  = AMP::Array<double>::staticView( { 1 }, &value );
    // Get the arguments
    double tmp[128];
    memcpy( tmp, d_defaults.data(), d_defaults.size() * sizeof( double ) );
    auto args2 = AMP::Array<double>::staticView( { d_arguments.size() }, tmp );
    evalArgs( args2, args... );
    checkArgs( args2 );
    // Evaluate the property
    eval( result, args2 );
    // Convert units if required
    if ( !unit.isNull() )
        value *= d_units.convert( unit );
    return value;
}
template<class... Args>
void Property::evalv( std::vector<double> &r, const Units &unit, Args... args ) const
{
    AMP::Array<std::vector<double> *> r2( 1 );
    r2( 0 ) = &r;
    evalv( r2, unit, args... );
}
template<class... Args>
void Property::evalv( AMP::LinearAlgebra::Vector &r, Args... args ) const
{
    AMP::Array<AMP::LinearAlgebra::Vector *> r2( 1 );
    r2( 0 ) = &r;
    evalv( r2, args... );
}
template<class... Args>
void Property::evalv( std::vector<std::shared_ptr<std::vector<double>>> &r,
                      const Units &unit,
                      Args... args ) const
{
    AMP::Array<std::vector<double> *> r2( r.size() );
    for ( size_t i = 0; i < r2.length(); i++ )
        r2( i ) = r[i].get();
    evalv( r2, unit, args... );
}
template<class... Args>
void Property::evalv( std::vector<std::shared_ptr<AMP::LinearAlgebra::Vector>> &r,
                      Args... args ) const
{
    AMP::Array<AMP::LinearAlgebra::Vector *> r2( r.size() );
    for ( size_t i = 0; i < r2.length(); i++ )
        r2( i ) = r[i].get();
    evalv( r2, args... );
}
template<class... Args>
void Property::evalv( AMP::Array<std::shared_ptr<std::vector<double>>> &r,
                      const Units &unit,
                      Args... args ) const
{
    AMP::Array<std::vector<double> *> r2( r.size() );
    for ( size_t i = 0; i < r2.length(); i++ )
        r2( i ) = r( i ).get();
    evalv( r2, unit, args... );
}
template<class... Args>
void Property::evalv( AMP::Array<std::shared_ptr<AMP::LinearAlgebra::Vector>> &r,
                      Args... args ) const
{
    AMP::Array<AMP::LinearAlgebra::Vector *> r2( r.size() );
    for ( size_t i = 0; i < r2.length(); i++ )
        r2( i ) = r( i ).get();
    evalv( r2, args... );
}
template<class... Args>
void Property::evalv( AMP::Array<std::vector<double> *> &r, const Units &unit, Args... args ) const
{
    // Load the arguments
    auto args2 = defaultArgs( r( 0 )->size() );
    evalArgs( args2, args... );
    checkArgs( args2 );

    // Allocate temporary output array
    size_t N     = args2.size( 1 );
    ArraySize rs = d_dim;
    rs.setNdim( d_dim.ndim() + 1 );
    rs.resize( d_dim.ndim(), N );
    Array<double> r2( rs );
    r2.fill( 0 );

    // Call eval
    eval( r2, args2 );

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
template<class... Args>
void Property::evalv( AMP::Array<AMP::LinearAlgebra::Vector *> &r, Args... args ) const
{
    // Load the arguments
    auto args2 = defaultArgs( r( 0 )->getLocalSize() );
    evalArgs( args2, args... );
    checkArgs( args2 );

    // Allocate temporary output array
    size_t N     = args2.size( 1 );
    ArraySize rs = d_dim;
    rs.setNdim( d_dim.ndim() + 1 );
    rs.resize( d_dim.ndim(), N );
    Array<double> r2( rs );
    r2.fill( 0 );

    // Call eval
    eval( r2, args2 );

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
 *  Evaluate the input arguments                                         *
 ************************************************************************/
template<class VEC>
void Property::evalArgs( AMP::Array<double> &args2, const std::map<std::string, VEC> &args ) const
{
    for ( const auto &arg : args )
        evalArgs( args2, arg.first, arg.second );
}
template<class... Args>
void Property::evalArgs( AMP::Array<double> &args2,
                         const std::string &name,
                         double v,
                         Args... args ) const
{
    size_t N = args2.size( 1 );
    int i    = get_arg_index( name );
    if ( i >= 0 ) {
        for ( size_t j = 0; j < N; j++ )
            args2( i, j ) = v;
    }
    evalArgs( args2, args... );
}
template<class... Args>
void Property::evalArgs( AMP::Array<double> &args2,
                         const std::string &name,
                         const Units &unit,
                         double v,
                         Args... args ) const
{
    size_t N = args2.size( 1 );
    int i    = get_arg_index( name );
    if ( i >= 0 ) {
        if ( !unit.isNull() )
            v *= unit.convert( d_argUnits[i] );
        for ( size_t j = 0; j < N; j++ )
            args2( i, j ) = v;
    }
    evalArgs( args2, args... );
}
template<class... Args>
void Property::evalArgs( AMP::Array<double> &args2,
                         const std::string &name,
                         const Units &unit,
                         const std::vector<double> &v,
                         Args... args ) const
{
    size_t N = args2.size( 1 );
    AMP_INSIST( v.size() == N, "Argument " + name + " size does not match input" );
    int i = get_arg_index( name );
    if ( i >= 0 ) {
        double scale = 1.0;
        if ( !unit.isNull() )
            scale = unit.convert( d_argUnits[i] );
        for ( size_t j = 0; j < N; j++ )
            args2( i, j ) = scale * v[j];
    }
    evalArgs( args2, args... );
}
template<class... Args>
void Property::evalArgs( AMP::Array<double> &args2,
                         const std::string &name,
                         const Units &unit,
                         const AMP::LinearAlgebra::Vector &v,
                         Args... args ) const
{
    size_t N = args2.size( 1 );
    AMP_INSIST( v.getLocalSize() == N, "Argument " + name + " size does not match input" );
    int i = get_arg_index( name );
    if ( i >= 0 ) {
        double scale = 1.0;
        if ( !unit.isNull() )
            scale = unit.convert( d_argUnits[i] );
        auto it = v.begin();
        for ( size_t j = 0; j < N; j++, ++it )
            args2( i, j ) = scale * *it;
    }
    evalArgs( args2, args... );
}
template<class... Args>
void Property::evalArgs( AMP::Array<double> &args2,
                         const std::string &name,
                         const std::vector<double> &v,
                         Args... args ) const
{
    evalArgs( args2, name, Units(), v, args... );
}
template<class... Args>
void Property::evalArgs( AMP::Array<double> &args2,
                         const std::string &name,
                         const AMP::LinearAlgebra::Vector &v,
                         Args... args ) const
{
    evalArgs( args2, name, v.getUnits(), v, args... );
}
template<class VEC, class... Args>
void Property::evalArgs( AMP::Array<double> &args2,
                         const std::string &name,
                         const Units &unit,
                         const std::shared_ptr<VEC> &v,
                         Args... args ) const
{
    evalArgs( args2, name, unit, *v, args... );
}
template<class VEC, class... Args>
void Property::evalArgs( AMP::Array<double> &args2,
                         const std::string &name,
                         const std::shared_ptr<VEC> &v,
                         Args... args ) const
{
    evalArgs( args2, name, *v, args... );
}
template<class... Args>
void Property::evalArgs( AMP::Array<double> &args2,
                         const std::vector<double> &args,
                         const std::vector<std::string> &names,
                         const std::vector<Units> &argUnits ) const
{
    size_t N = args2.size( 1 );
    AMP_ASSERT( args.size() == names.size() );
    AMP_ASSERT( argUnits.empty() || argUnits.size() == names.size() );
    for ( size_t k = 0; k < args.size(); k++ ) {
        int i = get_arg_index( names[k] );
        if ( i >= 0 ) {
            double scale = 1.0;
            if ( !argUnits.empty() ) {
                if ( !argUnits[k].isNull() )
                    scale = argUnits[k].convert( d_argUnits[i] );
            }
            for ( size_t j = 0; j < N; j++ )
                args2( i, j ) = scale * args[k];
        }
    }
}


/************************************************************************
 *  Determine if a set of values are all within range or not             *
 ************************************************************************/
inline bool
Property::in_range( const std::string &name, double value, Units unit, bool throwError ) const
{
    auto index = get_arg_index( name );
    if ( index == -1 )
        return true;
    double scale = 1.0;
    if ( !unit.isNull() && !d_argUnits[index].isNull() )
        scale = unit.convert( d_argUnits[index] );
    auto range = d_ranges[index];
    bool pass  = scale * value >= range[0] && scale * value <= range[1];
    if ( throwError && !pass ) {
        std::stringstream ss;
        ss << "Property '" + name + "' out of range in function '" + d_name + "'\n";
        ss << "Value is " << value << " ";
        ss << std::endl << "Valid range is [" << range[0] << "," << range[1] << "]" << std::endl;
        AMP_ERROR( ss.str() );
    }
    return pass;
}
template<class INPUT_VTYPE>
inline bool Property::in_range( const std::string &name,
                                const INPUT_VTYPE &values,
                                Units unit,
                                bool throwError ) const
{
    auto index = get_arg_index( name );
    if ( index == -1 )
        return true;
    double scale = 1.0;
    if ( !unit.isNull() && !d_argUnits[index].isNull() )
        scale = unit.convert( d_argUnits[index] );
    auto range = d_ranges[index];
    bool pass  = true;
    for ( auto value : values )
        pass = pass && scale * value >= range[0] && scale * value <= range[1];
    if ( throwError && !pass ) {
        std::stringstream ss;
        ss << "Property '" + name + "' out of range in function '" + d_name + "'\n";
        ss << "Values are ";
        for ( auto &value : values )
            ss << value << " ";
        ss << std::endl << "Valid range is [" << range[0] << "," << range[1] << "]" << std::endl;
        AMP_ERROR( ss.str() );
    }
    return pass;
}


/************************************************************************
 *  Get/Set auxillary data                                               *
 ************************************************************************/
template<class TYPE>
TYPE Property::getAuxiliaryData( const std::string &key ) const
{
    return d_auxiliaryData.getScalar<TYPE>( key );
}
template<class TYPE>
void Property::setAuxiliaryData( const std::string &key, const TYPE &data )
{
    return d_auxiliaryData.putScalar( key, data );
}


} // namespace AMP::Materials
