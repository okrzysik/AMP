#include <algorithm>
#include <sstream>

#include "AMP/vectors/Vector.h"
#include "AMP/vectors/data/VectorDataIterator.h"

namespace AMP::Materials {


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
 *  evalv                                                                *
 ************************************************************************/
template<class VEC, class... Args>
void Property::convertArgs1( std::vector<argumentDataStruct<VEC>> &args2,
                             const std::string &name,
                             const VEC &v,
                             Args... args )
{
    constexpr size_t n = sizeof...( args );
    if constexpr ( n > 0 ) {
        auto t = std::make_tuple( std::forward<Args>( args )... );
        if constexpr ( std::is_same<Units, decltype( std::get<0>( t ) )>::value ) {
            convertArgs2( args2, name, args... );
            return;
        }
    }
    args2.emplace_back( name, v );
    if constexpr ( n > 0 )
        convertArgs1( args2, args... );
}
template<class VEC, class... Args>
void Property::convertArgs2( std::vector<argumentDataStruct<VEC>> &args2,
                             const std::string &name,
                             const VEC &v,
                             const Units &u,
                             Args... args )
{
    constexpr size_t n = sizeof...( args );
    args2.emplace_back( name, v, u );
    if constexpr ( n > 0 )
        convertArgs1( args2, args... );
}
template<class VEC, class... Args>
std::vector<Property::argumentDataStruct<VEC>> Property::convertArgs( Args... args )
{
    constexpr size_t n = sizeof...( args );
    if ( n == 0 )
        return {};
    std::vector<argumentDataStruct<VEC>> args2;
    convertArgs1( args2, args... );
    return args2;
}
/*template<class... Args>
void Property::evalv( std::vector<double> &r, Units u, Args... args ) const
{
    auto args2 = convertArgs<std::vector<double>>( args... );
    evalv( r, u, args2 );
}*/


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
