#include <algorithm>
#include <sstream>

#include "AMP/vectors/Vector.h"
#include "AMP/vectors/data/VectorDataIterator.h"

namespace AMP::Materials {


/************************************************************************
 *  Determine if a set of values are all within range or not             *
 ************************************************************************/
inline bool Property::in_range( const std::string &name, double value, Units unit, bool throwError )
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
                                bool throwError )
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


} // namespace AMP::Materials
