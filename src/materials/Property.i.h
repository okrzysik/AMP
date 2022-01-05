#include <algorithm>
#include <sstream>

#include "AMP/vectors/Vector.h"
#include "AMP/vectors/data/VectorDataIterator.h"

namespace AMP::Materials {


/************************************************************************
 *  Determine if a set of values are all within range or not             *
 ************************************************************************/
inline bool Property::in_range( const std::string &argname, const double value )
{
    if ( !is_argument( argname ) )
        return true;
    auto range = get_arg_range( argname );
    return value >= range[0] && value <= range[1];
}
template<class INPUT_VTYPE>
inline bool Property::in_range( const std::string &argname, const INPUT_VTYPE &values )
{
    if ( !is_argument( argname ) )
        return true;
    auto range  = get_arg_range( argname );
    bool result = true;
    auto pos    = values.begin();
    auto end    = values.end();
    while ( pos != end ) {
        result = result && *pos >= range[0] && *pos <= range[1];
        ++pos;
    }
    return result;
}
template<class INPUT_VTYPE>
inline bool Property::in_range( const std::map<std::string, std::shared_ptr<INPUT_VTYPE>> &values )
{
    bool result = true;
    for ( const auto &value : values ) {
        result = result && in_range( value.first, *( value.second ) );
    }
    return result;
}


} // namespace AMP::Materials
