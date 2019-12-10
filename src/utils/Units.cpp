#include "AMP/utils/Units.h"
#include "AMP/utils/Utilities.h"


namespace AMP {


constexpr double Units::d_pow10[22];
constexpr char Units::d_prefixSymbol[];

static_assert( sizeof( Units ) == 2, "Unexpected size for Units" );


/********************************************************************
 * Write a string for the units                                      *
 ********************************************************************/
std::string Units::str() const
{
    AMP_ASSERT( !isNull() );
    std::string name( str( d_prefix ).data() );
    auto unit = str( d_unit );
    name.append( unit.data(), unit.size() );
    return name;
}
std::array<char, 3> Units::str( UnitPrefix p ) noexcept
{
    std::array<char, 3> str;
    str[0] = d_prefixSymbol[static_cast<int8_t>( p )];
    str[1] = 0;
    str[2] = 0;
    if ( p == UnitPrefix::deca )
        str[1] = 'a';
    return str;
}
AMP::string_view Units::str( UnitValue u ) noexcept
{
    if ( u == UnitValue::meter ) {
        return "m";
    } else if ( u == UnitValue::gram ) {
        return "g";
    } else if ( u == UnitValue::second ) {
        return "s";
    } else if ( u == UnitValue::ampere ) {
        return "A";
    } else if ( u == UnitValue::kelvin ) {
        return "K";
    } else if ( u == UnitValue::joule ) {
        return "J";
    } else if ( u == UnitValue::erg ) {
        return "erg";
    } else if ( u == UnitValue::degree ) {
        return "degree";
    } else if ( u == UnitValue::radian ) {
        return "radian";
    }
    return "unknown";
}


} // namespace AMP
