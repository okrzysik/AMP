#include "AMP/utils/Units.h"
#include "AMP/utils/Utilities.h"


namespace AMP {


constexpr double Units::d_pow10[22];


/********************************************************************
 * Write a string for the units                                      *
 ********************************************************************/
std::string Units::str() const
{
    if ( d_unit[0] != 0 ) {
        return printUnit();
    } else {
        return printSI();
    }
}
std::string Units::printUnit() const { return std::string( d_unit.data(), 0, d_unit.size() ); }
AMP::string_view Units::getPrefixStr( UnitPrefix p ) noexcept
{
    static constexpr const char *d_prefixSymbol[] = { "y", "z",  "a",  "f", "p", "n", "u", "m",
                                                      "c", "da", "\0", "d", "h", "k", "M", "G",
                                                      "T", "P",  "E",  "Z", "Y", "u" };
    auto i                                        = static_cast<int8_t>( p );
    return d_prefixSymbol[i];
}
std::string Units::printSI() const
{
    std::string s = Utilities::stringf( "%e", d_scale );
    for ( size_t i = 0; i < d_SI.size(); i++ ) {
        if ( d_SI[i] != 0 ) {
            s += " ";
            s += d_SI_units[i];
            if ( d_SI[i] != 1 )
                s += "^" + std::to_string( d_SI[i] );
        }
    }
    return s;
}
std::string Units::printFull() const
{
    return std::string( d_unit.data(), 0, d_unit.size() ) + " - " + printSI();
}


/********************************************************************
 * Run some compile-time tests                                       *
 ********************************************************************/
static_assert( sizeof( Units ) == 48 );
static_assert( std::is_final<Units>::value );
static_assert( std::is_trivially_copyable<Units>::value );
static_assert( !std::is_arithmetic<Units>::value );
static_assert( Units::convert( Units::getUnitPrefix( "y" ) ) == 1e-24 );
static_assert( Units::convert( Units::getUnitPrefix( "z" ) ) == 1e-21 );
static_assert( Units::convert( Units::getUnitPrefix( "a" ) ) == 1e-18 );
static_assert( Units::convert( Units::getUnitPrefix( "f" ) ) == 1e-15 );
static_assert( Units::convert( Units::getUnitPrefix( "p" ) ) == 1e-12 );
static_assert( Units::convert( Units::getUnitPrefix( "n" ) ) == 1e-9 );
static_assert( Units::convert( Units::getUnitPrefix( "u" ) ) == 1e-6 );
static_assert( Units::convert( Units::getUnitPrefix( "m" ) ) == 1e-3 );
static_assert( Units::convert( Units::getUnitPrefix( "c" ) ) == 1e-2 );
static_assert( Units::convert( Units::getUnitPrefix( "d" ) ) == 0.1 );
static_assert( Units::convert( Units::getUnitPrefix( "" ) ) == 1 );
static_assert( Units::convert( Units::getUnitPrefix( "da" ) ) == 10 );
static_assert( Units::convert( Units::getUnitPrefix( "h" ) ) == 100 );
static_assert( Units::convert( Units::getUnitPrefix( "k" ) ) == 1e3 );
static_assert( Units::convert( Units::getUnitPrefix( "M" ) ) == 1e6 );
static_assert( Units::convert( Units::getUnitPrefix( "G" ) ) == 1e9 );
static_assert( Units::convert( Units::getUnitPrefix( "T" ) ) == 1e12 );
static_assert( Units::convert( Units::getUnitPrefix( "P" ) ) == 1e15 );
static_assert( Units::convert( Units::getUnitPrefix( "E" ) ) == 1e18 );
static_assert( Units::convert( Units::getUnitPrefix( "Z" ) ) == 1e21 );
static_assert( Units::convert( Units::getUnitPrefix( "Y" ) ) == 1e24 );
static_assert( Units( "ergs/(s*cm^2)" ).convert( Units( "W/(m^2)" ) ) == 1e-3 );
static_assert( Units( "meter" ).getType() == UnitType::length );
static_assert( Units( "gram" ).getType() == UnitType::mass );
static_assert( Units( "second" ).getType() == UnitType::time );
static_assert( Units( "ampere" ).getType() == UnitType::current );
static_assert( Units( "kelvin" ).getType() == UnitType::temperature );
static_assert( Units( "joule" ).getType() == UnitType::energy );
static_assert( Units( "erg" ).getType() == UnitType::energy );
static_assert( Units( "watt" ).getType() == UnitType::power );
static_assert( Units( "mole" ).getType() == UnitType::mole );
static_assert( Units( "candela" ).getType() == UnitType::intensity );
static_assert( Units( "degree" ).getType() == UnitType::angle );
static_assert( Units( "radian" ).getType() == UnitType::angle );
static_assert( Units( "V" ) * Units( "A" ) == Units( "W" ) );
static_assert( Units( "W/m^2" ) == Units( "uW/mm^2" ) );


} // namespace AMP
