#ifndef included_AMP_Units_hpp
#define included_AMP_Units_hpp

#include "AMP/utils/Units.h"


namespace AMP {


/********************************************************************
 * Constructors                                                      *
 ********************************************************************/
constexpr Units::Units() : d_prefix( UnitPrefix::unknown ), d_unit( UnitValue::unknown ) {}
constexpr Units::Units( UnitPrefix p, UnitValue u ) : d_prefix( p ), d_unit( u ) {}
constexpr Units::Units( const AMP::string_view &str )
    : d_prefix( UnitPrefix::unknown ), d_unit( UnitValue::unknown )
{
    // Remove trailing/preceeding whitespace
    int i1 = 0, i2 = str.size() - 1;
    for ( ; i1 < (int) str.size() && str[i1] == ' '; i1++ ) {}
    for ( ; i2 > 0 && str[i2] == ' '; i2-- ) {}
    auto unit = i1 <= i2 ? str.substr( i1, i2 - i1 + 1 ) : AMP::string_view();
    // Check if the character '-' is present indicating a seperation between the prefix and unit
    size_t index = unit.find( '-' );
    if ( index != std::string::npos ) {
        d_prefix = getUnitPrefix( unit.substr( 0, index ) );
        d_unit   = getUnitValue( unit.substr( index + 1 ) );
    } else {
        if ( unit.size() <= 1 ) {
            d_prefix = UnitPrefix::none;
            d_unit   = getUnitValue( unit );
        } else if ( unit.substr( 0, 2 ) == "da" ) {
            d_prefix = UnitPrefix::deca;
            d_unit   = getUnitValue( unit.substr( 2 ) );
        } else {
            d_prefix = getUnitPrefix( unit.substr( 0, 1 ) );
            d_unit   = getUnitValue( unit.substr( 1 ) );
            if ( d_prefix == UnitPrefix::unknown || d_unit == UnitValue::unknown ) {
                d_prefix = UnitPrefix::none;
                d_unit   = getUnitValue( unit );
            }
        }
    }
}


/********************************************************************
 * Get prefix                                                        *
 ********************************************************************/
constexpr Units::UnitPrefix Units::getUnitPrefix( const AMP::string_view &str ) noexcept
{
    Units::UnitPrefix value = UnitPrefix::unknown;
    if ( str.empty() ) {
        value = UnitPrefix::none;
    } else if ( str == "yotta" || str == "Y" ) {
        value = UnitPrefix::yotta;
    } else if ( str == "zetta" || str == "Z" ) {
        value = UnitPrefix::zetta;
    } else if ( str == "exa" || str == "E" ) {
        value = UnitPrefix::exa;
    } else if ( str == "peta" || str == "P" ) {
        value = UnitPrefix::peta;
    } else if ( str == "tera" || str == "T" ) {
        value = UnitPrefix::tera;
    } else if ( str == "giga" || str == "G" ) {
        value = UnitPrefix::giga;
    } else if ( str == "mega" || str == "M" ) {
        value = UnitPrefix::mega;
    } else if ( str == "kilo" || str == "k" ) {
        value = UnitPrefix::kilo;
    } else if ( str == "hecto" || str == "h" ) {
        value = UnitPrefix::hecto;
    } else if ( str == "deca" || str == "da" ) {
        value = UnitPrefix::deca;
    } else if ( str == "deci" || str == "d" ) {
        value = UnitPrefix::deci;
    } else if ( str == "centi" || str == "c" ) {
        value = UnitPrefix::centi;
    } else if ( str == "milli" || str == "m" ) {
        value = UnitPrefix::milli;
    } else if ( str == "micro" || str == "u" ) {
        value = UnitPrefix::micro;
    } else if ( str == "nano" || str == "n" ) {
        value = UnitPrefix::nano;
    } else if ( str == "pico" || str == "p" ) {
        value = UnitPrefix::pico;
    } else if ( str == "femto" || str == "f" ) {
        value = UnitPrefix::femto;
    } else if ( str == "atto" || str == "a" ) {
        value = UnitPrefix::atto;
    } else if ( str == "zepto" || str == "z" ) {
        value = UnitPrefix::zepto;
    } else if ( str == "yocto" || str == "y" ) {
        value = UnitPrefix::yocto;
    }
    return value;
}


/********************************************************************
 * Get unit value                                                    *
 ********************************************************************/
constexpr Units::UnitValue Units::getUnitValue( const AMP::string_view &str ) noexcept
{
    Units::UnitValue value = UnitValue::unknown;
    if ( str == "meter" || str == "m" ) {
        value = UnitValue::meter;
    } else if ( str == "gram" || str == "g" ) {
        value = UnitValue::gram;
    } else if ( str == "second" || str == "s" ) {
        value = UnitValue::second;
    } else if ( str == "ampere" || str == "A" ) {
        value = UnitValue::ampere;
    } else if ( str == "kelvin" || str == "K" ) {
        value = UnitValue::kelvin;
    } else if ( str == "joule" || str == "J" ) {
        value = UnitValue::joule;
    } else if ( str == "ergs" || str == "erg" ) {
        value = UnitValue::erg;
    } else if ( str == "degree" || str == "degrees" ) {
        value = UnitValue::degree;
    } else if ( str == "radian" || str == "radians" ) {
        value = UnitValue::radian;
    }
    return value;
}


/********************************************************************
 * Get unit type                                                     *
 ********************************************************************/
constexpr Units::UnitType Units::getUnitType( UnitValue u ) noexcept
{
    switch ( u ) {
    case UnitValue::meter:
        return UnitType::length;
    case UnitValue::gram:
        return UnitType::mass;
    case UnitValue::second:
        return UnitType::time;
    case UnitValue::ampere:
        return UnitType::current;
    case UnitValue::kelvin:
        return UnitType::temperature;
    case UnitValue::joule:
    case UnitValue::erg:
        return UnitType::energy;
    case UnitValue::degree:
    case UnitValue::radian:
        return UnitType::angle;
    default:
        return UnitType::unknown;
    }
}


/********************************************************************
 * Convert to another unit system                                    *
 ********************************************************************/
constexpr double Units::convert( const Units &rhs ) const noexcept
{
    if ( this->operator==( rhs ) )
        return 1;
    // Convert the prefix
    double cp = convert( d_prefix ) / convert( rhs.d_prefix );
    if ( d_unit == rhs.d_unit )
        return cp; // Only need to convert prefix
    // Convert the unit
    if ( getUnitType( d_unit ) != getUnitType( rhs.d_unit ) )
        return 0; // Invalid conversion
    double cu = 0;
    if ( d_unit == UnitValue::joule && rhs.d_unit == UnitValue::erg )
        cu = 1e7;
    else if ( d_unit == UnitValue::erg && rhs.d_unit == UnitValue::joule )
        cu = 1e-7;
    else if ( d_unit == UnitValue::degree && rhs.d_unit == UnitValue::radian )
        cu = 0.017453292519943;
    else if ( d_unit == UnitValue::radian && rhs.d_unit == UnitValue::degree )
        cu = 57.295779513082323;
    // Return the total conversion
    return cp * cu;
}


} // namespace AMP


#endif
