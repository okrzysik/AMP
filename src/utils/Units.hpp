#ifndef included_AMP_Units_hpp
#define included_AMP_Units_hpp

#include "AMP/utils/Units.h"


namespace AMP {


/********************************************************************
 * Constructors                                                      *
 ********************************************************************/
constexpr Units::Units() : d_unit( { 0 } ), d_SI( { 0 } ), d_scale( 0 ) {}
constexpr Units::Units( const AMP::string_view &str ) : d_unit( { 0 } ), d_SI( { 0 } ), d_scale( 0 )
{
    // Remove trailing/preceeding whitespace
    int i1 = 0, i2 = str.size() - 1;
    for ( ; i1 < (int) str.size() && ( str[i1] == ' ' || str[i1] == '"' ); i1++ ) {}
    for ( ; i2 > 0 && ( str[i2] == ' ' || str[i2] == '"' ); i2-- ) {}
    // Copy the unit to internal memory (removed extra spaces)
    if ( i1 > i2 )
        return;
    int k     = 0;
    char last = ';';
    for ( int i = i1; i <= i2; i++ ) {
        if ( str[i] == '*' || str[i] == ' ' ) {
            if ( last != '*' && last != '/' && last != '(' && last != ')' ) {
                d_unit[k] = '*';
                k++;
            }
        } else if ( str[i] == '/' || str[i] == '(' || str[i] == ')' ) {
            if ( last == '*' ) {
                d_unit[k - 1] = str[i];
            } else {
                d_unit[k] = str[i];
                k++;
            }
        } else {
            d_unit[k] = str[i];
            k++;
        }
        last = d_unit[k - 1];
        if ( k >= (int) d_unit.size() )
            throw std::logic_error( "Unit size" );
    }
    AMP::string_view unit( d_unit.data() );
    // Convert the string to SI units
    k                    = 0;
    d_scale              = 1.0;
    auto findMatchingPar = [&unit]( size_t i ) {
        for ( size_t j = i, count = 1; j < unit.size(); j++ ) {
            if ( unit[j] == '(' )
                count++;
            else if ( unit[j] == ')' )
                count--;
            if ( count == 0 )
                return j;
        }
        return unit.size();
    };
    auto findSpace = [&unit]( size_t i ) {
        for ( size_t j = i + 1; j < unit.size(); j++ ) {
            if ( unit[j] == '*' || unit[j] == '/' || unit[j] == '(' )
                return j;
        }
        return unit.size();
    };
    for ( size_t i = 0; i < unit.size(); ) {
        if ( unit[i] == '*' )
            i++;
        size_t j  = 0;
        double s  = 0;
        SI_type u = { 0 };
        if ( unit[i] == '(' ) {
            j = findMatchingPar( i );
            Units u2( unit.substr( i + 1, j - i - 2 ) );
            s = u2.d_scale;
            u = u2.d_SI;
        } else if ( unit[i] == '/' ) {
            if ( unit[i + 1] == '(' ) {
                j = findMatchingPar( i + 1 );
                Units u2( unit.substr( i + 2, j - i - 3 ) );
                s = 1.0 / u2.d_scale;
                u = u2.d_SI;
                for ( size_t i = 0; i < u.size(); i++ )
                    u[i] = -u[i];
            } else {
                j             = findSpace( i );
                auto [u2, s2] = read( unit.substr( i + 1, j - i - 1 ) );
                s             = 1.0 / s2;
                for ( size_t i = 0; i < u.size(); i++ )
                    u[i] = -u2[i];
            }
        } else {
            j             = findSpace( i );
            auto [u2, s2] = read( unit.substr( i, j - i ) );
            u             = u2;
            s             = s2;
        }
        d_scale *= s;
        d_SI = combine( d_SI, u );
        i    = j;
    }
}
constexpr std::tuple<Units::SI_type, double> Units::read( const AMP::string_view &str )
{
    auto i = str.find( '^' );
    if ( i == std::string::npos ) {
        // Check if the character '-' is present indicating a seperation between the prefix and unit
        size_t index = str.find( '-' );
        if ( index != std::string::npos ) {
            auto prefix = getUnitPrefix( str.substr( 0, index ) );
            auto [u, s] = read2( str.substr( index + 1 ) );
            s *= convert( prefix );
            return std::tie( u, s );
        } else {
            return read2( str );
        }
        throw std::logic_error( "Error reading" );
    } else {
        auto [u, s] = read2( str.substr( 0, i ) );
        auto tmp    = str.substr( i + 1 );
        if ( tmp[0] == '-' ) {
            for ( auto &u2 : u )
                u2 = -u2;
            s   = 1.0 / s;
            tmp = tmp.substr( 1 );
        }
        if ( tmp.size() != 1u )
            throw std::logic_error( "Conversion error (1)" );
        int p = tmp[0] - 48;
        if ( p <= 0 || p > 3 )
            throw std::logic_error( "Conversion error (2)" );
        double s0 = s;
        for ( int j = 1; j < p; j++ )
            s *= s0;
        for ( size_t i = 0; i < u.size(); i++ )
            u[i] *= p;
        return std::tie( u, s );
    }
}
constexpr Units::SI_type Units::combine( const SI_type &a, const SI_type &b )
{
    SI_type c = { 0 };
    for ( size_t i = 0; i < a.size(); i++ )
        c[i] = a[i] + b[i];
    return c;
}


/********************************************************************
 * Get prefix                                                        *
 ********************************************************************/
constexpr UnitPrefix Units::getUnitPrefix( const AMP::string_view &str ) noexcept
{
    UnitPrefix value = UnitPrefix::unknown;
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
constexpr std::tuple<Units::SI_type, double> Units::read2( const AMP::string_view &str )
{
    if ( str.size() <= 1 ) {
        auto [u, s] = readUnit( str );
        return std::tie( u, s );
    } else if ( str.substr( 0, 2 ) == "da" ) {
        auto prefix = UnitPrefix::deca;
        auto [u, s] = readUnit( str.substr( 2 ) );
        s *= convert( prefix );
        return std::tie( u, s );
    } else {
        auto prefix = getUnitPrefix( str.substr( 0, 1 ) );
        auto [u, s] = readUnit( str.substr( 1 ), false );
        if ( prefix == UnitPrefix::unknown || s == 0 ) {
            return readUnit( str );
        } else {
            s *= convert( prefix );
            return std::tie( u, s );
        }
    }
    SI_type u = { 0 };
    double s  = 0.0;
    return std::tie( u, s );
}
constexpr std::tuple<Units::SI_type, double> Units::readUnit( const AMP::string_view &str,
                                                              bool throwErr )
{
    auto create = []( UnitType type, double s = 1.0 ) {
        auto u = getSI( type );
        return std::tuple<Units::SI_type, double>( u, s );
    };
    // Check base SI units
    if ( str == "second" || str == "s" )
        return create( UnitType::time );
    if ( str == "meter" || str == "m" )
        return create( UnitType::length );
    if ( str == "gram" || str == "g" )
        return create( UnitType::mass, 1e-3 );
    if ( str == "ampere" || str == "A" )
        return create( UnitType::current );
    if ( str == "kelvin" || str == "K" )
        return create( UnitType::temperature );
    if ( str == "mole" || str == "mol" )
        return create( UnitType::mole );
    if ( str == "candela" || str == "cd" )
        return create( UnitType::intensity );
    if ( str == "radian" || str == "radians" || str == "rad" )
        return create( UnitType::angle );
    if ( str == "steradian" || str == "sr" )
        return create( UnitType::solidAngle );
    // Check derived SI units
    if ( str == "degree" || str == "degrees" )
        return create( UnitType::angle, 0.017453292519943 );
    if ( str == "joule" || str == "J" )
        return create( UnitType::energy );
    if ( str == "watt" || str == "W" )
        return create( UnitType::power );
    if ( str == "hertz" || str == "Hz" )
        return create( UnitType::frequency );
    if ( str == "newton" || str == "N" )
        return create( UnitType::force );
    if ( str == "pascal" || str == "Pa" )
        return create( UnitType::pressure );
    if ( str == "coulomb" || str == "C" )
        return create( UnitType::electricCharge );
    if ( str == "volt" || str == "V" )
        return create( UnitType::electricalPotential );
    if ( str == "farad" || str == "F" )
        return create( UnitType::capacitance );
    if ( str == "ohm" )
        return create( UnitType::resistance );
    if ( str == "siemens" || str == "S" )
        return create( UnitType::electricalConductance );
    if ( str == "weber" || str == "Wb" )
        return create( UnitType::magneticFlux );
    if ( str == "tesla" || str == "T" )
        return create( UnitType::magneticFluxDensity );
    if ( str == "henry" || str == "H" )
        return create( UnitType::inductance );
    if ( str == "lumen" || str == "lm" )
        return create( UnitType::luminousFlux );
    if ( str == "lux" || str == "lx" )
        return create( UnitType::illuminance );
    if ( str == "litre" || str == "L" )
        return std::tuple<Units::SI_type, double>( { 0, 3, 0, 0, 0, 0, 0, 0, 0 }, 1e-3 );
    // Check cgs units
    if ( str == "ergs" || str == "erg" )
        return create( UnitType::energy, 1e-7 );
    if ( str == "eV" )
        return create( UnitType::energy, 1.602176634e-19 );
    // Check English units
    if ( str == "inch" || str == "in" || str == "\"" )
        return create( UnitType::length, 0.0254 );
    if ( str == "foot" || str == "ft" || str == "\'" )
        return create( UnitType::length, 0.3048 );
    if ( str == "yard" || str == "yd" )
        return create( UnitType::length, 0.9144 );
    if ( str == "furlong" || str == "fur" )
        return create( UnitType::length, 201.168 );
    if ( str == "mile" || str == "mi" )
        return create( UnitType::length, 1609.344 );
    if ( str == "acre" )
        return std::tuple<Units::SI_type, double>( { 0, 2, 0, 0, 0, 0, 0, 0, 0 }, 4046.8564224 );
    if ( str == "pint" || str == "pt" )
        return std::tuple<Units::SI_type, double>( { 0, 3, 0, 0, 0, 0, 0, 0, 0 }, 0.56826125 );
    if ( str == "quart" || str == "qt" )
        return std::tuple<Units::SI_type, double>( { 0, 3, 0, 0, 0, 0, 0, 0, 0 }, 1.1365225 );
    if ( str == "gallon" || str == "gal" )
        return std::tuple<Units::SI_type, double>( { 0, 3, 0, 0, 0, 0, 0, 0, 0 }, 4.54609 );
    if ( str == "ounce" || str == "oz" )
        return create( UnitType::mass, 0.028349523125 );
    if ( str == "pound" || str == "lb" )
        return create( UnitType::mass, 0.45359237 );
    if ( str == "ton" )
        return create( UnitType::mass, 1016.0469088 );
    // Check atomic units
    if ( str == "hartree" )
        return create( UnitType::energy, 4.359744722207185e-18 );
    if ( str == "bohr" )
        return create( UnitType::length, 5.2917721090380e-11 );
    // No success
    SI_type u = { 0 };
    double s  = 0;
    if ( throwErr && s == 0 ) {
        // We have an error intpreting the string, try to interpret it as a double
        // Note: this is only valid when we are not using constexpr (also true for throw)
        const char *begin = str.begin();
        const char *end   = str.end();
        s                 = strtod( (char *) begin, (char **) &end );
        if ( s != 0 )
            return std::tie( u, s );
        // Still unable, throw an error
        char err[64] = "Unknown unit: ";
        for ( size_t i = 0; i < str.size(); i++ )
            err[14 + i] = str[i];
        throw std::logic_error( err );
    }
    return std::tie( u, s );
}


/********************************************************************
 * Get unit type                                                     *
 ********************************************************************/
constexpr Units::SI_type Units::getSI( UnitType type )
{
    // s, m, kg, A, K, mol, cd, rad, sr
    constexpr SI_type id[25] = {
        { 0, 0, 0, 0, 0, 0, 0, 0, 0 },    // 0: unknown
        { 1, 0, 0, 0, 0, 0, 0, 0, 0 },    // 1: time (s)
        { 0, 1, 0, 0, 0, 0, 0, 0, 0 },    // 2: length (m)
        { 0, 0, 1, 0, 0, 0, 0, 0, 0 },    // 3: mass (kg)
        { 0, 0, 0, 1, 0, 0, 0, 0, 0 },    // 4: current (A)
        { 0, 0, 0, 0, 1, 0, 0, 0, 0 },    // 5: temperature (K)
        { 0, 0, 0, 0, 0, 1, 0, 0, 0 },    // 6: mole (mol)
        { 0, 0, 0, 0, 0, 0, 1, 0, 0 },    // 7: intensity (cd)
        { 0, 0, 0, 0, 0, 0, 0, 1, 0 },    // 8: angle (rad)
        { 0, 0, 0, 0, 0, 0, 0, 0, 1 },    // 9: solid angle (sr)
        { -2, 2, 1, 0, 0, 0, 0, 0 },      // 10: energy (kg m2 s-2)
        { -3, 2, 1, 0, 0, 0, 0, 0 },      // 11: power (kg m2 s-3)
        { -1, 0, 0, 0, 0, 0, 0, 0 },      // 12: frequency (s-1)
        { -2, 1, 1, 0, 0, 0, 0, 0, 0 },   // 13: force (kg m s-2)
        { -2, -1, -2, 0, 0, 0, 0, 0, 0 }, // 14: pressure (kg m-1 s-2)
        { 1, 0, 0, 1, 0, 0, 0, 0, 0 },    // 15: electric charge (s A)
        { -3, 2, 1, -1, 0, 0, 0, 0, 0 },  // 16: electrical potential (kg m2 s-3 A-1)
        { 4, -2, -1, 2, 0, 0, 0, 0, 0 },  // 17: capacitance (kg-1 m-2 s4 A2)
        { -3, 2, 1, -2, 0, 0, 0, 0, 0 },  // 18: resistance (kg m2 s-3 A-2)
        { 3, -2, -1, 2, 0, 0, 0, 0, 0 },  // 19: electrical conductance (kg-1 m-2 s3 A2)
        { -2, 2, 1, -1, 0, 0, 0, 0, 0 },  // 20: magnetic flux (kg m2 s-2 A-1)
        { -2, 0, 1, -1, 0, 0, 0, 0, 0 },  // 21: magnetic flux density (kg s-2 A-1)
        { -2, 2, 1, -2, 0, 0, 0, 0, 0 },  // 22: inductance (kg m2 s-2 A-2)
        { 0, 0, 0, 0, 0, 0, 1, 0, 1 },    // 23: luminous flux (cd sr)
        { 0, -2, 0, 0, 0, 0, 1, 0, 1 }    // 24: illuminance (cd sr m-2)
    };
    return id[static_cast<int>( type )];
}
constexpr UnitType Units::getType() const noexcept
{
    auto compare = []( const SI_type &a, const SI_type &b ) {
        bool test = true;
        for ( size_t i = 0; i < a.size(); i++ )
            test = test && a[i] == b[i];
        return test;
    };
    for ( int i = 0; i <= 24; i++ ) {
        auto type = static_cast<UnitType>( i );
        auto id   = getSI( type );
        if ( compare( id, d_SI ) )
            return type;
    }
    return UnitType::unknown;
}


/********************************************************************
 * Convert to another unit system                                    *
 ********************************************************************/
constexpr bool operator==( const std::array<int8_t, 9> &a, const std::array<int8_t, 9> &b )
{
    bool test = true;
    for ( size_t i = 0; i < a.size(); i++ )
        test = test && a[i] == b[i];
    return test;
}
constexpr bool Units::compatible( const Units &rhs ) noexcept
{
    constexpr SI_type energy     = { -2, 2, 1, 0, 0, 0, 0, 0 };
    constexpr SI_type tmperature = { 0, 0, 0, 0, 1, 0, 0, 0, 0 };
    if ( d_SI == rhs.d_SI )
        return true;
    if ( d_SI == energy && rhs.d_SI == tmperature )
        return true; // Convert from energy to temperature (using boltzmann constant)
    if ( d_SI == tmperature && rhs.d_SI == energy )
        return true; // Convert from temperature to energy (using boltzmann constant)
    return false;
}
constexpr double Units::convert( const Units &rhs ) const
{
    constexpr SI_type energy     = { -2, 2, 1, 0, 0, 0, 0, 0 };
    constexpr SI_type tmperature = { 0, 0, 0, 0, 1, 0, 0, 0, 0 };
    if ( d_SI == rhs.d_SI ) {
        // The SI units match
        return d_scale / rhs.d_scale;
    } else if ( d_SI == energy && rhs.d_SI == tmperature ) {
        // Convert from energy to temperature (using bolztmann constant)
        return 7.242971666663E+22 * d_scale / rhs.d_scale;
    } else if ( d_SI == tmperature && rhs.d_SI == energy ) {
        // Convert from temperature to energy (using bolztmann constant)
        return 1.38064878066922E-23 * d_scale / rhs.d_scale;
    } else {
        throw std::logic_error( "Incompatible units: " + printSIBase() + " - " +
                                rhs.printSIBase() );
    }
}


/********************************************************************
 * Compare two units                                                 *
 ********************************************************************/
constexpr bool Units::operator==( const Units &rhs ) const noexcept
{
    bool test = fabs( d_scale - rhs.d_scale ) < 1e-6 * d_scale;
    for ( size_t i = 0; i < d_SI.size(); i++ )
        test = test && d_SI[i] == rhs.d_SI[i];
    return test;
}


/********************************************************************
 * Multiply/Divide units                                             *
 ********************************************************************/
constexpr void Units::operator*=( const Units &rhs ) noexcept
{
    d_unit = { 0 };
    d_SI   = combine( d_SI, rhs.d_SI );
    d_scale *= rhs.d_scale;
}
constexpr void Units::operator/=( const Units &rhs ) noexcept
{
    SI_type tmp = { 0 };
    for ( size_t i = 0; i < tmp.size(); i++ )
        tmp[i] = -rhs.d_SI[i];
    d_unit = { 0 };
    d_SI   = combine( d_SI, tmp );
    d_scale *= rhs.d_scale;
}
constexpr Units operator*( const Units &a, const Units &b )
{
    Units c = a;
    c *= b;
    return c;
}
constexpr Units operator/( const Units &a, const Units &b )
{
    Units c = a;
    c /= b;
    return c;
}


} // namespace AMP

#endif
