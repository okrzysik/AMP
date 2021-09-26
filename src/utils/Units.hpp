#ifndef included_AMP_Units_hpp
#define included_AMP_Units_hpp

#include "AMP/utils/Units.h"

#include <iostream>
#include <tuple>


namespace AMP {


/********************************************************************
 * Helper functions                                                  *
 ********************************************************************/
constexpr std::string_view deblank( const std::string_view &str )
{
    if ( str.empty() )
        return std::string_view();
    int i1 = 0, i2 = str.size() - 1;
    for ( ; i1 < (int) str.size() && ( str[i1] == ' ' || str[i1] == '\t' ); i1++ ) {}
    for ( ; i2 > 0 && ( str[i2] == ' ' || str[i2] == '\t' || str[i2] == '\r' ); i2-- ) {}
    if ( i2 == 0 && ( str[i2] == ' ' || str[i2] == '\t' || str[i2] == '\r' ) )
        return std::string_view();
    return str.substr( i1, i2 - i1 + 1 );
}
// Find next token
constexpr std::pair<size_t, char> Units::findToken( const std::string_view &str, size_t i )
{
    char op = 0;
    while ( i < str.size() && op == 0 ) {
        if ( str[i] == '*' || str[i] == '(' || str[i] == '^' || str[i] == '/' || str[i] == ')' ) {
            op = str[i];
        } else if ( str[i] == ' ' ) {
            op = '*';
            while ( str[i] == ' ' ) {
                if ( i == str.size() )
                    return std::make_pair( i, op );
                i++;
            }
            if ( str[i] == '*' || str[i] == '(' || ( str[i] == 'x' && str[i + 1] == ' ' ) ) {
                op = '*';
            } else if ( str[i] == '^' || str[i] == '/' || str[i] == ')' ) {
                op = str[i];
            } else {
                i--;
            }
        } else {
            i++;
        }
    }
    return std::make_pair( i, op );
}
// Find matching parenthesis and evaluate substring
constexpr size_t Units::findPar( const std::string_view &str, size_t i )
{
    for ( size_t j = i + 1, c = 1; j < str.size(); j++ ) {
        if ( str[j] == '(' )
            c++;
        if ( str[j] == ')' )
            c--;
        if ( c == 0 )
            return j;
    }
    return str.size();
}


/********************************************************************
 * Constructors                                                      *
 ********************************************************************/
constexpr Units::Units() : d_unit( { 0 } ), d_SI( { 0 } ), d_scale( 0.0 ) {}
constexpr Units::Units( const char *str ) : Units( std::string_view( str ) ) {}
constexpr Units::Units( const std::string_view &str )
    : d_unit( { 0 } ), d_SI( { 0 } ), d_scale( 0.0 )
{
    if ( !str.empty() ) {
        Units tmp = read( str );
        d_SI      = tmp.d_SI;
        d_scale   = tmp.d_scale;
        if ( str.length() < d_unit.size() - 1 ) {
            for ( size_t i = 0; i < str.length(); i++ )
                d_unit[i] = str[i];
        }
    }
}
constexpr Units::Units( const SI_type &SI, double s ) : d_unit( { 0 } ), d_SI( SI ), d_scale( s ) {}
constexpr std::array<int8_t, 9> operator+( const std::array<int8_t, 9> &a,
                                           const std::array<int8_t, 9> &b )
{
    std::array<int8_t, 9> c = { 0 };
    for ( size_t i = 0; i < a.size(); i++ )
        c[i] = a[i] + b[i];
    return c;
}
constexpr int Units::atoi( std::string_view str, bool throw_error )
{
    str = deblank( str );
    if ( str.empty() )
        return 0;
    bool neg = false;
    if ( str[0] == '-' ) {
        neg = true;
        str = deblank( str.substr( 1 ) );
    } else if ( str[0] == '+' ) {
        str = deblank( str.substr( 1 ) );
    }
    int i = 0;
    for ( size_t j = 0; j < str.size(); j++ ) {
        if ( str[j] < '0' || str[j] > '9' ) {
            if ( throw_error )
                throw std::logic_error( "Error calling atoi: " + std::string( str ) );
            else
                return std::numeric_limits<int>::min();
        }
        i = i * 10 + static_cast<int>( str[j] - '0' );
    }
    return neg ? -i : i;
}
constexpr double Units::strtod( std::string_view str, bool throw_error )
{
    str = deblank( str );
    if ( str.empty() )
        return 0;
    bool neg = false;
    if ( str[0] == '-' ) {
        neg = true;
        str = deblank( str.substr( 1 ) );
    }
    size_t i1 = str.find( '.' );
    size_t i2 = std::min( str.find( 'e' ), str.find( 'E' ) );
    auto s1   = str.substr( 0, std::min( i1, i2 ) );
    std::string_view s2, s3;
    if ( i1 != std::string::npos )
        s2 = str.substr( i1 + 1, std::min( i2 - i1 - 1, str.size() ) );
    if ( i2 != std::string::npos )
        s3 = str.substr( i2 + 1 );
    constexpr int errInt = std::numeric_limits<int>::min();
    double x             = atoi( s1, throw_error );
    if ( x == errInt )
        return std::numeric_limits<double>::quiet_NaN();
    if ( !s2.empty() ) {
        double f = 0.1;
        for ( size_t i = 0; i < s2.size(); i++, f *= 0.1 ) {
            if ( s2[i] < '0' || s2[i] > '9' ) {
                if ( throw_error )
                    throw std::logic_error( "Error calling strtod" );
                else
                    return std::numeric_limits<double>::quiet_NaN();
            }
            x += f * static_cast<int>( s2[i] - '0' );
        }
    }
    if ( !s3.empty() ) {
        int p = atoi( s3, throw_error );
        if ( x == errInt )
            return std::numeric_limits<double>::quiet_NaN();
        double s = 10.0;
        if ( p < 0 ) {
            s = 0.1;
            p = -p;
        }
        for ( int i = 0; i < p; i++ )
            x *= s;
    }
    return neg ? -x : x;
}
constexpr Units Units::read( std::string_view str )
{
    str = deblank( str );
    if ( str.empty() )
        return Units( { 0 }, 1.0 );
    // Break the string into value/operator sets
    int N        = 0;
    char op[100] = { 0 };
    std::string_view v[100];
    for ( size_t i = 0; i < str.size(); ) {
        auto tmp  = findToken( str, i );
        size_t i2 = tmp.first;
        op[N]     = tmp.second;
        if ( tmp.second == '(' ) {
            bool test = false;
            for ( size_t j = i; j < i2; j++ )
                test = test || str[j];
            if ( test ) {
                op[N]  = '*';
                v[N++] = str.substr( i, i2 - i );
                i      = i2;
            } else {
                op[N]     = 0;
                i2        = findPar( str, i2 );
                size_t i3 = i2;
                for ( ; i3 < str.size() && str[i3] == ' '; i3++ ) {}
                if ( i3 + 1 < str.size() ) {
                    op[N] = '*';
                    if ( str[i3 + 1] == '^' || str[i3 + 1] == '*' || str[i3 + 1] == '/' ) {
                        i3++;
                        op[N] = str[i3];
                    }
                }
                v[N++] = str.substr( i + 1, i2 - i - 1 );
                i      = i3 + 1;
            }
        } else {
            v[N++] = str.substr( i, i2 - i );
            i      = i2 + 1;
        }
    }
    // Evaluate if we only have one entry
    if ( N == 1 )
        return read2( v[0] );
    // Evaluate and apply operators
    Units u( { 0 }, 1.0 );
    char last_op = '*';
    for ( int i = 0; i < N; i++ ) {
        Units u2 = read( v[i] );
        if ( op[i] == '^' ) {
            u2 = u2.pow( atoi( v[i + 1] ) );
            i++;
        }
        if ( last_op == '*' ) {
            u *= u2;
        } else if ( last_op == '/' ) {
            u /= u2;
        } else {
            throw std::logic_error( "Invalid operator" );
        }
        last_op = op[i];
    }
    return u;
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
constexpr UnitPrefix Units::getUnitPrefix( const std::string_view &str ) noexcept
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
    } else if ( str == "micro" || str == "u" || str == "μ" ) {
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
inline std::vector<std::string> Units::getAllPrefixes()
{
    return { "Y", "Z", "E", "P", "T", "G", "M", "k", "h", "da", "",
             "d", "c", "m", "u", "μ", "n", "p", "f", "a", "z",  "y" };
}


/********************************************************************
 * Get unit value                                                    *
 ********************************************************************/
constexpr Units Units::read2( std::string_view str )
{
    // Check for special prefixes
    if ( str.substr( 0, 2 ) == "da" ) {
        Units u = readUnit( str.substr( 2 ), false );
        u.d_scale *= 10.0;
        if ( u.d_scale != 0 )
            return u;
    }
    if ( str.substr( 0, 2 ) == "μ" ) {
        Units u = readUnit( str.substr( 2 ), false );
        u.d_scale *= 1e-6;
        if ( u.d_scale != 0 )
            return u;
    }
    // Try reading a prefix and then the unit
    auto prefix = getUnitPrefix( str.substr( 0, 1 ) );
    Units u     = readUnit( str.substr( 1 ), false );
    if ( prefix == UnitPrefix::unknown || u.d_scale == 0 ) {
        return readUnit( str );
    } else {
        u.d_scale *= convert( prefix );
        return u;
    }
    return Units( { 0 }, 0.0 );
}
constexpr Units Units::readUnit( const std::string_view &str, bool throwErr )
{
    auto create = []( UnitType type, double s = 1.0 ) {
        auto u = getSI( type );
        return Units( u, s );
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
    if ( str == "ohm" || str == "Ω" )
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
    if ( str == "becquerel" || str == "Bq" )
        return create( UnitType::frequency );
    if ( str == "gray" || str == "Gy" )
        return Units( { -2, 2, 0, 0, 0, 0, 0, 0, 0 }, 1 );
    if ( str == "sievert" || str == "Sv" )
        return Units( { -2, 2, 0, 0, 0, 0, 0, 0, 0 }, 1 );
    if ( str == "katal" || str == "kat" )
        return Units( { -1, 0, 0, 0, 0, 1, 0, 0, 0 }, 1 );
    if ( str == "litre" || str == "L" )
        return Units( { 0, 3, 0, 0, 0, 0, 0, 0, 0 }, 1e-3 );
    if ( str == "milliliters" || str == "ml" || str == "mL" )
        return Units( { 0, 3, 0, 0, 0, 0, 0, 0, 0 }, 1e-6 );
    // Non-SI units accepted for use with SI
    if ( str == "minute" || str == "minutes" )
        return create( UnitType::time, 60 );
    if ( str == "hour" || str == "hr" )
        return create( UnitType::time, 3600 );
    if ( str == "day" )
        return create( UnitType::time, 86400 );
    if ( str == "week" )
        return create( UnitType::time, 604800 );
    if ( str == "mmHg" )
        return create( UnitType::pressure, 133.322387415 );
    // Check cgs units
    if ( str == "ergs" || str == "erg" )
        return create( UnitType::energy, 1e-7 );
    if ( str == "eV" )
        return create( UnitType::energy, 1.602176634e-19 );
    if ( str == "dyn" )
        return create( UnitType::force, 1e-5 );
    if ( str == "barye" || str == "Ba" )
        return create( UnitType::pressure, 0.1 );
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
        return Units( { 0, 2, 0, 0, 0, 0, 0, 0, 0 }, 4046.8564224 );
    if ( str == "teaspoon" || str == "tsp" )
        return Units( { 0, 3, 0, 0, 0, 0, 0, 0, 0 }, 4.92892e-6 );
    if ( str == "tablespoon" || str == "tbsp" )
        return Units( { 0, 3, 0, 0, 0, 0, 0, 0, 0 }, 1.47868e-5 );
    if ( str == "cup" || str == "cp" )
        return Units( { 0, 3, 0, 0, 0, 0, 0, 0, 0 }, 2.36588e-4 );
    if ( str == "pint" || str == "pt" )
        return Units( { 0, 3, 0, 0, 0, 0, 0, 0, 0 }, 4.7317647274592e-4 );
    if ( str == "quart" || str == "qt" )
        return Units( { 0, 3, 0, 0, 0, 0, 0, 0, 0 }, 9.4635294549184e-4 );
    if ( str == "gallon" || str == "gal" )
        return Units( { 0, 3, 0, 0, 0, 0, 0, 0, 0 }, 3.78541178196736e-3 );
    if ( str == "ounce" || str == "oz" )
        return create( UnitType::mass, 0.028349523125 );
    if ( str == "pound" || str == "lb" )
        return create( UnitType::mass, 0.45359237 );
    if ( str == "ton" )
        return create( UnitType::mass, 1016.0469088 );
    if ( str == "lbf" )
        return create( UnitType::force, 4.4482216152605 );
    if ( str == "Rankine" || str == "R" )
        return create( UnitType::temperature, 5.0 / 9.0 );
    // Check atomic units
    if ( str == "hartree" )
        return create( UnitType::energy, 4.359744722207185e-18 );
    if ( str == "bohr" )
        return create( UnitType::length, 5.2917721090380e-11 );
    // Check special units/characters
    if ( str == "percent" || str == "%" )
        return create( UnitType::unitless, 0.01 );
    // No success
    SI_type u = { 0 };
    double s  = 0;
    if ( throwErr && s == 0 ) {
        // We have an error interpreting the string, try to interpret it as a double
        // Note: this is only valid when we are not using constexpr (also true for throw)
        s = strtod( str, false );
        if ( s == s && s != 0 )
            return Units( u, s );
        // Still unable, throw an error
        char err[256] = "Unknown unit: ";
        for ( size_t i = 0; i < str.size(); i++ )
            err[14 + i] = str[i];
        throw std::logic_error( err );
    }
    return Units( u, s );
}
inline std::vector<std::string> Units::getAllUnits()
{
    return { "second",   "s",           "meter",      "m",         "gram",    "g",
             "ampere",   "A",           "kelvin",     "K",         "mole",    "mol",
             "candela",  "cd",          "radian",     "radians",   "rad",     "steradian",
             "sr",       "degree",      "degrees",    "joule",     "J",       "watt",
             "W",        "hertz",       "Hz",         "newton",    "N",       "pascal",
             "Pa",       "coulomb",     "C",          "volt",      "V",       "farad",
             "F",        "ohm",         "Ω",          "siemens",   "S",       "weber",
             "Wb",       "tesla",       "T",          "henry",     "H",       "lumen",
             "lm",       "lux",         "lx",         "becquerel", "Bq",      "gray",
             "Gy",       "sievert",     "Sv",         "katal",     "kat",     "litre",
             "L",        "milliliters", "ml",         "mL",        "minute",  "minutes",
             "hour",     "hr",          "day",        "week",      "mmHg",    "ergs",
             "erg",      "eV",          "dyn",        "barye",     "Ba",      "inch",
             "in",       "\"",          "foot",       "ft",        "\'",      "yard",
             "yd",       "furlong",     "fur",        "mile",      "mi",      "acre",
             "teaspoon", "tsp",         "tablespoon", "tbsp",      "cup",     "cp",
             "pint",     "pt",          "quart",      "qt",        "gallon",  "gal",
             "ounce",    "oz",          "pound",      "lb",        "ton",     "lbf",
             "Rankine",  "R",           "hartree",    "bohr",      "percent", "%" };
}


/********************************************************************
 * Get unit type                                                     *
 ********************************************************************/
constexpr Units::SI_type Units::getSI( UnitType type )
{
    if ( type == UnitType::unknown )
        return { 0, 0, 0, 0, 0, 0, 0, 0, 0 };
    // s, m, kg, A, K, mol, cd, rad, sr
    constexpr SI_type id[25] = {
        { 0, 0, 0, 0, 0, 0, 0, 0, 0 },   // 0: unitless
        { 1, 0, 0, 0, 0, 0, 0, 0, 0 },   // 1: time (s)
        { 0, 1, 0, 0, 0, 0, 0, 0, 0 },   // 2: length (m)
        { 0, 0, 1, 0, 0, 0, 0, 0, 0 },   // 3: mass (kg)
        { 0, 0, 0, 1, 0, 0, 0, 0, 0 },   // 4: current (A)
        { 0, 0, 0, 0, 1, 0, 0, 0, 0 },   // 5: temperature (K)
        { 0, 0, 0, 0, 0, 1, 0, 0, 0 },   // 6: mole (mol)
        { 0, 0, 0, 0, 0, 0, 1, 0, 0 },   // 7: intensity (cd)
        { 0, 0, 0, 0, 0, 0, 0, 1, 0 },   // 8: angle (rad)
        { 0, 0, 0, 0, 0, 0, 0, 0, 1 },   // 9: solid angle (sr)
        { -2, 2, 1, 0, 0, 0, 0, 0 },     // 10: energy (kg m2 s-2)
        { -3, 2, 1, 0, 0, 0, 0, 0 },     // 11: power (kg m2 s-3)
        { -1, 0, 0, 0, 0, 0, 0, 0 },     // 12: frequency (s-1)
        { -2, 1, 1, 0, 0, 0, 0, 0, 0 },  // 13: force (kg m s-2)
        { -2, -1, 1, 0, 0, 0, 0, 0, 0 }, // 14: pressure (kg m-1 s-2)
        { 1, 0, 0, 1, 0, 0, 0, 0, 0 },   // 15: electric charge (s A)
        { -3, 2, 1, -1, 0, 0, 0, 0, 0 }, // 16: electrical potential (kg m2 s-3 A-1)
        { 4, -2, -1, 2, 0, 0, 0, 0, 0 }, // 17: capacitance (kg-1 m-2 s4 A2)
        { -3, 2, 1, -2, 0, 0, 0, 0, 0 }, // 18: resistance (kg m2 s-3 A-2)
        { 3, -2, -1, 2, 0, 0, 0, 0, 0 }, // 19: electrical conductance (kg-1 m-2 s3 A2)
        { -2, 2, 1, -1, 0, 0, 0, 0, 0 }, // 20: magnetic flux (kg m2 s-2 A-1)
        { -2, 0, 1, -1, 0, 0, 0, 0, 0 }, // 21: magnetic flux density (kg s-2 A-1)
        { -2, 2, 1, -2, 0, 0, 0, 0, 0 }, // 22: inductance (kg m2 s-2 A-2)
        { 0, 0, 0, 0, 0, 0, 1, 0, 1 },   // 23: luminous flux (cd sr)
        { 0, -2, 0, 0, 0, 0, 1, 0, 1 }   // 24: illuminance (cd sr m-2)
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
    constexpr SI_type energy      = { -2, 2, 1, 0, 0, 0, 0, 0 };
    constexpr SI_type temperature = { 0, 0, 0, 0, 1, 0, 0, 0, 0 };
    if ( d_SI == rhs.d_SI )
        return true;
    if ( d_SI == energy && rhs.d_SI == temperature )
        return true; // Convert from energy to temperature (using boltzmann constant)
    if ( d_SI == temperature && rhs.d_SI == energy )
        return true; // Convert from temperature to energy (using boltzmann constant)
    return false;
}
constexpr double Units::convert( const Units &rhs ) const
{
    constexpr SI_type energy      = { -2, 2, 1, 0, 0, 0, 0, 0 };
    constexpr SI_type temperature = { 0, 0, 0, 0, 1, 0, 0, 0, 0 };
    if ( d_SI == rhs.d_SI ) {
        // The SI units match
        return d_scale / rhs.d_scale;
    } else if ( d_SI == energy && rhs.d_SI == temperature ) {
        // Convert from energy to temperature (using bolztmann constant)
        return 7.242971666663E+22 * d_scale / rhs.d_scale;
    } else if ( d_SI == temperature && rhs.d_SI == energy ) {
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
    double err = d_scale >= rhs.d_scale ? d_scale - rhs.d_scale : rhs.d_scale - d_scale;
    bool test  = err < 1e-10 * d_scale;
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
    d_scale /= rhs.d_scale;
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
constexpr Units Units::pow( int exponent ) const noexcept
{
    if ( exponent == 0 )
        return Units( { 0 }, 1.0 );
    auto base = *this;
    if ( exponent < 0 ) {
        base     = Units( { 0 }, 1.0 ) / base;
        exponent = -exponent;
    }
    Units u = base;
    while ( exponent > 1 ) {
        u *= base;
        exponent -= 1;
    }
    return u;
}


} // namespace AMP

#endif
