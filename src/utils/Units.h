#ifndef included_AMP_Units
#define included_AMP_Units

#include <array>
#include <math.h>
#include <string_view>
#include <utility>


namespace AMP {

//! Enum to hold prefix
enum class UnitPrefix : int8_t {
    yocto   = 0,
    zepto   = 1,
    atto    = 2,
    femto   = 3,
    pico    = 4,
    nano    = 5,
    micro   = 6,
    milli   = 7,
    centi   = 8,
    deci    = 9,
    none    = 10,
    deca    = 11,
    hecto   = 12,
    kilo    = 13,
    mega    = 14,
    giga    = 15,
    tera    = 16,
    peta    = 17,
    exa     = 18,
    zetta   = 19,
    yotta   = 20,
    unknown = 21
};


//! Enum to hold type
enum class UnitType : int8_t {
    unitless              = 0,
    time                  = 1,
    length                = 2,
    mass                  = 3,
    current               = 4,
    temperature           = 5,
    mole                  = 6,
    intensity             = 7,
    angle                 = 8,
    solidAngle            = 9,
    energy                = 10,
    power                 = 11,
    frequency             = 12,
    force                 = 13,
    pressure              = 14,
    electricCharge        = 15,
    electricalPotential   = 16,
    capacitance           = 17,
    resistance            = 18,
    electricalConductance = 19,
    magneticFlux          = 20,
    magneticFluxDensity   = 21,
    inductance            = 22,
    luminousFlux          = 23,
    illuminance           = 24,
    unknown               = -1
};


/**
 * \class Units
 * \brief  Provides a class for storing units
 * \details  This class provides routines for creating and storing units.
 *    A user can specify a unit and get the scaling factor to a compatible unit.
 *    All user-provided units are converted to SI base units with a scaling factor.
 *    Currently we support all major SI units and some additional units with prefixes.
 *    Note: currently the majority of the routines are constexpr so that the creation
 *       of a unit can be done at compile-time.
 */
class alignas( 8 ) Units final
{
public:
    //! Empty constructor
    constexpr Units();

    /**
     * \brief  Construct the units from a const char array
     * \details  This is the default constructor for a const char array.
     *    It can create a unit from a string of the format "W/(m^2)"
     * \param unit      Input string
     */
    constexpr Units( const char *unit );

    /**
     * \brief  Construct the units from a const char array
     * \details  This is the default constructor for a string view.
     *    It can create a unit from a string of the format "W/(m^2)"
     * \param unit      Input string
     */
    constexpr Units( const std::string_view &unit );

    /**
     * \brief  Get the unit type
     * \details  This returns the type of unit (if it is reducible).
     *   e.g. time, length, energy, power, etc.
     */
    constexpr UnitType getType() const noexcept;

    /**
     * \brief  Get the prefix from a string
     * \details  This returns an enum representing the given input string
     * \param str       Input string
     */
    constexpr static UnitPrefix getUnitPrefix( const std::string_view &str ) noexcept;

    /**
     * \brief  Check if two units are compatible
     * \details  Check if two units are compatible with each other.
     *     I.e. if convert will succeed.
     * \param unit      Unit to compare
     */
    constexpr bool compatible( const Units &unit ) noexcept;

    /**
     * \brief  Convert the unit to a new type
     * \details  This function returns the scaling factor to convert between
     *    two different units.  For example, converting "J" to "ergs" returns 1e-7.
     *    If the two units are not compatible an exception will be thrown.
     * \param unit      Desired unit
     */
    constexpr double convert( const Units &unit ) const;

    /**
     * \brief  Convert the prefix to a double
     * \details  This returns the value of the given prefix
     * \param x         Enum for the prefix
     */
    constexpr static double convert( UnitPrefix x ) noexcept
    {
        return d_pow10[static_cast<int8_t>( x )];
    }

    //! Operator ==
    constexpr bool operator==( const Units &rhs ) const noexcept;

    //! Operator !=
    constexpr bool operator!=( const Units &rhs ) const noexcept { return !operator==( rhs ); }

    //! Operator *=
    constexpr void operator*=( const Units &rhs ) noexcept;

    //! Operator /=
    constexpr void operator/=( const Units &rhs ) noexcept;

    //! Raise the unit to the given power
    constexpr Units pow( int p ) const noexcept;

    //! Check if unit is null
    constexpr bool isNull() const { return d_scale == 0; }


public:
    //! Get a string representation of the units
    std::string str() const;

    //! Get a string representation for the prefix
    static std::string_view getPrefixStr( UnitPrefix ) noexcept;

    //! Get a string representation of the units in SI units (with scaling factor)
    std::string printSI() const;

    //! Get a string representation of the units in SI units (with scaling factor)
    std::string printUnit() const;

    //! Get the full unit and conversion string
    std::string printFull() const;


public:
    static constexpr int atoi( std::string_view str, bool throw_error = true );
    static constexpr double strtod( std::string_view str, bool throw_error = true );

protected:
    using unit_type = std::array<char, 31>;
    using SI_type   = std::array<int8_t, 9>;

    unit_type d_unit;
    SI_type d_SI;
    double d_scale;

protected:
    constexpr Units( const SI_type &u, double s );

    static constexpr Units read( std::string_view str );
    static constexpr Units read2( std::string_view str );
    static constexpr Units readUnit( const std::string_view &str, bool throwErr = true );
    static constexpr SI_type combine( const SI_type &a, const SI_type &b );
    static constexpr SI_type getSI( UnitType );

    std::string printSIBase() const;

    static constexpr std::pair<size_t, char> findToken( const std::string_view &, size_t );
    static constexpr size_t findPar( const std::string_view &, size_t );

private:
    static constexpr double d_pow10[22] = { 1e-24, 1e-21, 1e-18, 1e-15, 1e-12, 1e-9, 1e-6, 1e-3,
                                            1e-2,  0.1,   1,     10,    100,   1000, 1e6,  1e9,
                                            1e12,  1e15,  1e18,  1e21,  1e24,  0 };
    static constexpr const char *d_SI_units[] = {
        "s", "m", "kg", "A", "K", "mol", "cd", "rad", "sr"
    };

protected:
    friend constexpr Units pow( Units base, int exponent );
};


} // namespace AMP


#include "AMP/utils/Units.hpp"


#endif
