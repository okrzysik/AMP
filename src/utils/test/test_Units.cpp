#include "AMP/utils/AMPManager.h"
#include "AMP/utils/UnitTest.h"
#include "AMP/utils/Units.h"

using namespace AMP;

using prefix = Units::UnitPrefix;
using unit   = Units::UnitValue;
using type   = Units::UnitType;

inline void testUnit( const std::string &s, const Units &u2, double prefix, UnitTest &ut )
{
    Units u1( s );
    bool pass = u1 == u2;
    double x  = u1.convert( Units( Units::UnitPrefix::none, u1.getUnit() ) );
    double y  = Units::convert( u1.getPrefix() );
    auto p2   = Units::getUnitPrefix( Units::str( u1.getPrefix() ).data() );
    pass      = pass && x == y;
    pass      = pass && x == prefix;
    pass      = pass && p2 == u1.getPrefix();
    if ( pass )
        ut.passes( s );
    else
        ut.failure( s );
}
inline void testUnitValue( UnitTest &ut )
{
    bool pass     = true;
    auto getValue = Units::getUnitValue;
    auto getType  = []( const char *u ) { return Units::getUnitType( Units::getUnitValue( u ) ); };
    pass          = pass && getValue( "meter" ) == unit::meter;
    pass          = pass && getValue( "gram" ) == unit::gram;
    pass          = pass && getValue( "second" ) == unit::second;
    pass          = pass && getValue( "ampere" ) == unit::ampere;
    pass          = pass && getValue( "kelvin" ) == unit::kelvin;
    pass          = pass && getValue( "joule" ) == unit::joule;
    pass          = pass && getValue( "erg" ) == unit::erg;
    pass          = pass && getValue( "degree" ) == unit::degree;
    pass          = pass && getValue( "radian" ) == unit::radian;
    pass          = pass && getType( "meter" ) == type::length;
    pass          = pass && getType( "gram" ) == type::mass;
    pass          = pass && getType( "second" ) == type::time;
    pass          = pass && getType( "ampere" ) == type::current;
    pass          = pass && getType( "kelvin" ) == type::temperature;
    pass          = pass && getType( "joule" ) == type::energy;
    pass          = pass && getType( "erg" ) == type::energy;
    pass          = pass && getType( "degree" ) == type::angle;
    pass          = pass && getType( "radian" ) == type::angle;
    if ( pass )
        ut.passes( "Unit Types" );
    else
        ut.failure( "Unit Types" );
}


// Main
int main( int argc, char *argv[] )
{
    AMPManager::startup( argc, argv );
    UnitTest ut;

    // Test some different unit conversions
    testUnit( "ys", Units( prefix::yocto, unit::second ), 1e-24, ut );
    testUnit( "zs", Units( prefix::zepto, unit::second ), 1e-21, ut );
    testUnit( "as", Units( prefix::atto, unit::second ), 1e-18, ut );
    testUnit( "fs", Units( prefix::femto, unit::second ), 1e-15, ut );
    testUnit( "ps", Units( prefix::pico, unit::second ), 1e-12, ut );
    testUnit( "ns", Units( prefix::nano, unit::second ), 1e-9, ut );
    testUnit( "us", Units( prefix::micro, unit::second ), 1e-6, ut );
    testUnit( "ms", Units( prefix::milli, unit::second ), 1e-3, ut );
    testUnit( "cs", Units( prefix::centi, unit::second ), 1e-2, ut );
    testUnit( "ds", Units( prefix::deci, unit::second ), 0.1, ut );
    testUnit( "s", Units( prefix::none, unit::second ), 1, ut );
    testUnit( "das", Units( prefix::deca, unit::second ), 10, ut );
    testUnit( "hs", Units( prefix::hecto, unit::second ), 100, ut );
    testUnit( "ks", Units( prefix::kilo, unit::second ), 1e3, ut );
    testUnit( "Ms", Units( prefix::mega, unit::second ), 1e6, ut );
    testUnit( "Gs", Units( prefix::giga, unit::second ), 1e9, ut );
    testUnit( "Ts", Units( prefix::tera, unit::second ), 1e12, ut );
    testUnit( "Ps", Units( prefix::peta, unit::second ), 1e15, ut );
    testUnit( "Es", Units( prefix::exa, unit::second ), 1e18, ut );
    testUnit( "Zs", Units( prefix::zetta, unit::second ), 1e21, ut );
    testUnit( "Ys", Units( prefix::yotta, unit::second ), 1e24, ut );

    // Test the unit types
    testUnitValue( ut );

    // Return
    int N_errors = ut.NumFailGlobal();
    ut.report();
    ut.reset();
    AMPManager::shutdown();
    return N_errors;
}
