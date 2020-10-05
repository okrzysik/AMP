#include "AMP/utils/AMPManager.h"
#include "AMP/utils/Database.h"
#include "AMP/utils/UnitTest.h"
#include "AMP/utils/Utilities.h"

#include <chrono>
#include <complex>
#include <fstream>
#include <random>
#include <sstream>


using namespace AMP;


std::default_random_engine generator( 123 );


bool equal( double a, double b ) { return fabs( a - b ) <= 1e-12 * fabs( a + b ); }


// Generate random number
template<class TYPE>
static typename std::enable_if<std::is_floating_point<TYPE>::value, TYPE>::type random()
{
    static std::uniform_real_distribution<double> dist( 0, 1 );
    return static_cast<TYPE>( dist( generator ) );
}
template<class TYPE>
static typename std::enable_if<std::is_integral<TYPE>::value, TYPE>::type random()
{
    static std::uniform_int_distribution<int> dist( 0, 1 );
    return static_cast<TYPE>( dist( generator ) );
}
template<class TYPE>
static typename std::enable_if<std::is_same<TYPE, std::complex<double>>::value,
                               std::complex<double>>::type
random()
{
    return std::complex<double>( random<double>(), random<double>() );
}
template<class TYPE>
static typename std::enable_if<std::is_same<TYPE, std::complex<float>>::value,
                               std::complex<float>>::type
random()
{
    return std::complex<float>( random<float>(), random<float>() );
}


// Test a type in the database
template<class TYPE>
static void addType( Database &db, UnitTest &ut )
{
    bool pass            = true;
    std::string typeName = typeid( TYPE ).name();
    TYPE rand            = random<TYPE>();
    db.putScalar<TYPE>( "scalar-" + typeName, rand );
    db.putVector<TYPE>( "vector-" + typeName, { rand } );
    auto v1 = db.getScalar<TYPE>( "scalar-" + typeName );
    auto v2 = db.getVector<TYPE>( "vector-" + typeName );
    pass    = pass && v1 == rand && v2.size() == 1 && v2[0] == rand;
    pass    = pass && db.isType<TYPE>( "scalar-" + typeName );
    pass    = pass && db.isType<TYPE>( "vector-" + typeName );
    if ( std::is_floating_point<TYPE>::value ) {
        auto v3 = db.getScalar<double>( "scalar-" + typeName );
        auto v4 = db.getVector<double>( "vector-" + typeName );
        pass    = pass && static_cast<TYPE>( v3 ) == rand;
        pass    = pass && v4.size() == 1 && static_cast<TYPE>( v4[0] ) == rand;
        pass    = pass && db.isType<double>( "scalar-" + typeName );
        pass    = pass && !db.isType<int>( "scalar-" + typeName );
    } else if ( std::is_integral<TYPE>::value ) {
        auto v3 = db.getScalar<int>( "scalar-" + typeName );
        auto v4 = db.getVector<int>( "vector-" + typeName );
        pass    = pass && static_cast<TYPE>( v3 ) == rand;
        pass    = pass && v4.size() == 1 && static_cast<TYPE>( v4[0] ) == rand;
        pass    = pass && db.isType<double>( "scalar-" + typeName );
        pass    = pass && db.isType<int>( "scalar-" + typeName );
    }
    if ( !pass )
        ut.failure( typeName );
}


// Run some basic tests
void runBasicTests( UnitTest &ut )
{
    // Create a database with some different types of data
    Database db;
    db.putScalar<std::string>( "string", "test" );
    db.putScalar<bool>( "true", true );
    db.putScalar<bool>( "false", false );
    db.putScalar<double>( "double", 3.14 );
    db.putScalar<int>( "int", -2 );
    db.putScalar<double>( "inf", std::numeric_limits<double>::infinity() );
    db.putScalar<double>( "nan", std::numeric_limits<double>::quiet_NaN() );
    db.putVector<double>( "x", { 1.0, 2.0, 3.0 } );
    db.putScalar<double>( "x1", 1.5, "cm" );
    db.putVector<double>( "x2", { 2.5 }, "mm" );

    // Test adding some different types
    addType<uint8_t>( db, ut );
    addType<uint16_t>( db, ut );
    addType<uint32_t>( db, ut );
    addType<uint64_t>( db, ut );
    addType<int8_t>( db, ut );
    addType<int16_t>( db, ut );
    addType<int32_t>( db, ut );
    addType<int64_t>( db, ut );
    addType<float>( db, ut );
    addType<double>( db, ut );
    addType<long double>( db, ut );
    addType<std::complex<double>>( db, ut );

    // Try to read/add a database that ends in a comment and has units
    const char beamDatabaseText[] =
        "prof = \"gaussian\"     // Temporal profile\n"
        "z_shape = \"hard_hat\"  // Shape along the line\n"
        "y_shape = \"gaussian\"  // Shape perpendicular to the line\n"
        "geometry = 1            // Beam geometry\n"
        "x      = 0 cm           // Position\n"
        "FWHM   = 120 ps         // Full Width Half Max\n"
        "E      = 0.4 J          // Energy in the beam\n"
        "delay  = 0 ps           // Delay the beam respect to zero\n"
        "diam   = 30 um          // Beam diameter\n"
        "length = 0.4 cm         // Beam length\n"
        "lambda = 0.8 um         // Wavelength of laser\n"
        "angle  = 0 degrees      // Angle of beam with repect to normal\n";
    auto beam = Database::createFromString( beamDatabaseText );
    db.putDatabase( "beam", std::move( beam ) );

    // Write the database to a file
    std::ofstream inputfile;
    inputfile.open( "test_Database.out" );
    db.print( inputfile );
    inputfile.close();

    // Read the database and check that everything matches
    auto db2 = Database::parseInputFile( "test_Database.out" );
    if ( db2->isDatabase( "string" ) )
        ut.failure( "string is a database?" );
    if ( !db2->isType<std::string>( "string" ) ||
         db2->getScalar<std::string>( "string" ) != "test" )
        ut.failure( "string" );
    if ( !db2->isType<bool>( "true" ) || !db2->getScalar<bool>( "true" ) )
        ut.failure( "true" );
    if ( !db2->isType<bool>( "false" ) || db2->getScalar<bool>( "false" ) )
        ut.failure( "false" );
    if ( !db2->isType<double>( "double" ) || db2->getScalar<double>( "double" ) != 3.14 )
        ut.failure( "double" );
    if ( !db2->isType<int>( "int" ) || db2->getScalar<int>( "int" ) != -2 )
        ut.failure( "int" );
    if ( db2->getScalar<double>( "inf" ) != std::numeric_limits<double>::infinity() )
        ut.failure( "inf" );
    if ( db2->getScalar<double>( "nan" ) == db2->getScalar<double>( "nan" ) )
        ut.failure( "nan" );
    if ( !equal( db2->getScalar<double>( "x1", "m" ), 0.015 ) ||
         !equal( db2->getScalar<double>( "x2", "m" ), 0.0025 ) ||
         !equal( db2->getVector<double>( "x1", "um" )[0], 15000 ) ||
         !equal( db2->getVector<double>( "x2", "um" )[0], 2500 ) )
        ut.failure( "units" );
    if ( db == *db2 )
        ut.passes( "print" );
    else
        ut.failure( "print" );

    // Add more types that are not compatible with print
    addType<std::complex<float>>( db, ut );

    // Try to clone the database
    auto db3 = db.cloneDatabase();
    if ( db == *db3 )
        ut.passes( "clone" );
    else
        ut.failure( "clone" );

    // Test copy/move operators
    bool pass = true;
    Database db4( std::move( *db3 ) );
    pass = pass && db == db4;
    Database db5;
    db5  = std::move( db4 );
    pass = pass && db == db5;
    Database db6;
    db6.copy( db5 );
    pass = pass && db == db5;
    if ( pass )
        ut.passes( "copy/move" );
    else
        ut.failure( "copy/move" );
}


// Run tests using an input file
void runFileTests( UnitTest &ut, const std::string &filename )
{
    // Read the file
    auto t1 = std::chrono::system_clock::now();
    auto db = Database::parseInputFile( filename );
    auto t2 = std::chrono::system_clock::now();
    int us  = std::chrono::duration_cast<std::chrono::microseconds>( t2 - t1 ).count();
    printf( "Time to read %s: %i us\n", filename.c_str(), us );
    // Create database from string
    std::ifstream ifstream( filename );
    std::stringstream buffer;
    buffer << ifstream.rdbuf();
    t1       = std::chrono::system_clock::now();
    auto db2 = Database::createFromString( buffer.str() );
    t2       = std::chrono::system_clock::now();
    us       = std::chrono::duration_cast<std::chrono::microseconds>( t2 - t1 ).count();
    printf( "Time to create from string: %i us\n", us );
    // Clone the database
    t1       = std::chrono::system_clock::now();
    auto db3 = db->cloneDatabase();
    t2       = std::chrono::system_clock::now();
    us       = std::chrono::duration_cast<std::chrono::microseconds>( t2 - t1 ).count();
    printf( "Time to clone database: %i us\n", us );
    // Check the database
    bool pass = !db->empty();
    t1        = std::chrono::system_clock::now();
    pass      = pass && *db2 == *db;
    pass      = pass && *db3 == *db;
    t2        = std::chrono::system_clock::now();
    us        = std::chrono::duration_cast<std::chrono::microseconds>( t2 - t1 ).count();
    printf( "Time to compare database: %i us\n", us / 2 );
    // Check that data loaded
    if ( filename == "laser_plasma_input.txt" ) {
        pass = !db->getVector<std::string>( "material" ).empty();
    }
    // Finished test
    if ( pass )
        ut.passes( filename );
    else
        ut.failure( filename );
}


// Main
int main( int argc, char *argv[] )
{
    AMPManager::startup( argc, argv );
    UnitTest ut;

    // Run the tests
    runBasicTests( ut );
    for ( int i = 1; i < argc; i++ )
        runFileTests( ut, argv[i] );

    // Return
    int N_errors = ut.NumFailGlobal();
    ut.report();
    ut.reset();
    AMPManager::shutdown();
    return N_errors;
}
