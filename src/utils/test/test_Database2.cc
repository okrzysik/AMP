#include "AMP/utils/AMPManager.h"
#include "AMP/utils/AMP_MPI.h"
#include "AMP/utils/Database.h"
#include "AMP/utils/PIO.h"
#include "AMP/utils/UnitTest.h"
#include "AMP/utils/Utilities.h"
#include "AMP/utils/shared_ptr.h"

#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <sys/stat.h>
#include <vector>


/************************************************************************
 *                                                                       *
 * This tests whether we can create and use an InputManager object       *
 *                                                                       *
 ************************************************************************/
void readInputDatabase( AMP::UnitTest &ut )
{
    std::string input_file = "input_Database";
    std::string log_file   = "output_Database";

    // Create input database and parse all data in input file.
    auto input_db = AMP::Database::parseInputFile( input_file );

    auto tmp_db = input_db->getDatabase( "Try" );
    int number  = tmp_db->getScalar<int>( "number" );

    if ( number > 0 ) {
        auto intArray = tmp_db->getVector<int>( "intArray" );
        if ( (int) intArray.size() != number )
            ut.failure( "intArray was the wrong size" );
        auto doubleArray = tmp_db->getVector<double>( "doubleArray" );
        if ( (int) doubleArray.size() != number )
            ut.failure( "doubleArray was the wrong size" );
    }

    std::cout << "tmp_db created and destroyed successfully." << std::endl;

    input_db.reset();

    ut.passes( "Get Database Works and the Destructor of AMP::Database works." );
}


/************************************************************************
 *                                                                       *
 * This tests whether we can put/get keys with a database                *
 *                                                                       *
 ************************************************************************/
void testCreateDatabase( AMP::UnitTest &ut )
{
    auto db = AMP::make_shared<AMP::Database>( "database" );

    std::complex<double> zero( 0, 0 );
    std::complex<double> onetwo( 1, 2 );

    db->putScalar( "scalar_int", (int) 1 );
    db->putScalar( "scalar_float", (float) 1 );
    db->putScalar( "scalar_double", (double) 1 );
    db->putScalar( "scalar_complex", onetwo );
    db->putScalar( "scalar_char", (char) 1 );
    db->putScalar( "scalar_bool", true );

    AMP_ASSERT( db->keyExists( "scalar_int" ) );

    AMP_ASSERT( db->isType<int>( "scalar_int" ) );
    AMP_ASSERT( db->isType<float>( "scalar_float" ) );
    AMP_ASSERT( db->isType<double>( "scalar_double" ) );
    AMP_ASSERT( db->isType<std::complex<double>>( "scalar_complex" ) );
    AMP_ASSERT( db->isType<char>( "scalar_char" ) );
    AMP_ASSERT( db->isType<bool>( "scalar_bool" ) );

    AMP_ASSERT( db->getScalar<int>( "scalar_int" ) == 1 );
    AMP_ASSERT( db->getScalar<float>( "scalar_float" ) == 1.0 );
    AMP_ASSERT( db->getScalar<double>( "scalar_double" ) == 1.0 );
    AMP_ASSERT( db->getScalar<std::complex<double>>( "scalar_complex" ) == onetwo );
    AMP_ASSERT( db->getScalar<char>( "scalar_char" ) == 1 );
    AMP_ASSERT( db->getScalar<bool>( "scalar_bool" ) == true );

    AMP_ASSERT( db->getWithDefault<int>( "scalar_int", 0 ) == 1 );
    AMP_ASSERT( db->getWithDefault<float>( "scalar_float", 0.0 ) == 1.0 );
    AMP_ASSERT( db->getWithDefault<double>( "scalar_double", 0 ) == 1.0 );
    AMP_ASSERT( db->getWithDefault<std::complex<double>>( "scalar_complex", zero ) == onetwo );
    AMP_ASSERT( db->getWithDefault<char>( "scalar_char", 0 ) == 1 );
    AMP_ASSERT( db->getWithDefault<bool>( "scalar_bool", false ) == true );

    AMP_ASSERT( db->getVector<int>( "scalar_int" ).size() == 1 );
    AMP_ASSERT( db->getVector<float>( "scalar_float" ).size() == 1 );
    AMP_ASSERT( db->getVector<double>( "scalar_double" ).size() == 1 );
    AMP_ASSERT( db->getVector<std::complex<double>>( "scalar_complex" ).size() == 1 );
    AMP_ASSERT( db->getVector<char>( "scalar_char" ).size() == 1 );
    AMP_ASSERT( db->getVector<bool>( "scalar_bool" ).size() == 1 );

    // Test base class functions
    db->AMP::Database::putScalar( "int2", (int) 1 );
    db->AMP::Database::putScalar( "float2", (float) 1 );
    db->AMP::Database::putScalar( "double2", (double) 1 );
    db->AMP::Database::putScalar( "complex2", onetwo );
    db->AMP::Database::putScalar( "char2", (char) 1 );
    db->AMP::Database::putScalar( "bool2", true );
    AMP_ASSERT( db->AMP::Database::getScalar<int>( "int2" ) == 1 );
    AMP_ASSERT( db->AMP::Database::getScalar<float>( "float2" ) == 1.0 );
    AMP_ASSERT( db->AMP::Database::getScalar<double>( "double2" ) == 1.0 );
    AMP_ASSERT( db->AMP::Database::getScalar<std::complex<double>>( "complex2" ) == onetwo );
    AMP_ASSERT( db->AMP::Database::getScalar<char>( "char2" ) == 1 );
    AMP_ASSERT( db->AMP::Database::getScalar<bool>( "bool2" ) == true );
    AMP_ASSERT( db->AMP::Database::getWithDefault<int>( "scalar_int", 0 ) == 1 );
    AMP_ASSERT( db->AMP::Database::getWithDefault<float>( "scalar_float", 0 ) == 1.0 );
    AMP_ASSERT( db->AMP::Database::getWithDefault<double>( "scalar_double", 0 ) == 1.0 );
    AMP_ASSERT( db->AMP::Database::getWithDefault<std::complex<double>>( "scalar_complex", zero ) ==
                onetwo );
    AMP_ASSERT( db->AMP::Database::getWithDefault<char>( "scalar_char", 0 ) == 1 );
    AMP_ASSERT( db->AMP::Database::getWithDefault<bool>( "scalar_bool", false ) == true );

    // Finished
    ut.passes( "Create database passes" );
}


//---------------------------------------------------------------------------//

int main( int argc, char *argv[] )
{
    AMP::AMP_MPI::start_MPI( argc, argv );
    AMP::UnitTest ut;

    readInputDatabase( ut );
    testCreateDatabase( ut );

    ut.report();

    int num_failed = ut.NumFailGlobal();
    AMP::AMP_MPI::stop_MPI();
    return num_failed;
}
