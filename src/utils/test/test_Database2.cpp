#include "AMP/utils/AMPManager.h"
#include "AMP/utils/AMP_MPI.h"
#include "AMP/utils/Database.h"
#include "AMP/utils/PIO.h"
#include "AMP/utils/UnitTest.h"
#include "AMP/utils/Utilities.h"

#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <memory>
#include <sstream>
#include <string>
#include <sys/stat.h>
#include <vector>

#ifdef USE_SAMRAI
    #include "SAMRAI/tbox/InputManager.h"
    #include "SAMRAI/tbox/MemoryDatabase.h"
#endif


/************************************************************************
 * This tests whether we can read a basic input file                     *
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
 * This tests whether we can put/get keys with a database                *
 ************************************************************************/
void testCreateDatabase( AMP::UnitTest &ut )
{
    auto db = std::make_shared<AMP::Database>( "database" );

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


/************************************************************************
 * This tests converting to/from SAMRAI                                  *
 ************************************************************************/
#if USE_SAMRAI
bool compare_SAMRAI( SAMRAI::tbox::Database &db1, SAMRAI::tbox::Database &db2 )
{
    // Check that the keys match
    auto keys1 = db1.getAllKeys();
    auto keys2 = db2.getAllKeys();
    std::sort( keys1.begin(), keys1.end() );
    std::sort( keys2.begin(), keys2.end() );
    if ( keys1 != keys2 )
        return false;
    for ( const auto &key : keys1 ) {
        auto type1     = db1.getArrayType( key );
        auto type2     = db2.getArrayType( key );
        using DataType = SAMRAI::tbox::Database::DataType;
        if ( type1 != type2 )
            return false;
        if ( type1 == DataType::SAMRAI_DATABASE ) {
            compare_SAMRAI( *db1.getDatabase( key ), *db2.getDatabase( key ) );
        } else if ( type1 == DataType::SAMRAI_BOOL ) {
            return db1.getBoolVector( key ) == db2.getBoolVector( key );
        } else if ( type1 == DataType::SAMRAI_CHAR ) {
            return db1.getCharVector( key ) == db2.getCharVector( key );
        } else if ( type1 == DataType::SAMRAI_INT ) {
            return db1.getIntegerVector( key ) == db2.getIntegerVector( key );
        } else if ( type1 == DataType::SAMRAI_COMPLEX ) {
            return db1.getComplexVector( key ) == db2.getComplexVector( key );
        } else if ( type1 == DataType::SAMRAI_DOUBLE ) {
            return db1.getDoubleVector( key ) == db2.getDoubleVector( key );
        } else if ( type1 == DataType::SAMRAI_FLOAT ) {
            return db1.getFloatVector( key ) == db2.getFloatVector( key );
        } else if ( type1 == DataType::SAMRAI_STRING ) {
            return db1.getStringVector( key ) == db2.getStringVector( key );
        } else if ( type1 == DataType::SAMRAI_BOX ) {
            return db1.getDatabaseBoxVector( key ) == db2.getDatabaseBoxVector( key );
        } else {
            AMP_ERROR( "Unknown type" );
        }
    }
    return true;
}
void testSAMRAI( AMP::UnitTest &ut, const std::string &inputFile )
{
    bool pass = true;

    // Read the input database from SAMRAI
    auto SAMRAI_db = std::make_shared<SAMRAI::tbox::MemoryDatabase>( "input_SAMRAI" );
    SAMRAI::tbox::InputManager::getManager()->parseInputFile( inputFile, SAMRAI_db );

    // Read the input database through the default reader
    auto AMP_db = AMP::Database::parseInputFile( inputFile );

    // Convert AMP to SAMRAI to AMP and compare
    AMP::Database db1( AMP_db->cloneToSAMRAI() );
    if ( db1 != *AMP_db ) {
        pass = false;
        ut.failure( "Read AMP dbdatabase, convert to SAMRAI then back to AMP" );
    }

    // Convert SAMRAI to AMP to SAMRAI and compare
    auto db2 = std::make_shared<AMP::Database>( *SAMRAI_db )->cloneToSAMRAI();
    if ( !compare_SAMRAI( *db2, *SAMRAI_db ) ) {
        pass = false;
        ut.failure( "Read SAMRAI database, convert to AMP then back to SAMRAI" );
    }

    // Convert SAMRAI's database to AMP and compare to AMP's reader
    AMP::Database db3( SAMRAI_db );
    if ( db3 != *AMP_db ) {
        pass = false;
        ut.failure( "Read AMP and SAMRAI databases, compare in AMP" );
    }

    // Convert the AMP database to SAMRAI and compare to AMP's reader
    auto db4 = AMP_db->cloneToSAMRAI();
    if ( !compare_SAMRAI( *db4, *SAMRAI_db ) ) {
        pass = false;
        ut.failure( "Read AMP and SAMRAI databases, compare in SAMRAI" );
    }

    // Check readThroughSAMRAI
    auto db5 = AMP::Database::readThroughSAMRAI( inputFile )->cloneToSAMRAI();
    if ( !compare_SAMRAI( *db5, *SAMRAI_db ) ) {
        pass = false;
        ut.failure( "readThroughSAMRAI" );
    }

    if ( pass )
        ut.passes( "Converting between AMP and SAMRAI databases" );
}
#else
void testSAMRAI( AMP::UnitTest &, const std::string & ) {}
#endif


/************************************************************************
 * Main                                                                  *
 ************************************************************************/
int main( int argc, char *argv[] )
{
    AMP::AMPManager::startup( argc, argv );
    AMP::UnitTest ut;

    readInputDatabase( ut );
    testCreateDatabase( ut );
    testSAMRAI( ut, "input_SAMRAI" );
    testSAMRAI( ut, "input_NRDF" );
    testSAMRAI( ut, "input_NRDF2" );

    ut.report();

    int num_failed = ut.NumFailGlobal();
    AMP::AMPManager::shutdown();
    return num_failed;
}
