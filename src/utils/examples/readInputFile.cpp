#include "AMP/utils/AMPManager.h"
#include "AMP/utils/Database.h"
#include "AMP/utils/PIO.h"

int main( int argc, char **argv )
{
    // Every AMP program starts with a startup call
    // Internally AMP then initializes all the
    // third party libraries it is linking in
    AMP::AMPManager::startup( argc, argv );

    std::string input_file = "input_readInputFile";

    // AMP reads an input file into an AMP::Database object
    // Every database consists of a set of key value pairs
    // The key values can be scalar types, arrays of scalar types,
    // or other databases
    auto inputDB = AMP::Database::parseInputFile( input_file );

    if ( inputDB->keyExists( "a" ) ) {
        auto a = inputDB->getScalar<int>( "a" );
        AMP::pout << "a: " << a << std::endl;
    }

    if ( inputDB->keyExists( "b" ) ) {
        auto b = inputDB->getScalar<int>( "b" );
        AMP::pout << "b: " << b << std::endl;
    }

    if ( inputDB->keyExists( "myFirstDatabase" ) ) {
        auto myFirstDatabase = inputDB->getDatabase( "myFirstDatabase" );
        AMP::pout << "myFirstDatabase contains " << myFirstDatabase->getString( "greeting" )
                  << std::endl;
    }

    if ( inputDB->keyExists( "myArray" ) ) {
        auto myArray = inputDB->getVector<int>( "myArray" );
        AMP::pout << "myArray: ";
        for ( auto &v : myArray )
            AMP::pout << v << " ";
        AMP::pout << std::endl;
    }

    // Every AMP program needs to call shutdown
    // to properly release resources used by AMP
    // and third party libraries it is using
    AMP::AMPManager::shutdown();

    return 0;
}
