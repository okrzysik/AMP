#include "AMP/IO/PIO.h"
#include "AMP/operators/OperatorParameters.h"
#include "AMP/utils/AMPManager.h"
#include "AMP/utils/AMP_MPI.h"
#include "AMP/utils/Database.h"
#include "AMP/utils/UnitTest.h"
#include <memory>

#include <fstream>
#include <iostream>
#include <string>
#include <sys/stat.h>


class TestParameters : public AMP::Operator::OperatorParameters
{
public:
    explicit TestParameters( std::shared_ptr<AMP::Database> db )
        : AMP::Operator::OperatorParameters( db )
    {
    }
};


static void runTest( AMP::UnitTest *ut )
{
    std::string input_file = "inputOperatorParameters1";
    std::string log_file   = "outputOperatorParameters1";

    AMP::logOnlyNodeZero( log_file );

    auto input_db = AMP::Database::parseInputFile( input_file );
    input_db->print( AMP::plog );

    std::shared_ptr<AMP::Database> test_db = input_db->getDatabase( "Test" );

    if ( ( test_db.get() ) == NULL ) {
        ut->failure( "Testing getDatabase" );
    } else {
        ut->passes( "Testing getDatabase" );
    }

    if ( test_db->keyExists( "print_info_level" ) ) {
        ut->passes( "Testing keyExists" );
    } else {
        ut->failure( "Testing keyExists" );
    }

    auto discreteOpParams = std::make_shared<AMP::Operator::OperatorParameters>( test_db );

    if ( ( discreteOpParams.get() ) == NULL ) {
        ut->failure( "Testing OperatorParameters' Constructor" );
    } else {
        ut->passes( "Testing OperatorParameters' Constructor" );
    }

    if ( ( ( discreteOpParams->d_db ).get() ) != ( test_db.get() ) ) {
        ut->failure( "Testing OperatorParameters::d_db" );
    } else {
        ut->passes( "Testing OperatorParameters::d_db" );
    }

    if ( ( discreteOpParams->d_db )->keyExists( "print_info_level" ) ) {
        ut->passes( "Testing OperatorParameters::d_db keyExists" );
    } else {
        ut->failure( "Testing OperatorParameters::d_db keyExists" );
    }

    auto testParams = std::make_shared<TestParameters>( test_db );

    if ( ( testParams.get() ) == NULL ) {
        ut->failure( "Testing TestParameters' Constructor" );
    } else {
        ut->passes( "Testing TestParameters' Constructor" );
    }

    if ( ( ( testParams->d_db ).get() ) != ( test_db.get() ) ) {
        ut->failure( "Testing TestParameters::d_db" );
    } else {
        ut->passes( "Testing TestParameters::d_db" );
    }

    if ( ( testParams->d_db )->keyExists( "print_info_level" ) ) {
        ut->passes( "Testing TestParameters::d_db keyExists" );
    } else {
        ut->failure( "Testing TestParameters::d_db keyExists" );
    }

    std::shared_ptr<AMP::Operator::OperatorParameters> testParamCopy( testParams );

    if ( ( testParamCopy.get() ) == NULL ) {
        ut->failure( "Testing Copy-1" );
    } else {
        ut->passes( "Testing Copy-1" );
    }

    if ( ( ( testParamCopy->d_db ).get() ) != ( test_db.get() ) ) {
        ut->failure( "Testing Copy-2" );
    } else {
        ut->passes( "Testing Copy-2" );
    }

    if ( ( testParamCopy->d_db )->keyExists( "print_info_level" ) ) {
        ut->passes( "Testing Copy-3" );
    } else {
        ut->failure( "Testing Copy-3" );
    }

    auto secondTestParam =
        std::dynamic_pointer_cast<TestParameters, AMP::Operator::OperatorParameters>(
            testParamCopy );

    if ( ( secondTestParam.get() ) == NULL ) {
        ut->failure( "Testing Dynamic Cast-1" );
    } else {
        ut->passes( "Testing Dynamic Cast-1" );
    }

    if ( ( ( secondTestParam->d_db ).get() ) != ( test_db.get() ) ) {
        ut->failure( "Testing Dynamic Cast-2" );
    } else {
        ut->passes( "Testing Dynamic Cast-2" );
    }

    if ( ( secondTestParam->d_db )->keyExists( "print_info_level" ) ) {
        ut->passes( "Testing Dynamic Cast-3" );
    } else {
        ut->failure( "Testing Dynamic Cast-3" );
    }

    AMP::AMPManager::shutdown();
}

int testDiscreteOperatorParameters( int argc, char **argv )
{
    AMP::AMPManager::startup( argc, argv );
    AMP::UnitTest ut;

    runTest( &ut );

    ut.report();

    int num_failed = ut.NumFailGlobal();
    AMP::AMPManager::shutdown();
    return num_failed;
}
