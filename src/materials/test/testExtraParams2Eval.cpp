#include "AMP/materials/Material.h"
#include "AMP/utils/AMPManager.h"
#include "AMP/utils/UnitTest.h"
#include "AMP/utils/Utilities.h"

#include <algorithm>
#include <iostream>
#include <memory>
#include <string>
#include <vector>

void myTest( AMP::UnitTest *ut, const std::string &exeName )
{
    auto material = AMP::Materials::getMaterial( "UO2_MSRZC_09" );

    std::string temperatureString = "temperature";
    std::string burnupString      = "burnup";
    std::string oxygenString      = "concentration";

    auto tempVec   = std::make_shared<std::vector<double>>();
    auto burnupVec = std::make_shared<std::vector<double>>();
    auto oxygenVec = std::make_shared<std::vector<double>>();

    tempVec->push_back( 310.0 );
    burnupVec->push_back( 0.0 );
    oxygenVec->push_back( 0.0 );

    std::vector<double> YM( 1 );
    std::vector<double> PR( 1 );

    std::string ymString = "YoungsModulus";
    std::string prString = "PoissonRatio";

    for ( int i = 0; i < 2; i++ ) {
        if ( i == 0 ) {
            material->property( ymString )
                ->evalv( YM, {}, "temperature", tempVec, "concentration", oxygenVec );
            material->property( prString )
                ->evalv( PR, {}, "temperature", tempVec, "concentration", oxygenVec );
            std::cout << exeName << ": Passed if burnup is NOT used." << std::endl;
        } else {
            material->property( ymString )
                ->evalv( YM,
                         {},
                         "temperature",
                         tempVec,
                         "concentration",
                         oxygenVec,
                         "burnup",
                         burnupVec );
            material->property( prString )
                ->evalv( PR,
                         {},
                         "temperature",
                         tempVec,
                         "concentration",
                         oxygenVec,
                         "burnup",
                         burnupVec );
            std::cout << exeName << ": Passed if burnup is used." << std::endl;
        }
    } // end for i

    ut->passes( exeName );
}

int main( int argc, char *argv[] )
{
    AMP::AMPManager::startup( argc, argv );
    AMP::UnitTest ut;

    std::string exeName = "testExtraParams2Eval";

    myTest( &ut, exeName );

    ut.report();

    int num_failed = ut.NumFailGlobal();
    AMP::AMPManager::shutdown();
    return num_failed;
}
