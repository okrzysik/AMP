#include "AMP/utils/AMPManager.h"
#include "AMP/utils/UnitTest.h"

#include "testHex8ElementLinearInterpolation.h"
#include "testHex8ElementMapping.h"
#include "testHex8ElementNormalToFaces.h"


int main( int argc, char *argv[] )
{
    AMP::AMPManager::startup( argc, argv );
    AMP::UnitTest ut;

    testHex8ElementLinearInterpolation( ut );
    testHex8ElementMapping( ut );
    testHex8ElementNormalToFaces( ut );

    ut.report();
    int num_failed = ut.NumFailGlobal();

    AMP::AMPManager::shutdown();
    return num_failed;
}
