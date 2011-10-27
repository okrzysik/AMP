#include "string.h"
#include "utils/AMPManager.h"
#include "utils/UnitTest.h"

#include "test_ContactSearchTests.h"

using namespace AMP::unit_test;


int main ( int argc , char **argv )
{

    AMP::AMPManager::startup(argc, argv);
    AMP::UnitTest ut;

    InstantiateCoarseSearch<>::run_test ( &ut );
    SimpleCoarseSearch<>::run_test ( &ut );
    FindNodesInElementsCoarseSearch<>::run_test ( &ut );
    CallContactSearchSerialPenaltyInterface<>::run_test ( &ut );
    CallContactSearchSerialPenaltyInterface<ExodusReaderGenerator>::run_test ( &ut );
    VerifyContactSearchSerialPenaltyInterface<ExodusReaderGenerator>::run_test ( &ut );

    ut.report ();

    int num_failed = ut.NumFailGlobal();
    AMP::AMPManager::shutdown();
    return num_failed;

}
