#include "AMP/time_integrators/TimeIntegratorParameters.h"
#include "AMP/utils/AMPManager.h"
#include "AMP/utils/InputDatabase.h"
#include "AMP/utils/UnitTest.h"
#include "AMP/utils/Utilities.h"


class InstantiateTimeIntegratorParameter
{
public:
    static const char *get_test_name() { return "instantiate TimeIntegratorParameters"; }

    template<typename UTILS>
    static void run_test( UTILS *utils )
    {
        AMP::shared_ptr<AMP::InputDatabase> new_db( new AMP::InputDatabase( "Dummy db" ) );
        AMP::TimeIntegrator::TimeIntegratorParameters params( new_db );
        utils->passes( "instantiate TimeIntegratorParameters" );
    }
};


int testTimeIntegratorParameters( int argc, char **argv )
{
    AMP::AMPManager::startup( argc, argv );
    AMP::UnitTest ut;

    InstantiateTimeIntegratorParameter::run_test( &ut );

    ut.report();

    int num_failed = ut.NumFailGlobal();
    AMP::AMPManager::shutdown();
    return num_failed;
}
