#include "AMP/utils/AMPManager.h"
#include "AMP/utils/Database.h"
#include "AMP/utils/UnitTest.h"

#ifdef AMP_USE_SAMRAI
    #include "SAMRAI/tbox/InputManager.h"
    #include "SAMRAI/tbox/MemoryDatabase.h"
#endif

/************************************************************************
 * Main                                                                  *
 ************************************************************************/
int main( int argc, char *argv[] )
{
    AMP::AMPManager::startup( argc, argv );
    AMP::UnitTest ut;

#ifdef AMP_USE_SAMRAI
    auto amp_db = std::make_shared<AMP::Database>( "AMPDB" );
    std::vector<double> v{ 1.0, -1.0 };
    amp_db->putVector<double>( "v", v );
    auto samrai_db = amp_db->cloneToSAMRAI();
    auto vc        = samrai_db->getDoubleVector( "v" );
    if ( v == vc ) {
        ut.passes( "No conversion to ints" );
    } else {
        ut.passes( "Conversion from double to int happening" );
    }
#endif
    ut.report();

    int num_failed = ut.NumFailGlobal();
    AMP::AMPManager::shutdown();
    return num_failed;
}
