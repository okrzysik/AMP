#include "AMP/utils/AMPManager.h"
#include "AMP/vectors/MultiVariable.h"
#include "AMP/vectors/MultiVector.h"
#include "AMP/vectors/testHelpers/VectorTests.h"

#include "test_VectorHelpers.h"


int test_VectorClone( int argc, char **argv )
{
    AMP::AMPManager::startup( argc, argv );
    AMP::UnitTest ut;

    
#if defined( USE_EXT_TRILINOS ) &&  defined( USE_EXT_PETSC )
    AMP::pout << "Testing Iterator" << std::endl;
    VectorIteratorTests( ut, CloneMVFactory1 );
    AMP::pout << std::endl;

    AMP::pout << "Testing NativePetscVector" << std::endl;
    testManagedVector( ut, CloneSNPVFactory );
    AMP::pout << std::endl;

    AMP::pout << "Testing simple multivector" << std::endl;
    testManagedVector( ut, CloneMVFactory1 );
    AMP::pout << std::endl;

    AMP::pout << "Testing bigger multivector" << std::endl;
    testManagedVector( ut, CloneMVFactory2 );
    AMP::pout << std::endl;

    AMP::pout << "Testing multivector of multivector" << std::endl;
    testManagedVector( ut, CloneMVFactory3 );
    AMP::pout << std::endl;
#else
    ut.expected_failure( "Compiled without petsc or trilinos" );
#endif

    ut.report();

    int num_failed = ut.NumFailGlobal();
    AMP::AMPManager::shutdown();
    return num_failed;
}
