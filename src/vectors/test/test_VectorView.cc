#include "vectors/MultiVariable.h"
#include "vectors/MultiVector.h"

#include "utils/AMPManager.h"

#include "test_VectorHelpers.h"


int main( int argc, char **argv )
{
    AMP::AMPManager::startup( argc, argv );
    AMP::UnitTest ut;

#if defined( USE_EXT_PETSC ) && defined( USE_EXT_TRILINOS )
    AMP::pout << "Testing Iterator" << std::endl;
    VectorIteratorTests( ut, ViewMVFactory1 );
    AMP::pout << std::endl;

    AMP::pout << "Testing ManagedEpetraVector" << std::endl;
    VectorIteratorTests( ut, ViewSMEVFactory );
    AMP::pout << std::endl;

    AMP::pout << "Testing NativePetscVector" << std::endl;
    VectorIteratorTests( ut, ViewSNPVFactory );
    AMP::pout << std::endl;

    AMP::pout << "Testing simple multivector" << std::endl;
    VectorIteratorTests( ut, ViewMVFactory1 );
    AMP::pout << std::endl;

    AMP::pout << "Testing bigger multivector" << std::endl;
    VectorIteratorTests( ut, ViewMVFactory2 );
    AMP::pout << std::endl;

    AMP::pout << "Testing multivector of multivector" << std::endl;
    VectorIteratorTests( ut, ViewMVFactory3 );
    AMP::pout << std::endl;
#else
    ut.expected_failure( "Compiled without petsc or trilinos" );
#endif

    ut.report();

    int num_failed = ut.NumFailGlobal();
    AMP::AMPManager::shutdown();
    return num_failed;
}
