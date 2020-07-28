#include "AMP/utils/AMPManager.h"
#include "AMP/vectors/MultiVariable.h"
#include "AMP/vectors/MultiVector.h"

#include "test_VectorHelpers.h"


int test_VectorView( int argc, char **argv )
{
    AMP::AMPManager::startup( argc, argv );
    AMP::UnitTest ut;

#if defined( USE_EXT_TRILINOS )
    AMP::pout << "Testing EpetraVector::view<SimpleVector>" << std::endl;
    VectorIteratorTests( ut, "ViewFactory<EpetraVector," + SimpleFactory1 + ">" );
    AMP::pout << std::endl;
#endif

#if defined( USE_EXT_PETSC )
    AMP::pout << "Testing PetscVector::view<SimpleVector>" << std::endl;
    VectorIteratorTests( ut, "ViewFactory<PetscVector," + SimpleFactory1 + ">" );
    AMP::pout << std::endl;
#endif

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
#endif

    ut.report();

    int num_failed = ut.NumFailGlobal();
    AMP::AMPManager::shutdown();
    return num_failed;
}
